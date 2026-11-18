"""Independent orbit expansion checks of continuum residual matrices."""
from fractions import Fraction
import unittest

import numpy as np
from flint import acb, arb, ctx

from purposeful.validated_linear import BallArray as BA, quadratic
from purposeful.validated_residual import ResidualMatrices, ContinuumValidator, PROFILES
from purposeful.validated_spatial import SpatialTerm


def enclosing_arb(array):
    return arb(float(array.center),float(array.radius))


class ValidatedResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=ResidualMatrices(2,2)

    def basis(self,index,source,frequency_shift=(0,0,0),scale=1):
        name,derivative=source
        terms=[];factor=arb(scale)*arb(2)**(-int(self.m.nz[index]))
        for vector in self.m.orbits[index]:
            k=tuple(int(x) for x in vector)
            frequency=tuple(x-y for x,y in zip(k,frequency_shift))
            if not derivative:
                terms.append(SpatialTerm(PROFILES[name],frequency,coefficient=acb(factor)))
            else:
                for j in range(3):
                    if k[j]:
                        terms.append(SpatialTerm(PROFILES[name],frequency,(j,),
                            acb(0,2*arb.pi()*k[j]*factor)))
        return terms

    def test_transform_matches_direct_exponential_orbits(self):
        with ctx.workprec(120):
            for source in (('one',0),('w',0),('on',0),('drift',1)):
                matrix=self.m.transform(source)
                for output in (1,2,5):
                    shift=tuple(int(x) for x in self.m.outmodes[output])
                    for input_index in (1,2,5):
                        direct=self.m.oracle.integral(self.basis(input_index,source,shift))
                        self.assertTrue(direct.imag.contains(0))
                        self.assertTrue(enclosing_arb(matrix[output,input_index]).contains(direct.real),
                            str((source,output,input_index)))

    def test_gram_matches_direct_cartesian_expansion(self):
        with ctx.workprec(120):
            for left,right in ((('on',0),('w',0)),(('w',0),('drift',1)),
                               (('drift',1),('drift',1))):
                matrix=self.m.gram(left,right)
                for i,j in ((1,2),(2,1),(5,5)):
                    direct=self.m.oracle.inner(self.basis(i,left),self.basis(j,right))
                    self.assertTrue(direct.imag.contains(0))
                    self.assertTrue(enclosing_arb(matrix[i,j]).contains(direct.real),
                        str((left,right,i,j)))

    def test_random_dyadic_gram_contraction(self):
        a=np.zeros(len(self.m.modes)); b=a.copy()
        a[[0,1,5]]=[.5,-.75,.625];b[[0,1,5]]=[-.25,.875,.125]
        with ctx.workprec(120):
            left=('on',0);right=('drift',1)
            direct_left=[];direct_right=[]
            for index in (0,1,5):
                direct_left.extend(self.basis(index,left,scale=float(a[index])))
                direct_right.extend(self.basis(index,right,scale=float(b[index])))
            direct=self.m.oracle.inner(direct_left,direct_right)
            contraction=quadratic(BA.exact(a),self.m.gram(left,right),BA.exact(b))
            self.assertTrue(enclosing_arb(contraction).contains(direct.real))

    def test_chemical_ledger_constant_trial(self):
        v=ContinuumValidator(self.m)
        coeff=np.zeros((20,len(self.m.modes)));coeff[:,0]=1
        P,Q,S,A=v.sources(BA.exact(coeff),BA.exact(np.zeros_like(coeff)),"spectator")
        # Constant one has no repair or organization difference. Killing arises
        # only from the h=1 damage exit, at rate1; all live health jumps cancel.
        self.assertTrue(np.all(Q.lower()<=0) and np.all(Q.upper()>=0))
        self.assertTrue(np.all(S.lower()<=0) and np.all(S.upper()>=0))
        for i,(h,b,q) in enumerate(v.states):
            expected=np.zeros(len(self.m.modes));expected[0]=int(h==1)
            self.assertTrue(np.all(P[i].lower()<=expected))
            self.assertTrue(np.all(P[i].upper()>=expected))
        # Independent exact equilibrium mass, with healthy ratio c1/c2=4
        # and fuel ratio c(b+1)/c(b)=(4-b)/(8(b+1)).
        for h,b,q in v.states:
            i=v.ids[(h,b,q)]
            if h==1:
                self.assertEqual(v.chemical[i],4*v.chemical[v.ids[(2,b,q)]])
            if b<4:
                self.assertEqual(v.chemical[v.ids[(h,b+1,q)]]/v.chemical[i],
                                 Fraction(4-b,8*(b+1)))

    def test_sources_against_independent_reaction_edge_ledger(self):
        v=ContinuumValidator(self.m)
        rng=np.random.default_rng(4811)
        exact=[[Fraction(int(n),8) for n in row] for row in
               rng.integers(-3,4,size=(20,len(self.m.modes))) ]
        coeff=np.array(exact,dtype=float)
        for policy in ("spectator",Fraction(3,7)):
            P,Q,S,A=v.sources(BA.exact(coeff),BA.exact(np.zeros_like(coeff)),policy)
            qref=[[Fraction(0) for _ in row] for row in exact]
            sref=[[Fraction(0) for _ in row] for row in exact]
            href=[[Fraction(0) for _ in row] for row in exact]
            def edge(target,origin,destination,rate):
                i=v.ids[origin];j=v.ids[destination] if destination else None
                for mode in range(len(self.m.modes)):
                    target[i][mode]+=rate*(exact[i][mode]-(exact[j][mode] if j is not None else 0))
            for b in range(5):
                for q in range(2):
                    edge(href,(1,b,q),None,Fraction(1))
                    edge(href,(1,b,q),(2,b,q),Fraction(1,2))
                    edge(href,(2,b,q),(1,b,q),Fraction(2))
            for b in range(1,5):
                for q in range(2):
                    edge(qref,(1,b,q),(2,b-1,q),Fraction(16*b))
                    edge(qref,(2,b-1,q),(1,b,q),Fraction(8*(5-b)))
            for h in (1,2):
                for b in range(5):
                    edge(sref,(h,b,0),(h,b,1),Fraction(1))
                    edge(sref,(h,b,1),(h,b,0),Fraction(1))
                for b in range(1,5):
                    conduct=Fraction(b-1,3) if policy=="spectator" else policy
                    edge(sref,(h,b,0),(h,b-1,1),4*b*conduct)
                    edge(sref,(h,b-1,1),(h,b,0),Fraction(5-b,2)*conduct)
            with ctx.workprec(120):
                for i in range(20):
                    for j,k in enumerate(self.m.modes):
                        def ab(q):return arb(q.numerator)/q.denominator
                        expected=ab(href[i][j])+ab(exact[i][j])*sum(int(n*n) for n in k)*arb.pi()**2/10
                        self.assertTrue(enclosing_arb(P[i,j]).contains(expected))
                        self.assertTrue(enclosing_arb(Q[i,j]).contains(ab(qref[i][j])))
                        self.assertTrue(enclosing_arb(S[i,j]).contains(ab(sref[i][j])))


if __name__=="__main__":unittest.main()
