"""Validated continuum residuals in the raw symmetric cosine basis.

Arb encloses radial coefficients. Explicit midpoint/radius BLAS bounds enclose
finite contractions. Fourier Parseval bounds all omitted frequencies; no grid
quadrature or measured convergence rate enters an error radius.
"""
from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product
from pathlib import Path
import json
from hashlib import sha256
import numpy as np
from flint import arb, ctx

from .validated_linear import BallArray as BA, rational_ball, quadratic
from .validated_spatial import (SpatialOracle, ONE, CONTACT, WEIGHT, SWITCH,
                               DRIFT, midpoint_radius)

PROFILES = {'one':ONE, 'g':CONTACT, 'w':WEIGHT, 'wg':WEIGHT*CONTACT,
            'on':SWITCH, 'drift':DRIFT}
SOURCE_NAMES = [('one',0), ('g',0), ('on',0), ('w',0), ('wg',0), ('drift',1)]


def orbit(k):
    return np.array(sorted(set(tuple(p[i]*s[i] for i in range(3))
        for p in permutations(k) for s in product((-1,1), repeat=3))), dtype=int)


class ResidualMatrices:
    def __init__(self, cutoff=12, flux_cutoff=24, cache=None):
        self.K, self.M = cutoff, flux_cutoff
        self.modes = np.array(list(combinations_with_replacement(range(cutoff+1),3)))
        self.outmodes = np.array(list(combinations_with_replacement(range(flux_cutoff+1),3)))
        self.orbits = [orbit(k) for k in self.modes]
        self.sizes = np.array([len(o) for o in self.orbits])
        self.nz = np.count_nonzero(self.modes, axis=1)
        self.outsizes = np.array([len(orbit(k)) for k in self.outmodes])
        self.oracle = SpatialOracle(bits=100,target_bits=70)
        self.cache = Path(cache) if cache else None
        if self.cache:
            self.cache.mkdir(parents=True,exist_ok=True)
            # Reusing a cache produced by another formula/version is forbidden.
            digest=sha256(Path(__file__).read_bytes()+Path(__file__).with_name('validated_linear.py').read_bytes()
                          +Path(__file__).with_name('validated_spatial.py').read_bytes()).hexdigest()
            identity=self.cache/'SOURCE-SHA256.txt'
            if identity.exists() and identity.read_text().strip()!=digest:
                raise ArithmeticError('matrix cache source mismatch; choose a new cache directory')
            identity.write_text(digest+'\n')
        self._tables, self._matrices = {}, {}

    def table(self, p, max_n, d=0):
        key = (p,max_n,d)
        if key not in self._tables:
            self._tables[key] = BA(*self.oracle.table(p,max_n,d))
        return self._tables[key]

    def _cached(self, name, generate):
        if name in self._matrices:
            return self._matrices[name]
        target = self.cache/f'{name}-K{self.K}-M{self.M}.npz' if self.cache else None
        if target and target.exists():
            checksum=target.with_suffix('.sha256')
            if not checksum.exists() or checksum.read_text().strip()!=sha256(target.read_bytes()).hexdigest():
                raise ArithmeticError('matrix cache checksum mismatch')
            data=np.load(target)
            ans=BA(data['center'],data['radius'])
        else:
            ans=generate()
            if target:
                np.savez_compressed(target,center=ans.center,radius=ans.radius)
                target.with_suffix('.sha256').write_text(sha256(target.read_bytes()).hexdigest()+'\n')
        self._matrices[name]=ans
        return ans

    def gram(self, left, right):
        """Integral radial profiles times phi or x.grad(phi) on both sides."""
        nl,dl=left; nr,dr=right
        key=f'gram-{nl}{dl}-{nr}{dr}'
        def make():
            p=PROFILES[nl]*PROFILES[nr]
            nmax=12*self.K*self.K
            tables={d:self.table(p,nmax,d) for d in range(dl+dr+1)}
            c=np.empty((len(self.modes),len(self.modes))); r=np.empty_like(c)
            k=self.modes[:,None,:]
            for j,l0 in enumerate(self.orbits):
                l=l0[None,:,:]; v=k+l; n2=np.sum(v*v,axis=-1)
                if dl+dr==0:
                    values=tables[0][n2]
                elif dl+dr==1:
                    factors=2*np.sum((k if dl else l)*v,axis=-1)
                    values=tables[1][n2]*factors
                else:
                    f1=2*np.sum(k*l,axis=-1)
                    f2=4*np.sum(k*v,axis=-1)*np.sum(l*v,axis=-1)
                    values=tables[1][n2]*f1+tables[2][n2]*f2
                ans=values.sum(axis=1)*(self.sizes*2.**(-self.nz-self.nz[j]))
                c[:,j],r[:,j]=ans.center,ans.radius
            return BA(c,r)
        return self._cached(key,make)

    def transform(self, source):
        """Each row is one Fourier coefficient of profile*phi or profile*x.gradphi."""
        name,deriv=source
        def make():
            p=PROFILES[name]; nmax=3*(self.K+self.M)**2
            tab=self.table(p,nmax,deriv)
            c=np.empty((len(self.outmodes),len(self.modes))); r=np.empty_like(c)
            n=self.outmodes[:,None,:]
            for j,k0 in enumerate(self.orbits):
                k=k0[None,:,:]; v=n-k; n2=np.sum(v*v,axis=-1)
                values=tab[n2]
                if deriv:
                    values=values*(-2*np.sum(k*v,axis=-1))
                ans=values.sum(axis=1)*2.**(-self.nz[j])
                c[:,j],r[:,j]=ans.center,ans.radius
            return BA(c,r)
        return self._cached(f'transform-{name}{deriv}',make)


def power_to_legendre(coefficients):
    """Exact rational transformation enclosed in binary64 balls.

    Integral_0^1 s^k L_n(2s-1) ds=(k!)²/((k-n)!(k+n+1)!).
    The resulting component has time norm weight 1/(2*n+1).
    """
    from math import factorial
    if not isinstance(coefficients, BA):
        coefficients=BA.exact(coefficients)
    degree=coefficients.center.shape[0]-1
    rows=[]
    for n in range(degree+1):
        val=BA.exact(np.zeros(coefficients.center.shape[1:]))
        for k in range(n,degree+1):
            q=Fraction((2*n+1)*factorial(k)**2,factorial(k-n)*factorial(k+n+1))
            val=val+coefficients[k]*rational_ball(q)
        rows.append(val)
    return BA(np.stack([v.center for v in rows]),np.stack([v.radius for v in rows]))


class ContinuumValidator:
    def __init__(self, matrices):
        self.m=matrices
        self.states=[(h,b,q) for h in (1,2) for b in range(5) for q in (0,1)]
        self.ids={s:i for i,s in enumerate(self.states)}
        from math import comb
        self.chemical=[Fraction(comb(2,h)*2**(2-h),9)*Fraction(comb(4,b)*8**(4-b),9**4)
                       for h,b,q in self.states]
        with ctx.workprec(120):
            self.pi2=BA(*midpoint_radius((2*arb.pi())**2))
            self.inv_pi2=BA(*midpoint_radius(1/(2*arb.pi())**2))
            self.exp2=BA(*midpoint_radius(arb(2).exp()))
            self.inv_lambda=BA(*midpoint_radius(4/(7-arb(17).sqrt())))
        self.lap=self.pi2*BA.exact(np.sum(self.m.modes**2,axis=1))*rational_ball(Fraction(1,40))
        self.drift=rational_ball(Fraction(4,1)*Fraction(1,40)*6/Fraction(8,25)**2)

    def sources(self, coeff, deriv, policy):
        """R_i/rho_chemical, before radial multipliers; accepts batched times."""
        shape=coeff.center.shape
        P=deriv+coeff*self.lap
        pc=P.center.copy(); pr=P.radius.copy()
        Q=BA.exact(np.zeros(shape)); qc=Q.center.copy(); qr=Q.radius.copy()
        S=BA.exact(np.zeros(shape)); sc=S.center.copy(); sr=S.radius.copy()
        for i,(h,b,q) in enumerate(self.states):
            # Killed thermal health generator, including exits to h=0.
            health=coeff[...,i,:]*rational_ball(Fraction(3,2) if h==1 else 2)
            health=health-coeff[...,self.ids[(3-h,b,q)],:]*rational_ball(Fraction(1,2) if h==1 else 2)
            val=P[...,i,:]+health
            pc[...,i,:],pr[...,i,:]=val.center,val.radius
            repair=BA.exact(np.zeros(shape[:-2]+(shape[-1],)))
            if h==1 and b:
                repair=(coeff[...,i,:]-coeff[...,self.ids[(2,b-1,q)],:])*(16*b)
            if h==2 and b<4:
                repair=(coeff[...,i,:]-coeff[...,self.ids[(1,b+1,q)],:])*(8*(4-b))
            qc[...,i,:],qr[...,i,:]=repair.center,repair.radius
            switch=coeff[...,i,:]-coeff[...,self.ids[(h,b,1-q)],:]
            if q==0 and b:
                conduct=Fraction(b-1,3) if policy=='spectator' else Fraction(policy)
                switch=switch+(coeff[...,i,:]-coeff[...,self.ids[(h,b-1,1)],:])*rational_ball(4*b*conduct)
            if q==1 and b<4:
                conduct=Fraction(b,3) if policy=='spectator' else Fraction(policy)
                switch=switch+(coeff[...,i,:]-coeff[...,self.ids[(h,b+1,0)],:])*rational_ball(Fraction(4-b,2)*conduct)
            sc[...,i,:],sr[...,i,:]=switch.center,switch.radius
        return BA(pc,pr), BA(qc,qr), BA(sc,sr), coeff*self.drift

    def source_list(self, sources, i):
        P,Q,S,A=sources
        if self.states[i][2]==0:
            return [(('one',0),P[...,i,:]),(('g',0),Q[...,i,:]),(('on',0),S[...,i,:])]
        return [(('w',0),P[...,i,:]),(('wg',0),Q[...,i,:]),(('on',0),S[...,i,:]),(('drift',1),A[...,i,:])]

    def residual_energy(self, coeff, deriv, policy):
        """For one spatial coefficient vector, returns dual-norm-square upper.

        The bound uses sup reciprocal spatial density (1 or e²), unweighted
        resolved Fourier flux, and a Parseval L2 tail beyond the complete cube.
        """
        sources=self.sources(coeff,deriv,policy)
        batch=coeff.center.shape[:-2]
        total=BA.exact(np.zeros(batch)); diagnostics={'tail_upper':np.zeros(batch),'resolved_upper':np.zeros(batch),'mean_upper':np.zeros(batch)}
        freq2=np.sum(self.m.outmodes**2,axis=1)
        for i,state in enumerate(self.states):
            terms=self.source_list(sources,i)
            rhat=BA.exact(np.zeros(batch+(len(self.m.outmodes),)))
            strong=BA.exact(np.zeros(batch))
            for a,(sa,va) in enumerate(terms):
                rhat=rhat+va.matmul(self.m.transform(sa).T)
                for b in range(a,len(terms)):
                    sb,vb=terms[b]
                    strong=strong+quadratic(va,self.m.gram(sa,sb),vb)*(1 if a==b else 2)
            squares=rhat.squared()
            # Sum of true squares has nonnegative lower endpoint.
            lower=np.maximum(squares.lower(),0.)
            resolved=(BA.exact(self.m.outsizes)*BA(lower,0.)).sum(axis=-1)
            # Tail upper uses total upper minus a certified LOWER resolved norm.
            tail=strong-BA.exact(resolved.lower())
            if np.any(tail.upper() < 0):
                raise ArithmeticError('negative Parseval tail enclosure')
            tailupper=tail.upper()
            if np.any(strong.upper() < -1.e-30):
                raise ArithmeticError('negative strong norm enclosure')
            tailterm=BA.exact(tailupper)*rational_ball(Fraction(40,(self.m.M+1)**2))
            # pi² divisor common to resolved and tail terms.
            invfreq=rational_ball(Fraction(40))*BA.exact(self.m.outsizes[1:]/freq2[1:])
            # Rounded outsizes/freq2 is NOT an exact rational: add its division budget.
            invfreq=invfreq+BA(np.zeros(len(freq2)-1),np.abs(invfreq.center)*8*np.finfo(float).eps)
            flux=(squares[...,1:]*invfreq).sum(axis=-1)+tailterm
            invpi=self.inv_pi2
            flux=flux*invpi
            mean=squares[...,0]*self.inv_lambda
            pref=rational_ball(self.chemical[i])*(self.exp2 if state[2] else 1)
            total=total+(flux+mean)*pref
            diagnostics['tail_upper']+=(tailterm*invpi*pref).upper()
            diagnostics['resolved_upper']+=((flux-tailterm*invpi)*pref).upper()
            diagnostics['mean_upper']+=(mean*pref).upper()
        return total.upper(), diagnostics

    def pairing(self, left, residual_sources):
        answer=BA.exact(0.)
        for i in range(20):
            val=BA.exact(0.)
            for source,vec in self.source_list(residual_sources,i):
                val=val+quadratic(left[...,i,:],self.m.gram(('one',0),source),vec)
            answer=answer+val*rational_ball(self.chemical[i])
        return answer

    def observable(self, terminal):
        answer=BA.exact(0.)
        for i,(_,_,q) in enumerate(self.states):
            row=self.m.gram(('one',0),('w' if q else 'one',0))[0]
            answer=answer+row.matmul(terminal[i])*rational_ball(self.chemical[i])
        return answer
