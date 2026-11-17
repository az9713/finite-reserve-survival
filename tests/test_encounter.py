"""Check dimensional scaling and analytical diffusion formulas independently."""
from math import pi
import unittest

import numpy as np
from scipy.integrate import quad, solve_bvp

from purposeful.encounter import (shell_radius, mean_capture_time, capture_before_closure,
                                 uniform_capture_probability, partial_reactivity_mean)


class EncounterTests(unittest.TestCase):
    def test_mean_from_backward_equation_and_scaling(self):
        volume, diffusion, a = 1., 7.7, .005
        outer = shell_radius(volume, a)
        self.assertAlmostEqual(4*pi*(outer**3-a**3)/3,volume)
        mean_at = lambda r: (outer**3*(1/a-1/r)-(r*r-a*a)/2)/(3*diffusion)
        integrated = quad(lambda r: 4*pi*r*r*mean_at(r), a, outer)[0]/volume
        self.assertAlmostEqual(mean_capture_time(volume,diffusion,a),integrated,places=12)
        self.assertAlmostEqual(mean_capture_time(volume,2*diffusion,a),integrated/2)
        self.assertAlmostEqual(mean_capture_time(8*volume,diffusion,2*a),4*integrated)
        self.assertAlmostEqual(partial_reactivity_mean(volume,diffusion,a,10)-integrated,
                               volume/(4*pi*a*a*10))
        kappa=10.
        mesh=np.geomspace(a,outer,150)
        solution=solve_bvp(lambda r,y: np.array([y[1],-1/diffusion-2*y[1]/r]),
                           lambda ya,yb: np.array([diffusion*ya[1]-kappa*ya[0],yb[1]]),
                           mesh,np.array([np.full_like(mesh,300.),np.zeros_like(mesh)]),
                           tol=1e-9,max_nodes=15000)
        self.assertTrue(solution.success)
        robin_average=quad(lambda r:4*pi*r*r*solution.sol(r)[0],a,outer)[0]/volume
        self.assertAlmostEqual(partial_reactivity_mean(volume,diffusion,a,kappa),robin_average,delta=2e-8)

    def test_probability_against_independent_boundary_value_solver(self):
        volume,diffusion,a = 1.,7.7,.005
        outer=shell_radius(volume,a)
        mesh=np.geomspace(a,outer,150)
        for gamma in [.01,1.,100.]:
            solution=solve_bvp(lambda r,y: np.array([y[1],gamma/diffusion*y[0]-2*y[1]/r]),
                               lambda ya,yb: np.array([ya[0]-1,yb[1]]),
                               mesh,np.array([a/mesh,-a/mesh**2]),tol=1e-9,max_nodes=15000)
            self.assertTrue(solution.success)
            for r in np.linspace(a,outer,15):
                self.assertAlmostEqual(capture_before_closure(r,volume,diffusion,a,gamma),
                                       solution.sol(r)[0],delta=2e-8)
            average=quad(lambda r:4*pi*r*r*solution.sol(r)[0],a,outer,epsabs=1e-10)[0]/volume
            self.assertAlmostEqual(uniform_capture_probability(volume,diffusion,a,gamma),average,delta=2e-8)

    def test_small_rate_and_boundary_limits(self):
        volume,diffusion,a = 1.,7.7,.005
        mean=mean_capture_time(volume,diffusion,a)
        self.assertAlmostEqual(uniform_capture_probability(volume,diffusion,a,0),1.)
        gamma=1e-7
        probability=uniform_capture_probability(volume,diffusion,a,gamma)
        self.assertAlmostEqual((1-probability)/gamma,mean,delta=2e-6)
        for gamma in [.01,1,100]:
            self.assertEqual(capture_before_closure(a,volume,diffusion,a,gamma),1.)
            self.assertGreater(uniform_capture_probability(volume,diffusion,a,gamma),0.)
            self.assertLess(uniform_capture_probability(volume,diffusion,a,gamma),1.)

    def test_invalid_inputs(self):
        for volume,diffusion,a in [(0,1,1),(1,0,1),(1,1,-1)]:
            with self.assertRaises(ValueError):mean_capture_time(volume,diffusion,a)
        with self.assertRaises(ValueError):uniform_capture_probability(1,1,.01,-1)
        with self.assertRaises(ValueError):capture_before_closure(2,1,1,.01,1)


if __name__=='__main__':
    unittest.main()
