import unittest
import numpy as np
from enm import hessian, modes_and_covariance, scale_b_factors, fluctuations
from analysis_methods import (chain_matrix, cross_correlation, principal_frame,
                              directional_variances, equilibrium_matrix)


class AdditionalScientificChecks(unittest.TestCase):
    def test_free_chain_analytic_covariance(self):
        n = 41
        k = chain_matrix(n)
        _, _, c, zeros = modes_and_covariance(k)
        # Exact diagonal of the Moore-Penrose inverse of a path Laplacian.
        i = np.arange(n)
        expected = (i-(n-1)/2)**2/n + (n*n-1)/(12*n)
        np.testing.assert_allclose(np.diag(c), expected, atol=1e-10)
        self.assertEqual(zeros, 1)
        np.testing.assert_allclose(c@np.ones(n), 0, atol=1e-9)

    def test_boundary_interpretations_are_distinct(self):
        free, literal, ring = [modes_and_covariance(chain_matrix(41,b))[2]
                               for b in ('free','literal','ring')]
        self.assertGreater(free[0,0], free[20,20])
        self.assertLess(literal[0,0], literal[20,20])
        np.testing.assert_allclose(np.diag(ring), ring[0,0], atol=1e-10)

    def test_equilibrium_operator_independent_hessian(self):
        xyz = np.random.default_rng(51).normal(size=(8,3))
        b, _ = equilibrium_matrix(xyz, 10)
        np.testing.assert_allclose(b@b.T, hessian(xyz,10), atol=1e-12)

    def test_principal_axes_and_directional_trace(self):
        xyz = np.random.default_rng(4).normal(size=(9,3))
        center, axes, _ = principal_frame(xyz)
        np.testing.assert_allclose(axes.T@axes, np.eye(3), atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(axes), 1)
        q = (xyz-center)@axes
        np.testing.assert_allclose(q.T@q, np.diag(np.diag(q.T@q)), atol=1e-12)
        _, _, c, _ = modes_and_covariance(hessian(xyz,10))
        v = directional_variances(c, axes)
        np.testing.assert_allclose(v.sum(axis=1), fluctuations(c,True), atol=1e-12)
        np.testing.assert_allclose((v-v.mean(axis=1)[:,None]).sum(axis=1),0,atol=1e-12)

    def test_correlation_normalization_and_thermal_calibration(self):
        rng = np.random.default_rng(22)
        a = rng.normal(size=(15,15)); c = a@a.T
        for anisotropic in (True,False):
            cc = cross_correlation(c,anisotropic)
            np.testing.assert_allclose(np.diag(cc),1)
            self.assertLessEqual(np.max(np.abs(cc)),1+1e-12)
            np.testing.assert_allclose(cc,cc.T)
        f = fluctuations(c,True); measured=np.arange(1,6.)
        predicted, multiplier=scale_b_factors(f,measured)
        self.assertAlmostEqual(predicted.sum(),measured.sum())
        np.testing.assert_allclose(predicted/(8*np.pi**2/3), f*multiplier)


if __name__ == '__main__':
    unittest.main()
