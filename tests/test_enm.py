import unittest
import numpy as np
from scipy.linalg import pinvh
from enm import kirchhoff, hessian, modes_and_covariance, fluctuations, contacts


class ScientificChecks(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[1.,1.,1.]])

    def test_gnm_laplacian_and_single_zero(self):
        k = kirchhoff(self.x, 3)
        np.testing.assert_allclose(k.sum(axis=1), 0)
        self.assertEqual(modes_and_covariance(k)[3], 1)

    def test_anm_six_rigid_body_zeros(self):
        h = hessian(self.x, 3)
        values, _, _, zeros = modes_and_covariance(h)
        self.assertEqual(zeros, 6)
        self.assertGreaterEqual(values.min(), -1e-8)

    def test_covariances_against_independent_pinv(self):
        for matrix in (kirchhoff(self.x, 3), hessian(self.x, 3)):
            np.testing.assert_allclose(modes_and_covariance(matrix)[2], pinvh(matrix, atol=1e-8), atol=1e-12)

    def test_hessian_energy_finite_difference(self):
        h = hessian(self.x, 3)
        displacement = np.random.default_rng(42).normal(size=self.x.shape)
        _, _, d2, adjacency = contacts(self.x, 3)
        pairs = np.triu(adjacency, 1)
        def energy(y):
            lengths = np.sqrt(np.sum((y[:, None]-y[None, :])**2, axis=2))
            return 0.5*np.sum((lengths[pairs]-np.sqrt(d2[pairs]))**2)
        epsilon = 1e-4
        curvature = (energy(self.x+epsilon*displacement)+energy(self.x-epsilon*displacement))/epsilon**2
        expected = displacement.ravel() @ h @ displacement.ravel()
        self.assertAlmostEqual(curvature, expected, places=5)

    def test_rigid_translation_invariance(self):
        for builder in (kirchhoff, hessian):
            np.testing.assert_allclose(builder(self.x, 3), builder(self.x+[20,8,-5], 3), atol=1e-12)

    def test_rotation_covariance(self):
        q, _ = np.linalg.qr(np.random.default_rng(17).normal(size=(3,3)))
        transform = np.kron(np.eye(len(self.x)), q.T)
        original = modes_and_covariance(hessian(self.x,3))[2]
        rotated = modes_and_covariance(hessian(self.x@q,3))[2]
        np.testing.assert_allclose(rotated, transform@original@transform.T, atol=1e-12)
        np.testing.assert_allclose(fluctuations(original,True), fluctuations(rotated,True), atol=1e-12)

    def test_disconnected_network_zeros_are_visible(self):
        self.assertGreater(modes_and_covariance(kirchhoff(self.x, .1))[3], 1)

    def test_invalid_and_coincident_coordinates(self):
        for x, cutoff in [(self.x,0),(self.x,float('nan')),([[0,0,0],[0,0,0]],3),([[0,0,0],[1,float('nan'),0]],3)]:
            with self.assertRaises(ValueError):
                hessian(x, cutoff)


if __name__ == '__main__':
    unittest.main()
