"""Numerical and scientific checks for the isolated branch-variability study."""

import unittest

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.stats import multivariate_normal

from full_tree_pareto import publication_analysis as analysis
from full_tree_pareto.branch_variability_sensitivity import (
    GAUSSIAN, MIXTURE, fit_model, gaussian_mean_norm, paired_norms,
    radial_coordinates, radial_logpdf, second_moment, subtree_folds,
)
from full_tree_pareto.clock_sensitivity import edge_manifest


class BranchVariabilityTests(unittest.TestCase):
    def test_zero_variability_is_gaussian(self):
        cov = np.array([[2., .3], [.3, 1.]])
        x = np.array([[1., -2.], [0., 0.], [-.4, .8]])
        q, logdet = radial_coordinates(x, cov)
        np.testing.assert_allclose(radial_logpdf(q, 2, logdet, 0),
                                   multivariate_normal.logpdf(x, cov=cov))
        sim = paired_norms(cov, 0, 20, 100, 3)
        np.testing.assert_array_equal(sim[GAUSSIAN], sim[MIXTURE])
        other = paired_norms(cov, 1.2, 20, 100, 3)
        np.testing.assert_array_equal(sim[GAUSSIAN], other[GAUSSIAN])

    def test_quadrature_against_adaptive_integral(self):
        tau, dim = 1.1, 20
        for q in [.05, 5., 20., 200.]:
            def logkernel(z):
                y = tau*z - .5*tau*tau
                return -.5*z*z - .5*dim*y - .5*q*np.exp(-y) - .5*np.log(2*np.pi)
            mode = minimize_scalar(lambda z: -logkernel(z), bounds=(-15, 15), method="bounded").x
            shift = logkernel(mode)
            # Rescale before integration: absolute tolerances are otherwise
            # inappropriate for tiny tail densities.
            integral = quad(lambda z: np.exp(logkernel(z)-shift), -15, 15,
                            epsabs=1e-12, epsrel=1e-10)[0]
            expected = -.5*dim*np.log(2*np.pi) + np.log(integral) + shift
            self.assertAlmostEqual(radial_logpdf(np.array([q]), dim, 0, tau)[0], expected, places=7)

    def test_simulated_moments_and_determinism(self):
        cov = np.diag([1., 2., 4.])
        tau = .9
        draws = paired_norms(cov, tau, 400, 1_000, 47)
        again = paired_norms(cov, tau, 400, 1_000, 47)
        np.testing.assert_array_equal(draws[MIXTURE], again[MIXTURE])
        values = draws[MIXTURE]
        expected_sq = np.trace(cov)
        variance_sq = np.exp(tau*tau)*(np.trace(cov)**2+2*np.trace(cov@cov))-expected_sq**2
        self.assertLess(abs((values**2).mean()-expected_sq), 6*np.sqrt(variance_sq/values.size))
        expected_mean = gaussian_mean_norm(cov)*np.exp(-tau*tau/8)
        self.assertLess(abs(values.mean()-expected_mean), 6*np.sqrt((expected_sq-expected_mean**2)/values.size))

    def test_parameter_recovery_and_scale_equivariance(self):
        rng = np.random.default_rng(872)
        tau = .8
        delta = rng.normal(size=(3_000, 5))*np.exp(.5*(tau*rng.normal(size=3_000)-.5*tau*tau))[:, None]
        covariance, fitted, info = fit_model(delta)
        self.assertLess(abs(fitted-tau), .08)
        np.testing.assert_allclose(covariance, delta.T@delta/len(delta))
        scaled_cov, scaled_tau, _ = fit_model(delta*3)
        np.testing.assert_allclose(scaled_cov, covariance*9)
        self.assertAlmostEqual(fitted, scaled_tau, places=5)
        self.assertGreater(info["mixture_loglik"], info["gaussian_loglik"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            second_moment(np.empty((0, 3)))
        with self.assertRaises(np.linalg.LinAlgError):
            second_moment(np.zeros((5, 3)))
        with self.assertRaises(ValueError):
            paired_norms(np.eye(3), 1, 0, 10, 1)

    def test_publication_scope_and_disjoint_subtree_folds(self):
        context = analysis.load_ce_protein_context(max_workers=1)
        edges = edge_manifest(context)
        folds = list(subtree_folds(context, edges))
        self.assertEqual(len(edges), 1_000)
        self.assertTrue(edges.canonical_generation_gap.eq(1).all())
        self.assertEqual(len(folds), 8)
        covered = np.zeros(len(edges), dtype=int)
        for _, train, test in folds:
            covered += test
            self.assertFalse(np.any(train & test))
            train_nodes = set(edges.loc[train, "tree_id"]) | set(edges.loc[train, "parent_tree_id"])
            test_nodes = set(edges.loc[test, "tree_id"]) | set(edges.loc[test, "parent_tree_id"])
            self.assertFalse(train_nodes & test_nodes)
        self.assertEqual(int(covered.sum()), 992)
        self.assertLessEqual(covered.max(), 1)


if __name__ == "__main__":
    unittest.main()
