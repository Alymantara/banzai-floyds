

def test_gauss_hermite_analytic_derivatives_match_finite_differences():
    # the line fits of the arcs use these analytic Jacobians instead of finite differences (half the arc time);
    # the explicit Hermite polynomials must equal scipy's, and the Jacobians the numerical derivatives
    import numpy as np
    from scipy.optimize._numdiff import approx_derivative
    from scipy.special import eval_hermite, factorial
    from banzai_floyds.utils import fitting_utils
    from banzai_floyds.utils.wavelength_utils import gauss_hermite_residuals, gauss_hermite_jacobian
    from banzai_floyds.wavelengths import _blend_residuals, _blend_jacobian
    w = np.linspace(-6, 6, 301)
    for n in (3, 4):
        expected = eval_hermite(n, w) / np.sqrt(2.0 ** n * factorial(n))
        np.testing.assert_allclose(fitting_utils._normalized_hermite(n, w), expected, rtol=1e-12, atol=1e-10)
    rng = np.random.default_rng(1)
    x = np.arange(40.0)
    flux, error = rng.normal(size=40), np.full(40, 0.8)
    params = np.array([20.3, 150.0, 5.0, 2.4, 0.06, -0.04])
    numerical = approx_derivative(gauss_hermite_residuals, params, method='3-point', args=(x, flux, error))
    np.testing.assert_allclose(gauss_hermite_jacobian(params, x, flux, error), numerical, rtol=1e-6, atol=1e-6)
    args = (x, flux, error, np.array([-3.1, 0.0, 4.2]), 2.2, 0.04, -0.02)
    params = np.array([19.7, 4.0, 80.0, 200.0, 50.0])
    numerical = approx_derivative(_blend_residuals, params, method='3-point', args=args)
    np.testing.assert_allclose(_blend_jacobian(params, *args), numerical, rtol=1e-6, atol=1e-6)
