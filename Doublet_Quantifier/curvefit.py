from Doublet_Quantifier.models import *
from Doublet_Quantifier.continuum import *
from scipy.optimize import curve_fit
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt


def fit_single_curve(
    spec, cont, dip, function="gaussian", show=False, continuum="linear"
):
    """
    spec (array-like): spectrum
    function (string): takes "gaussian" or "voigt" as function to be fit to the spectrum
    cont (list):  first value is start of continuum and second is end of continuum
    dip (list): first value is start of absorption and second is end of absorption
    continuum (string): current only working option is "linear"


    Returns:
        continuum parameters (array-like): ex: line params = (m, b, σ_m, σ_b)
        curve parameters (array-like):  ex: Gaussian params = (A, a, mu)
        curve parameters standard deviations (array-like):  ex: Gaussian params = (σ_A, σ_a, σ_mu)
        if show == True then returns a graph of spectra and fit
    """
    # getting desired section of m
    continuum_mask = ma.masked_inside(spec[0, :], cont[0], cont[1])
    # making new array
    spec_section = np.full((3, ma.count_masked(continuum_mask)), np.NaN)
    spec_section[0, :] = spec[0, :][continuum_mask.mask]
    spec_section[1, :] = spec[1, :][continuum_mask.mask]
    spec_section[2, :] = spec[2, :][continuum_mask.mask]

    if continuum == "linear":
        # getting continuum best fit line
        [m, b], cov_lin = linear_continuum(spec_section, dip)
        cont_parameters = [m, b, cov_lin[0, 0], cov_lin[1, 1]]

        # normalizing the plot
        norm_spec_section = np.full(spec_section.shape, np.NaN)
        norm_spec_section[0, :] = spec_section[0, :]
        norm_spec_section[1, :] = spec_section[1, :] - (m * spec_section[0, :] + b)
        norm_spec_section[2, :] = spec_section[2, :]

    # getting each dip region expanded to include extra region on each side of dip
    center = (dip[0] + dip[1]) / 2
    dip_width = abs(dip[0] - dip[1])
    fitting_region = [center - 1.5 * dip_width / 2, center + 1.5 * dip_width / 2]
    fit_mask = ma.masked_inside(
        norm_spec_section[0, :], fitting_region[0], fitting_region[1]
    ).mask

    # fitting Gaussians
    if function == "gaussian":
        # initial guesses in process
        sigma0 = dip_width / 10
        A0 = 0.75 * np.min(norm_spec_section[1, :])
        mu0 = center

        # regions for each variable
        sigma_bounds = [dip_width / 20, dip_width / 2]
        A_bounds = [
            np.min(norm_spec_section[1, :]) * 2.5,
            np.min(norm_spec_section[1, :]) / 5,
        ]
        mu_bounds = [dip[0], dip[1]]

        # fitting using scipy
        params, covariance = curve_fit(
            gaussian,
            norm_spec_section[0, :][fit_mask],
            norm_spec_section[1, :][fit_mask],
            p0=[A0, sigma0, mu0],
            sigma=norm_spec_section[2, :][fit_mask],
            bounds=(
                (A_bounds[0], sigma_bounds[0], mu_bounds[0]),
                (A_bounds[1], sigma_bounds[1], mu_bounds[1]),
            ),
        )
        A_fit, sigma_fit, mu_fit = params
        params_sd = np.diagonal(covariance)

    # fitting Voigt
    if function == "voigt":
        # initial guesses in process
        sigma0 = 0.2 * dip_width / 2.35
        A0 = np.min(norm_spec_section[1, :])
        gamma0 = 0.2 * dip_width / 2
        mu0 = center

        # regions for each variable
        sigma_bounds = [0, dip_width / 2.35]
        A_bounds = [A0 * 2.5, A0 / 5]
        gamma_bounds = [0, dip_width / 2]
        mu_bounds = [dip[0], dip[1]]

        # fitting using scipy
        params, covariance = curve_fit(
            voigt,
            norm_spec_section[0, :][fit_mask],
            norm_spec_section[1, :][fit_mask],
            p0=[A0, sigma0, gamma0, mu0],
            sigma=norm_spec_section[2, :][fit_mask],
            bounds=(
                (A_bounds[0], sigma_bounds[0], gamma_bounds[0], mu_bounds[0]),
                (A_bounds[1], sigma_bounds[1], gamma_bounds[1], mu_bounds[1]),
            ),
        )
        A_fit, sigma_fit, gamma_fit, mu_fit = params
        params_sd = np.diagonal(covariance)

    if show == True:
        plt.plot(
            norm_spec_section[0, :],
            norm_spec_section[1, :],
            label="Normalized Data",
            color="k",
            alpha=0.8,
        )
        # plotting if show == True and gaussian
        if show and function == "gaussian":
            print("curve fit: A, σ, μ: ", A_fit, sigma_fit, mu_fit)
            #  print("my guess: A, σ, μ: ", A0, sigma0, mu0)
            plt.plot(
                spec_section[0],
                gaussian(spec_section[0], A0, sigma0, mu0),
                label="Gaussian start parameters",
                color="grey",
                alpha=0.5,
                linestyle="--",
            )
            plt.plot(
                spec_section[0],
                gaussian(spec_section[0], A_fit, sigma_fit, mu_fit),
                label="Fitted Gaussian Curve",
                color="red",
                alpha=0.9,
                linestyle="--",
            )

        # plotting if show == True and voigt
        if function == "voigt":
            print("curve fit: A, σ, γ, μ: ", A_fit, sigma_fit, gamma_fit, mu_fit)
            #print("my guess: A, σ,γ, μ: ", A0, sigma0, gamma0, mu0)
            plt.plot(
                spec_section[0],
                voigt(spec_section[0], A0, sigma0, gamma0, mu0),
                label="Voigt start parameters",
                color="grey",
                alpha=0.5,
                linestyle="--",
            )
            plt.plot(
                spec_section[0],
                voigt(spec_section[0], A_fit, sigma_fit, gamma_fit, mu_fit),
                label="Fitted Voigt Curve",
                color="blue",
                alpha=0.9,
                linestyle="--",
            )
        plt.legend()
        plt.show()
    return cont_parameters, params, params_sd


def two_curve_bounds(spec, cont, dip, continuum="linear", mus=None):
    """
    this function returns better constrains on where the absorption exists

    spec: spectrum
    cont: list where first value is start of continuum and second is end of continuum
    dip: list where first value is start of absorption and second is end of the second absorption
    continuum: default is "linear" fitting a linear function for the continuum, no other options as of now


    Returns:
        list of 4 continuum locations:
            [first absorption starts, first absorption stops, second absorption starts, second absorption stops]
    """

    # getting desired section of m
    continuum_mask = ma.masked_inside(spec[0, :], cont[0], cont[1])
    # making new array
    spec_section = np.full((3, ma.count_masked(continuum_mask)), np.NaN)
    spec_section[0, :] = spec[0, :][continuum_mask.mask]
    spec_section[1, :] = spec[1, :][continuum_mask.mask]
    spec_section[2, :] = spec[2, :][continuum_mask.mask]

    if continuum == "linear":
        # getting continuum best fit line
        [m, b], cov = linear_continuum(spec_section, dip)

        # normalizing the plot
        norm_spec_section = np.full(spec_section.shape, np.NaN)
        norm_spec_section[0, :] = spec_section[0, :]
        norm_spec_section[1, :] = spec_section[1, :] - (m * spec_section[0, :] + b)
        norm_spec_section[2, :] = spec_section[2, :]

    # getting each dip region expanded to include extra region on each side of dip
    center = (dip[0] + dip[1]) / 2
    dip_width = abs(dip[0] - dip[1])
    fitting_region = [center - 1.5 * dip_width / 2, center + 1.5 * dip_width / 2]
    fit_mask = ma.masked_inside(
        norm_spec_section[0, :], fitting_region[0], fitting_region[1]
    ).mask

    # fitting Gaussians
    # initial guesses in process
    sigma0_g = dip_width / 20
    A0_g = 0.75 * np.min(norm_spec_section[1, :])
    if mus == None:
        mu10_g = center - dip_width * 0.25
        mu20_g = center + dip_width * 0.25
    else:
        [mu10_g, mu20_g] = mus

    # regions for each variable
    sigma_bounds_g = [dip_width / 40, dip_width / 2]
    A_bounds_g = [
        np.min(norm_spec_section[1, :]) * 2.5,
        np.min(norm_spec_section[1, :]) / 10,
    ]
    A_bounds_g = [
        np.min(norm_spec_section[1, :]) * 100,
        np.min(norm_spec_section[1, :]) / 10,
    ]
    mu1_bounds_g = [dip[0], center]
    mu2_bounds_g = [center, dip[1]]

    # fitting using scipy
    params_g, covariance_g = curve_fit(
        two_gaussians,
        norm_spec_section[0, :][fit_mask],
        norm_spec_section[1, :][fit_mask],
        p0=[A0_g, A0_g, sigma0_g, sigma0_g, mu10_g, mu20_g],
        sigma=norm_spec_section[2, :][fit_mask],
        bounds=(
            (
                A_bounds_g[0],
                A_bounds_g[0],
                sigma_bounds_g[0],
                sigma_bounds_g[0],
                mu1_bounds_g[0],
                mu2_bounds_g[0],
            ),
            (
                A_bounds_g[1],
                A_bounds_g[1],
                sigma_bounds_g[1],
                sigma_bounds_g[1],
                mu1_bounds_g[1],
                mu2_bounds_g[1],
            ),
        ),
    )
    A1_fit_g, A2_fit_g, sigma1_fit_g, sigma2_fit_g, mu1_fit_g, mu2_fit_g = params_g

    return [
        mu1_fit_g - 4 * sigma1_fit_g,
        mu1_fit_g + 4 * sigma1_fit_g,
        mu2_fit_g - 4 * sigma2_fit_g,
        mu2_fit_g + 4 * sigma2_fit_g,
    ]


def fit_two_curves(
    spec,
    cont,
    dip,
    function="all",
    show=False,
    mus=None,
    bin_continuum=True,
    bin_size=10,
    bin_take="max",
):
    # needs editing copied from single function prior to further editing
    """
    spec (array-like): spectrum
    function (string): takes "gaussian", "voigt", or "pseudo-Voigt" as function to be fit to the spectrum
    cont (list):  first value is the starting point and second is the ending point of continuum region
    dip (list): first value is the starting point and second is the ending point of absorption region
    continuum (string): current only working option is "linear"
    mus: list of mu1, mu2 start parameters otherwise its at 25% of "dip" and 75% respectively
        Note: mu1 and m2 MUST be in numerical order
    bin_continuum (boolean value): True the continuum will be binned, False it will not be binned

    Returns:
        continuum parameters (array-like): ex: line params = (m, b, σ_m, σ_b)
        curve parameters (array-like):  ex: Gaussian params = (A, a, mu)
        curve parameters standard deviations (array-like):  ex: Gaussian params = (σ_A, σ_a, σ_mu)
        if show == True then returns a graph of spectra and fit
    """
    # getting 8 sigma region for where the absorption is
    fitting_regions = two_curve_bounds(spec, cont, dip, continuum="linear", mus=mus)

    # using to get detailed mu start parameter and bounds
    mu1_bounds = fitting_regions[0:2]
    mu2_bounds = fitting_regions[2:4]

    if mus == None:
        mu1_0 = np.average(fitting_regions[0:2])
        mu2_0 = np.average(fitting_regions[2:4])
    else:
        mu1_0, mu2_0 = mus
        if (
            mu1_bounds[0] >= mu1_0 >= mu1_bounds[1]
            or mu2_bounds[0] >= mu1_0 >= mu2_bounds[1]
        ):
            raise ValueError(
                f"The either: {mu1_0} is not within the detected range [{mu1_bounds[0]}, {mu1_bounds[1]}] \
                             \n OR \n \
                            {mu2_0} is not within the detected mu range [{mu2_bounds[0]}, {mu2_bounds[1]}]"
            )

    # getting desired section of m
    continuum_mask = ma.masked_inside(spec[0, :], cont[0], cont[1])
    # making new array
    spec_section = np.full((3, ma.count_masked(continuum_mask)), np.NaN)
    spec_section[0, :] = spec[0, :][continuum_mask.mask]
    spec_section[1, :] = spec[1, :][continuum_mask.mask]
    spec_section[2, :] = spec[2, :][continuum_mask.mask]

    # fitting linear continuum
    [m, b], cov_lin = linear_continuum(
        spec_section, dip, binned=bin_continuum, bin_size=bin_size, bin_take=bin_take
    )
    cont_parameters = [m, b, cov_lin[0, 0], cov_lin[1, 1]]

    # normalizing spectra
    norm_spec_section = np.full(spec_section.shape, np.NaN)
    norm_spec_section[0, :] = spec_section[0, :]
    norm_spec_section[1, :] = spec_section[1, :] - (m * spec_section[0, :] + b)
    norm_spec_section[2, :] = spec_section[2, :]

    # getting each dip region expanded to include extra region on each side of dip
    center = (dip[0] + dip[1]) / 2
    dip_width = abs(dip[0] - dip[1])
    fitting_region = [center - dip_width / 2, center + dip_width / 2]
    fit_mask = ma.masked_inside(
        norm_spec_section[0, :], fitting_region[0], fitting_region[1]
    ).mask

    # fitting Gaussians
    if function == "gaussian" or function == "all":
        # initial guesses in process
        sigma01_g = (fitting_regions[1] - fitting_regions[0]) / 8
        sigma02_g = (fitting_regions[3] - fitting_regions[2]) / 8
        A0_g = np.min(norm_spec_section[1, :])

        # regions for each variable
        sigma_bounds1_g = [sigma01_g / 3, sigma01_g * 3]
        sigma_bounds2_g = [sigma02_g / 3, sigma02_g * 3]

        # fix this
        if np.min(norm_spec_section[1, :]) < 0:
            A_bounds_g = [A0_g * 2, A0_g / 10]
        else:  # TODO: in functions throw an error
            print("minimum spectra value is positive")

        # fitting using scipy
        try:
            params_g, covariance_g = curve_fit(
                two_gaussians,
                norm_spec_section[0, :][fit_mask],
                norm_spec_section[1, :][fit_mask],
                p0=[A0_g, A0_g, sigma01_g, sigma02_g, mu1_0, mu2_0],
                sigma=norm_spec_section[2, :][fit_mask],
                bounds=(
                    (
                        A_bounds_g[0],
                        A_bounds_g[0],
                        sigma_bounds1_g[0],
                        sigma_bounds2_g[0],
                        mu1_bounds[0],
                        mu2_bounds[0],
                    ),
                    (
                        A_bounds_g[1],
                        A_bounds_g[1],
                        sigma_bounds1_g[1],
                        sigma_bounds2_g[1],
                        mu1_bounds[1],
                        mu2_bounds[1],
                    ),
                ),
            )
            (
                A1_fit_g,
                A2_fit_g,
                sigma1_fit_g,
                sigma2_fit_g,
                mu1_fit_g,
                mu2_fit_g,
            ) = params_g
            params_sd_g = np.diagonal(covariance_g)
        except:  # TODO: make into an error
            print(
                "error occurred here were the input bound ons A, sigma1, sigma2, mu1, and mu2 respectively",
                "\n",
                A_bounds_g,
                "\n",
                sigma_bounds1_g,
                "\n",
                sigma_bounds2_g,
                "\n",
                mu1_bounds,
                "\n",
                mu2_bounds,
            )
            return

    if function == "voigt" or function == "all":
        # testing new start parameters
        sigma0_v = dip_width / 20
        gamma0_v = dip_width / 1000
        A0_v = np.min(norm_spec_section[1, :]) * 20

        # regions for each variable
        sigma_bounds_v = [dip_width / 40, dip_width / 2]
        A_bounds_v = [A0_v * 500, A0_v / 2]
        gamma_bounds_v = [0, dip_width / 2]

        # fitting using scipy
        try:
            params_v, covariance_v = curve_fit(
                two_voigts,
                norm_spec_section[0, :][fit_mask],
                norm_spec_section[1, :][fit_mask],
                p0=[A0_v, A0_v, sigma0_v, sigma0_v, gamma0_v, gamma0_v, mu1_0, mu2_0],
                sigma=norm_spec_section[2, :][fit_mask],
                bounds=(
                    (
                        A_bounds_v[0],
                        A_bounds_v[0],
                        sigma_bounds_v[0],
                        sigma_bounds_v[0],
                        gamma_bounds_v[0],
                        gamma_bounds_v[0],
                        mu1_bounds[0],
                        mu2_bounds[0],
                    ),
                    (
                        A_bounds_v[1],
                        A_bounds_v[1],
                        sigma_bounds_v[1],
                        sigma_bounds_v[1],
                        gamma_bounds_v[1],
                        gamma_bounds_v[1],
                        mu1_bounds[1],
                        mu2_bounds[1],
                    ),
                ),
            )

            (
                A1_fit_v,
                A2_fit_v,
                sigma1_fit_v,
                sigma2_fit_v,
                gamma1_fit_v,
                gamma2_fit_v,
                mu1_fit_v,
                mu2_fit_v,
            ) = params_v
            params_sd_v = np.diagonal(covariance_v)
        except:  # TODO: make into an error
            print(
                "error occurred here were the input bound ons A, sigma, gamma, mu1, and mu2 respectively",
                "\n",
                A_bounds_v,
                "\n",
                sigma_bounds_v,
                "\n",
                gamma_bounds_v,
                "\n",
                mu1_bounds,
                "\n",
                mu2_bounds,
            )
            return

    # fitting pseudo-voigt
    if function == "pseudo-voigt" or function == "all":
        # initial guesses in process
        sigma01_p = (fitting_regions[1] - fitting_regions[0]) / 8
        sigma02_p = (fitting_regions[3] - fitting_regions[2]) / 8
        A0_p = np.min(norm_spec_section[1, :])
        nu0 = 0.5

        # regions for each variable
        sigma_bounds1_p = [sigma01_p / 3, sigma01_p * 3]
        sigma_bounds2_p = [sigma02_p / 3, sigma02_p * 3]
        nu_bounds = [0, 1]

        if np.min(norm_spec_section[1, :]) < 0:
            A_bounds_p = [A0_p * 2, A0_p / 10]
        else:  # TODO: make an error
            print("minimum spectra value is positive")
        # fitting using scipy

        try:
            params_p, covariance_p = curve_fit(
                two_pseudo_voigts,
                norm_spec_section[0, :][fit_mask],
                norm_spec_section[1, :][fit_mask],
                p0=[nu0, nu0, A0_p, A0_p, sigma01_p, sigma02_p, mu1_0, mu2_0],
                sigma=norm_spec_section[2, :][fit_mask],
                bounds=(
                    (
                        nu_bounds[0],
                        nu_bounds[0],
                        A_bounds_p[0],
                        A_bounds_p[0],
                        sigma_bounds1_p[0],
                        sigma_bounds2_p[0],
                        mu1_bounds[0],
                        mu2_bounds[0],
                    ),
                    (
                        nu_bounds[1],
                        nu_bounds[1],
                        A_bounds_p[1],
                        A_bounds_p[1],
                        sigma_bounds1_p[1],
                        sigma_bounds2_p[1],
                        mu1_bounds[1],
                        mu2_bounds[1],
                    ),
                ),
            )
            (
                nu1_fit_p,
                nu2_fit_p,
                A1_fit_p,
                A2_fit_p,
                sigma1_fit_p,
                sigma2_fit_p,
                mu1_fit_p,
                mu2_fit_p,
            ) = params_p
            params_sd_p = np.diagonal(covariance_p)

            # doing a secondary fit only adjusting mu
            eps = np.finfo(float).eps
            params_p, covariance_p = curve_fit(
                two_pseudo_voigts,
                norm_spec_section[0, :][fit_mask],
                norm_spec_section[1, :][fit_mask],
                p0=[
                    nu1_fit_p,
                    nu2_fit_p,
                    A1_fit_p,
                    A2_fit_p,
                    sigma1_fit_p,
                    sigma2_fit_p,
                    mu1_fit_p,
                    mu2_fit_p,
                ],
                sigma=norm_spec_section[2, :][fit_mask],
                bounds=(
                    (
                        nu_bounds[0],
                        nu_bounds[0],
                        A1_fit_p - abs(A1_fit_p * eps),
                        A2_fit_p - abs(A2_fit_p * eps),
                        sigma1_fit_p - eps,
                        sigma2_fit_p - eps,
                        mu1_fit_p - eps,
                        mu2_fit_p - eps,
                    ),
                    (
                        nu_bounds[1],
                        nu_bounds[1],
                        A1_fit_p + abs(A1_fit_p * eps),
                        A2_fit_p + abs(A2_fit_p * eps),
                        sigma1_fit_p + eps,
                        sigma2_fit_p + eps,
                        mu1_fit_p + eps,
                        mu2_fit_p + eps,
                    ),
                ),
            )
            (
                nu1_fit_p,
                nu2_fit_p,
                A1_fit_p,
                A2_fit_p,
                sigma1_fit_p,
                sigma2_fit_p,
                mu1_fit_p,
                mu2_fit_p,
            ) = params_p
        except:  # TODO: make an error
            print(
                "error occurred here were the input bound ons A, sigma1, sigma2, mu1, and mu2 respectively",
                "\n",
                A_bounds_p,
                "\n",
                sigma_bounds1_p,
                "\n",
                sigma_bounds2_p,
                "\n",
                mu1_bounds,
                "\n",
                mu2_bounds,
            )
            return

    # plotting if show == True
    if show == True:
        plt.plot(
            norm_spec_section[0, :],
            norm_spec_section[1, :],
            label="Normalized Data",
            color="k",
            alpha=0.8,
        )

        if function == "gaussian" or function == "all":
            print("curve fit: A1, A2, σ1, σ2, μ1, μ2: ", params_g)
            # Plotting Start Parameters
            # plt.plot(spec_section[0], two_gaussians(spec_section[0], A0_g, A0_g, sigma01_g, sigma02_g, mu1_0, mu2_0),
            #  label='Gaussian start parameters', color='grey', alpha=0.5, linestyle='--')
            plt.plot(
                spec_section[0],
                two_gaussians(
                    spec_section[0],
                    A1_fit_g,
                    A2_fit_g,
                    sigma1_fit_g,
                    sigma2_fit_g,
                    mu1_fit_g,
                    mu2_fit_g,
                ),
                label="Gaussian Curve fit",
                color="red",
                alpha=0.7,
                linestyle="--",
            )

        if function == "voigt" or function == "all":
            print("curve fit: A1, A2, σ1, σ2, γ1, γ2, μ1, μ2: ", params_v)
            # Plotting Start Parameters
            # plt.plot(spec_section[0],
            #          two_voigts(spec_section[0], A0_v, A0_v, sigma0_v,
            #                     sigma0_v, gamma0_v, gamma0_v, mu1_0, mu2_0),
            #          label='Voigt start parameters', color='grey', alpha=0.5, linestyle='--')
            plt.plot(
                spec_section[0],
                two_voigts(
                    spec_section[0],
                    A1_fit_v,
                    A2_fit_v,
                    sigma1_fit_v,
                    sigma2_fit_v,
                    gamma1_fit_v,
                    gamma2_fit_v,
                    mu1_fit_v,
                    mu2_fit_v,
                ),
                label="Fitted Voigt Curve",
                color="blue",
                alpha=0.7,
                linestyle="--",
            )

        if function == "pseudo-voigt" or function == "all":
            print("curve fit: η1, η2 A1, A2, Γ1, Γ2, μ1, μ2: ", params_p)
            # Plotting Start Parameters
            # plt.plot(spec_section[0], two_pseudo_voigts(spec_section[0], nu0, nu0, A0_p, A0_p, sigma01_p, sigma02_p, mu1_0, mu2_0),
            #  label='pseudo-voigt start parameters', color='grey', alpha=0.5, linestyle='--')
            plt.plot(
                spec_section[0],
                two_pseudo_voigts(
                    spec_section[0],
                    nu1_fit_p,
                    nu2_fit_p,
                    A1_fit_p,
                    A2_fit_p,
                    sigma1_fit_p,
                    sigma2_fit_p,
                    mu1_fit_p,
                    mu2_fit_p,
                ),
                label="pseudo-voigt fit",
                color="green",
                alpha=0.7,
                linestyle="--",
            )

        plt.legend()
    if show == False:
        plt.clf()  # clear all plots since not wanted

    if function == "all":
        return (
            cont_parameters,
            params_g,
            params_sd_g,
            params_v,
            params_sd_v,
            params_p,
            params_sd_p,
        )
    elif function == "voigt":
        return cont_parameters, params_v, params_sd_v
    elif function == "pseudo-voigt":
        return cont_parameters, params_p, params_sd_p
    elif function == "gaussian":
        return cont_parameters, params_g, params_sd_g
