import numpy as np
from scipy.special import wofz


def gaussian(x, A, sigma, mu):
    """
        Model a Gaussian distribution.

        Parameters:
            x (array-like): Independent variable (e.g., wavelength or frequency).
            A (float): Amplitude of the distribution
            mu (float): Central position of the distribution
            sigma (float): standard deviation
    .

        Returns:
            array-like: Gaussian distribution values at the specified x values.
    """
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def voigt(x, A, sigma, gamma, mu):
    """
    Model a Voigt profile scaled to Lorentz peak amplitude.
    (same scaling used as astropy)

    Parameters:
        x (array-like): Independent variable (e.g., wavelength or frequency).
        mu (float): Central position of the profile.
        sigma (float): Gaussian standard deviation.
        gamma (float): Lorentzian line width.
        A (float): Amplitude of the profile.

    Returns:
        array-like: Voigt profile values at the specified x values.
    """
    z = ((x - mu) + 1j * gamma) / (sigma * np.sqrt(2))
    V = A * wofz(z).real * np.sqrt(np.pi / 2) * gamma / sigma
    return V


def lorentzian(x, mu, gamma):
    """
    Computes the normalized Lorentzian distribution.

    Parameters:
    - x: NumPy array or scalar, input values.
    - mu: float, location parameter (peak position).
    - gamma: float, scale parameter.

    Returns:
    NumPy array, probability density function values corresponding to input x.
    """
    return ((gamma / 2) / ((x - mu) ** 2 + (gamma / 2) ** 2)) / np.pi


def pseudo_voigts(x, nu, A=1, FWHM=1, mu=0):
    """
    Model of two Voigt profile scaled to Lorentz peak amplitude.
    (scaled same as astropy)

    Parameters:
        x (array-like): Independent variable (e.g., wavelength or frequency).
        nu (float): Mixing coefficient (from 0-1) of the Gaussian and Lorentzian equation, 1 is full Gaussian and 0 is fully Lorentzian
        A (float): Amplitude of the distribution
        FWHM (float): Full width at half maximum of the profile
        mu (float): Central position of the distribution

    Returns:
        array-like: Voigt profile values at the specified x values.
    """
    sigma = FWHM / np.sqrt(2 * np.log(2))
    gamma = FWHM * 2
    ag = 1 / (sigma * np.sqrt(2 * np.pi))
    I = A / (nu * ag + (1 - nu) * 2 / (np.pi * gamma))

    return I * (
        (nu * gaussian(x, ag, sigma, mu)) + ((1 - nu) * lorentzian(x, mu, gamma))
    )


def two_gaussians(x, A1, A2, sigma1, sigma2, mu1, mu2):
    """
    Model of two Gaussian distributions.

    Parameters:
        x (array-like): Independent variable (e.g., wavelength or frequency).
        A1 (float): Amplitude of the  first distribution
        A2 (float): Amplitude of the  second distribution
        sigma1 (float): standard deviation of the first distribution
        sigma2 (float): standard deviation of the second distribution
        mu1 (float): Central position of the first distribution
        mu2 (float): Central position of the second distribution


    Returns:
        array-like: Gaussian distribution values at the specified x values.
    """
    y = (A1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)) + (
        A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2)
    )

    return y


def two_voigts(x, A1, A2, sigma1, sigma2, gamma1, gamma2, mu1, mu2):
    """
    Model of two Voigt profile scaled to Lorentz peak amplitude.
    (scaled same as Astropy)

    Parameters:
        x (array-like): Independent variable (e.g., wavelength or frequency).
        A1 (float): Amplitude of the  first distribution
        A2 (float): Amplitude of the  second distribution
        sigma1 (float): Gaussian standard deviation of the first distribution
        sigma2 (float): Gaussian standard deviation of the second distribution
        gamma1 (float): Lorentzian line width of the first distribution
        gamma2 (float): Lorentzian line width of the second distribution
        mu1 (float): Central position of the first distribution
        mu2 (float): Central position of the second distribution

    Returns:
        array-like: Voigt profile values at the specified x values.
    """
    z1 = ((x - mu1) + 1j * gamma1) / (sigma1 * np.sqrt(2))
    z2 = ((x - mu2) + 1j * gamma2) / (sigma2 * np.sqrt(2))
    y = (A1 * wofz(z1).real * np.sqrt(np.pi / 2) * gamma1 / sigma1) + (
        A2 * wofz(z2).real * np.sqrt(np.pi / 2) * gamma2 / sigma2
    )

    return y


def two_pseudo_voigts(x, nu1, nu2, A1, A2, FWHM1, FWHM2, mu1, mu2):
    """
    Model of two Voigt profile scaled to Lorentz peak amplitude.
    (scaled same as astropy)

    Parameters:
        x (array-like): Independent variable (e.g., wavelength or frequency).
        nu1 (float): Mixing coefficient (from 0-1) of the Gaussian and Lorentzian equation, 1 is fully Gaussian and 0 is fully Lorentzian
        nu2 (float): Mixing coefficient (from 0-1) of the Gaussian and Lorentzian equation, 1 is full Gaussian and 0 is fully Lorentzian
        A1 (float): Amplitude of the  first distribution
        A2 (float): Amplitude of the  second distribution
        FWHM1 (float): Full width at half maximum of the first distribution
        FWHM2 (float): Full width at half maximum of the second distribution
        mu1 (float): Central position of the first distribution
        mu2 (float): Central position of the second distribution

    Returns:
        array-like: Voigt profile values at the specified x values.
    """
    sigma1 = FWHM1 / np.sqrt(2 * np.log(2))
    sigma2 = FWHM2 / np.sqrt(2 * np.log(2))
    gamma1 = FWHM1 * 2
    gamma2 = FWHM2 * 2
    ag1 = 1 / (sigma1 * np.sqrt(2 * np.pi))
    ag2 = 1 / (sigma2 * np.sqrt(2 * np.pi))
    I1 = A1 / (nu1 * ag1 + (1 - nu1) * 2 / (np.pi * gamma1))
    I2 = A2 / (nu2 * ag2 + (1 - nu2) * 2 / (np.pi * gamma2))
    return I1 * (
        (nu1 * gaussian(x, ag1, sigma1, mu1)) + ((1 - nu1) * lorentzian(x, mu1, gamma1))
    ) + I2 * (
        (nu2 * gaussian(x, ag2, sigma2, mu2)) + ((1 - nu2) * lorentzian(x, mu2, gamma2))
    )
