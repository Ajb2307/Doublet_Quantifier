import numpy as np
import numpy.ma as ma


def bin_spectra(spectra, bin_size, take="max"):
    """
    Bin spectra data and compute the maximum (or median or minimum) value in each bin and return a new binned spectra.

    Parameters:
    - spectra: array,  of spectra
        spectra[0, :] (wavelengths or frequencies)
        spectra[1, :] (spectral intensities)
        spectra[2, :] (spectral intensities uncertainty)
    - y: array, y-values (spectral intensities)
    - bin_size: int, size of each bin
    - take: str ("max", "min", "med"), which value represents the bin

    Returns:
    - spectra now binned with the given value in take
    """
    # Determine the number of bins
    num_bins = len(spectra[0, :]) // bin_size

    # Reshape data to create bins
    binned_x = spectra[0, : num_bins * bin_size].reshape((num_bins, bin_size))
    binned_y = spectra[1, : num_bins * bin_size].reshape((num_bins, bin_size))
    binned_e = spectra[2, : num_bins * bin_size].reshape((num_bins, bin_size))

    if take == "max":
        # Compute the maximum value in each bin
        binned_y = np.max(binned_y, axis=1)
    elif take == "min":
        # Compute the maximum value in each bin
        binned_y = np.min(binned_y, axis=1)
    elif take == "med":
        # Compute the maximum value in each bin
        binned_y = np.median(binned_y, axis=1)

    binned_spectra = np.full((3, num_bins), np.NaN)
    binned_spectra[0, :] = np.median(binned_x, axis=1)
    binned_spectra[1, :] = binned_y
    binned_spectra[2, :] = np.median(binned_e, axis=1)

    return binned_spectra


def linear_continuum(spec, dip, binned=False, bin_take="max", bin_size=10):
    """

    fits line returns values of slope (m) and intercept (b)

    spec: spectrum where you want to fit line
    dips: list containing starting and stopping point of the absorption:
            [absorption starts, absorption stops]
        or list of 4 dip locations:
            [first absorption starts, first absorption stops, second absorption starts, second absorption stops]
            Note: number_of_dips = 2 for this to work


    Returns:
        [m, b], covariance matrix
    """

    continuum_mask = np.full(spec[0, :].shape, True)  # empty continuum mask

    # create mask for where continuum will be fit
    for i in range(len(dip) // 2):
        outside_dip = ma.masked_outside(spec[0, :], dip[i * 2], dip[i * 2 + 1]).mask
        continuum_mask = np.logical_and(continuum_mask, outside_dip)

    # creating continuum spectra array
    continuum_spectra = np.full((3, continuum_mask.sum()), np.NaN)
    continuum_spectra[0, :] = spec[0, :][continuum_mask]
    continuum_spectra[1, :] = spec[1, :][continuum_mask]
    continuum_spectra[2, :] = spec[2, :][continuum_mask]

    # binning if applicable
    if binned:
        continuum_spectra = bin_spectra(continuum_spectra, bin_size, take=bin_take)

    [m, b], cov = np.polyfit(
        continuum_spectra[0, :],
        continuum_spectra[1, :],
        1,
        w=1 / continuum_spectra[2, :],
        cov=True,
    )

    return [m, b], cov


def polynomial_continuum(spec, dip, deg=2, binned=False, bin_take="max", bin_size=10):
    """

    fits polynomial equations outside dip(s)
    returns coefficients (in descending degree order), and covariance matrix


    spec: spectrum where you want to fit line
    dips: list containing starting and stopping point of the absorption:
            [absorption starts, absorption stops]
        or list of 4 dip locations:
            [first absorption starts, first absorption stops, second absorption starts, second absorption stops]
            Note: number_of_dips = 2 for this to work
    deg: int degree of line fit


    Returns:
        ndarray, shape (deg + 1,) or (deg + 1, K) Polynomial coefficients, highest power first.
        covariance matrix
    """

    continuum_mask = np.full(spec[0, :].shape, True)  # empty continuum mask

    for i in len(dip) / 2:
        dip = ma.masked_outside(spec[0, :], dip[i * 2], dip[i * 2 + 1]).mask
        continuum_mask = np.logical_and(continuum_mask, dip.mask)

    # creating continuum spectra array
    continuum_spectra = np.full((3, continuum_mask.sum()), np.NaN)
    continuum_spectra[0, :] = spec[0, :][continuum_mask]
    continuum_spectra[1, :] = spec[1, :][continuum_mask]
    continuum_spectra[2, :] = spec[2, :][continuum_mask]

    # binning if applicable
    if binned:
        continuum_spectra = bin_spectra(continuum_spectra, bin_size, take=bin_take)

    coef, cov = np.polyfit(
        continuum_spectra[0, :],
        continuum_spectra[1, :],
        deg,
        w=1 / continuum_spectra[2, :],
        cov=True,
    )

    return coef, cov
