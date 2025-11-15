import numpy as np
import astropy.io.fits
import matplotlib.pyplot as plt
import scipy.ndimage
import astropy.stats
import warnings


def readfits(fitspath):
    hdu = astropy.io.fits.open(fitspath)[0]
    header, data = hdu.header, hdu.data.astype(np.float32)
    data[np.where(data == 0)] = np.nan
    return header, data


def _clippedmean(data, axis=None, sigma=3):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean, median, stdde = astropy.stats.sigma_clipped_stats(
            data, axis=axis, sigma=sigma, cenfunc=np.nanmedian
        )
    return mean


def channelcorrection(
    data, order=1, nyref=4, verbose=False, vrange=50, correctacn=True
):

    def correctone(sx):
        ref0 = _clippedmean(data[sy0, sx])
        ref1 = _clippedmean(data[sy1, sx])
        if verbose:
            print("ref0 = %+.2f ref1 = %+.2f" % (ref0, ref1))
        if order == 0:
            ref = 0.5 * (ref0 + ref1)
        else:
            fy = np.arange(0, data.shape[0]) / (data.shape[0] - 1)
            ref = ref0 * (1 - fy) + ref1 * fy
            ref = ref[:, np.newaxis]
        refdata[:, sx] -= ref

    assert order == 0 or order == 1
    sy0 = slice(0, nyref)
    sy1 = slice(-nyref - 1, -1)
    refdata = data.copy()
    for i in range(0, 32):
        if correctacn:
            correctone(slice(i * 64 + 0, (i + 1) * 64, 2))
            correctone(slice(i * 64 + 1, (i + 1) * 64, 2))
        else:
            correctone(slice(i * 64, (i + 1) * 64))
        if verbose:
            plt.title("Channel %d Lower Reference" % i)
            plt.imshow(data[sy0, sx] - ref0, vmin=-0.5 * vrange, vmax=+0.5 * vrange)
            plt.show()
            plt.title("Channel %d Upper Reference" % i)
            plt.imshow(data[sy1, sx] - ref1, vmin=-0.5 * vrange, vmax=+0.5 * vrange)
            plt.show()
    return refdata


def rowcorrection(data, verbose=False, nxref=4, nfilter=None):
    if verbose:
        iy = np.arange(0, data.shape[0])
        plt.plot(iy[600:900], np.nanmedian(data, axis=1)[600:900] - np.nanmedian(data))
        plt.title("Y Correction")
    if nxref < 0:
        refy = _clippedmean(data[:, 1024 + nxref : 1024 - nxref], axis=1)
    else:
        refy = _clippedmean(np.roll(data, nxref, axis=1)[:, 0 : 2 * nxref], axis=1)
    if nfilter is not None:
        refy = scipy.ndimage.median_filter(refy, nfilter)
    refy -= _clippedmean(refy)
    if verbose:
        plt.plot(iy[600:900], refy[600:900])
        # plt.plot(iy[600:900], np.nanmedian(data - refy[:, np.newaxis], axis=1)[600:900])
        plt.title("Y Correction")
        plt.show()
    return data - refy[:, np.newaxis]
