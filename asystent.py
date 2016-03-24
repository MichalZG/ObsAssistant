import time
import ds9
import os
import glob
import astropy.io.fits as fits
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
import subprocess as sub
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import warnings
import threading
import numpy as np
import scipy.optimize as optimize
import math
from photutils import CircularAnnulus, CircularAperture
from photutils import aperture_photometry
import matplotlib.pyplot as plt
import matplotlib
import sys
from photutils import daofind
from photutils.extern.imageutils.stats import sigma_clipped_stats
import matplotlib.image as mpimg
warnings.filterwarnings('ignore')
watchdog_img = mpimg.imread('watchdog.dif')

matplotlib.use('QT4Agg')

# photometry, backgrund estimation, snr etc.
radius = 30  # radius from object center (fits coordinates) for backgroun est.
sigma = 3.0  # The number of standard deviations to use as the clipping limit.
fwhm_daofind = 2.0  # estimate FWHM for daofind
fwhm_multiplier = 2  # fwhm * fwhm_multiplier = aperture
annulis_in = 2  # value add to aperture
annulis_out = 4  # value add to aperture

# solve-filed param
time_limit = 15  # time limit before giving up
scale_low = 1.11  # arcsec
scale_high = 1.13  # arcsec
solve_radius = 0.8  # deg
solve_depth = '40,80,100,160,250'


# head Key
ra_key = 'RA'
dec_key = 'DEC'

# watchdog param
file_to_watch = ["*.fits"]
sleep_time = 1.0
# plot parameters
x_plot = 10
y_plot = 8
imageSize = 1000
# clear
files_to_rm = ['*.axy', '*.corr', '*.xyls', '*.match',
               '*.new', '*.rdls', '*.solved', '*.wcs']

##############
flux_tab = []
snr_tab = []
fwhm_tab = []
pol_tab = []
pol_flux_tab = []
plt.ion()
avg_pol = 0
##############


class WatchObs(PatternMatchingEventHandler):
    patterns = file_to_watch

    def on_created(self, event):
        print "Got it!", event.src_path
        self.file_to_open = event.src_path
        global pol_tab
        global avg_pol
        fits_coo = open_file(self.file_to_open)
        if solve_field(fits_coo):
            solve_coo, solve_file_hdr, solve_file_data = open_solve_file(
                str(self.file_to_open).split(".")[0]+".new")
            show_ds9(fits_coo, solve_coo)

            star = Star()
            x, y = star_pix(solve_file_hdr, fits_coo)
            if ((x < imageSize) and (y < imageSize)):
                #fileName = self.file_to_open.split("/")[-1]
                if fits_coo.separation(solve_coo) < 5 * u.arcmin:
                    x, y = star.centroid(x, y, solve_file_data, radius)
                    fwhm_x, fwhm_y = star.get_fwhm(x, y, radius,
                                               solve_file_data, medv=None)
                    flux, aperture_area = star.phot(x, y, radius,
                                            solve_file_data, fwhm_x, fwhm_y)
                    signal_to_noise = star.snr(x, y, radius,
                                       solve_file_data,
                                       fwhm_x, fwhm_y,
                                       flux, aperture_area)
                    flux_tab.append(flux)
                    snr_tab.append(signal_to_noise)
                    fwhm_tab.append((fwhm_x + fwhm_y)/2.)
                    if filter_check(solve_file_hdr):
                        pol_tab, avg_pol = pol_calc(flux, solve_file_hdr)
                    else:
                        pass
                    plot(flux_tab, snr_tab, fwhm_tab, pol_tab, avg_pol)

    def on_modified(self, event):
        print "Got it!", event.src_path
        self.file_to_open = event.src_path
        fits_coo = open_file(self.file_to_open)
        if solve_field(fits_coo):
            solve_coo, solve_file_hdr, solve_file_data = open_solve_file(
                str(self.file_to_open).split(".")[0]+".new")
            show_ds9(fits_coo, solve_coo)






class Star:

#code copy from GINGA https://ginga.readthedocs.org/en/latest/
    def __init__(self):
        self.lock = threading.RLock()
        self.skylevel_magnification = 1.05
        self.skylevel_offset = 40.0

    def gaussian(self, x, p):
        y = (1.0 / (p[1] * np.sqrt(2*np.pi)) *
             np.exp(-(x - p[0])**2 / (2*p[1]**2))) * p[2]
        return y

    def calc_fwhm(self, arr1d, medv=None, gauss_fn=None):
        if not gauss_fn:
            gauss_fn = self.gaussian
        N = len(arr1d)
        X = np.array(range(N))
        Y = arr1d
        if medv is None:
            medv = np.median(Y)
        Y = Y - medv
        maxv = Y.max()
        Y = Y.clip(0, maxv)
        p0 = [0, N-1, maxv]
        errfunc = lambda p, x, y: gauss_fn(x, p) - y
        with self.lock:
            p1, success = optimize.leastsq(errfunc, p0[:], args=(X, Y))
        if not success:
            raise IQCalcError("FWHM gaussian fitting failed")
        mu, sdev, maxv = p1
        fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sdev

        return (float(fwhm), float(mu), float(sdev), maxv)

    def get_fwhm(self, x, y, radius, data, medv=None):

        if medv is None:
            medv = np.median(data)

        x0, y0, xarr, yarr = self.cut_cross(x, y, radius, data)

        fwhm_x, cx, sdx, maxx = self.calc_fwhm(xarr, medv=medv)
        fwhm_y, cy, sdy, maxy = self.calc_fwhm(yarr, medv=medv)

        print 'fwhm:', fwhm_x, fwhm_y

        return fwhm_x, fwhm_y

    def centroid(self, x, y, data, radius):
        dist = 1024
        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)
        sources = daofind(arr - median, fwhm=fwhm_daofind, threshold=5.*std)
        cx, cy = x, y
        for i in sources:
            dist_temp = math.sqrt(abs(i[1]+x0-x)**2 + abs(i[2]+y0-y)**2)
            if dist_temp < dist:
                dist = dist_temp
                cx = i[1] + x0
                cy = i[2] + y0

        print 'center: ', cx + 1, cy + 1
        return (cx, cy)

    def cut_region(self, x, y, radius, data):
        n = radius
        ht, wd = data.shape
        x0, x1 = max(0, x-n), min(wd-1, x+n)
        y0, y1 = max(0, y-n), min(ht-1, y+n)
        arr = data[y0:y1+1, x0:x1+1]

        return (x0, y0, arr)

    def cut_cross(self, x, y, radius, data):
        n = radius
        ht, wd = data.shape
        x0, x1 = max(0, x-n), min(wd-1, x+n)
        y0, y1 = max(0, y-n), min(ht-1, y+n)
        xarr = data[y, x0:x1+1]
        yarr = data[y0:y1+1, x]

        return (x0, y0, xarr, yarr)

    def snr(self, x, y, radius, data, fwhm_x, fwhm_y, flux, aperture_area):
        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)
        try:
            signal_to_noise = flux / math.sqrt(flux +
                                               aperture_area * math.pow(std, 2))
        except ValueError:
            signal_to_noise = 0
        print 'snr:', signal_to_noise

        return signal_to_noise

    def phot(self, x, y, radius, data, fwhm_x, fwhm_y):

        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)

        aperture_r = ((fwhm_x + fwhm_y)/2.0) * fwhm_multiplier
        print 'aper:', aperture_r
        # r_in = aperture_r + annulis_in
        # r_out = aperture_r + annulis_out
        apertures = CircularAperture((x, y), aperture_r)
        # annulus_apertures = CircularAnnulus((x, y), r_in, r_out)
        rawflux_table = aperture_photometry(data, apertures)
        # bkgflux_table = aperture_photometry(data, annulus_apertures)
        # phot_table = hstack([rawflux_table, bkgflux_table],
        #                     table_names=['raw', 'bkg'])
        aperture_area = np.pi * aperture_r ** 2
        # annulus_area = np.pi * (r_out ** 2 - r_in ** 2)
        # bkg_sum = phot_table['aperture_sum_bkg'] *\
        #                      aperture_area / annulus_area
        # print 'bkg', phot_table['aperture_sum_bkg'] / annulus_area
        # final_sum = phot_table['aperture_sum_raw'] - bkg_sum
        final_sum = rawflux_table['aperture_sum'] - aperture_area * median
        print 'flux:', final_sum[0]
        return final_sum[0], aperture_area


def show_ds9(fits_coo, solve_coo):
    d = ds9.ds9()
    d.set("file "+str(event_handler.file_to_open).split(".")[0]+".new")
    d.set('scale zscale')
    d.set('zoom to fit')
    d.set('match frame wcs')
    d.set('regions', 'fk5; line(' + str(solve_coo.ra.deg) + ',' +
          str(solve_coo.dec.deg) + ',' + str(fits_coo.ra.deg) +
          ',' + str(fits_coo.dec.deg) + ')')
    d.set('regions', 'fk5; circle(' + str(fits_coo.ra.deg) +
          ',' + str(fits_coo.dec.deg) + ',7")')


def solve_field(fits_coo):
    solve_field_command = ['solve-field',
                           '--ra', '%s' % (fits_coo.ra.deg),
			               '--dec', '%s' % (fits_coo.dec.deg),
			               '--radius', '%1.1f' % solve_radius,
			               '--depth', solve_depth,
			               '--cpulimit', '%f' % time_limit,
			               '--scale-units', 'arcsecperpix',
			               '--scale-low', '%.5f' % scale_low,
			               '--scale-high', '%.5f' % scale_high,
	                       '--overwrite',
			               '--no-verify',
			               '--no-plots',
			               '--no-fits2fits',
                           str(event_handler.file_to_open)]
    sub.Popen(solve_field_command, stdout=sub.PIPE,
              stderr=sub.PIPE).communicate()
    if os.path.exists(str(event_handler.file_to_open).split(".")[0]+'.new'):
        return True
    else:
        print 'solve error'
        return False


def open_file(file_to_open):
    hdr = fits.getheader(file_to_open)
    print 'fits coo:', str(hdr[ra_key])+" "+str(hdr[dec_key])
    fits_coo = SkyCoord(hdr[ra_key]+" "+hdr[dec_key],
                        'icrs', unit=(u.hour, u.deg))
    return fits_coo


def open_solve_file(file_to_open):
    f = fits.open(file_to_open)
    data = f[0].data
    hdr = f[0].header
    w = WCS(hdr)
    wx, wy = w.wcs_pix2world(hdr['NAXIS1']/2, hdr['NAXIS2']/2, 1)
    solve_coo = SkyCoord(wx, wy, unit='deg')
    return solve_coo, hdr, data.astype(int)


def star_pix(solve_file_hdr, fits_coo):
    w = WCS(solve_file_hdr)
    px, py = w.wcs_world2pix(fits_coo.ra.deg, fits_coo.dec.deg, 1)

    return int(px), int(py)


def pol_calc(flux, solve_file_hdr):
    global avg_pol

    pol_nr = [int(s) for s in solve_file_hdr['FILTER'] if s.isdigit()][0]
    if pol_nr == int(len(pol_flux_tab)+1):
        pol_flux_tab.append(flux)
    else:
        del pol_flux_tab[:]

    if len(pol_flux_tab) == 4:
        u_st = (pol_flux_tab[1] - pol_flux_tab[0]) /\
               (pol_flux_tab[1] + pol_flux_tab[0])
        q_st = (pol_flux_tab[3] - pol_flux_tab[2]) /\
               (pol_flux_tab[3] + pol_flux_tab[2])
        pd = np.sqrt(np.power(u_st, 2) + np.power(q_st, 2))
        pol_tab.append(pd * 100)
        avg_pol = reduce(lambda a, b: a + b, pol_tab) / len(pol_tab)
        del pol_flux_tab[:]

    return pol_tab, avg_pol


def filter_check(solve_file_hdr):
    filter_name = solve_file_hdr['FILTER']
    if 'P' in filter_name:
        return True
    else:
        return False


def plot(flux_tab, snr_tab, fwhm_tab, pol_tab, avg_pol):
    plt.figure(1, figsize=(x_plot, y_plot))
    ax0 = plt.subplot(511)
    ax0 = plt.imshow(watchdog_img)
    ax0.axes.get_xaxis().set_visible(False)
    ax0.axes.get_yaxis().set_visible(False)
    ax1 = plt.subplot(512)
    ax1.set_title('FLUX [counts] - %.2f' % (flux_tab[-1]))
    ax1.set_xlim(-0.2, len(flux_tab) - 0.8)
    ax1.set_ylim(min(flux_tab) * 0.9, max(flux_tab) * 1.1)
    ax1.plot(flux_tab, 'ro')
    ax2 = plt.subplot(513)
    ax2.set_title('SNR - %.2f' % (snr_tab[-1]))
    ax2.set_xlim(-0.2, len(snr_tab) - 0.8)
    ax2.set_ylim(min(snr_tab) * 0.9, max(snr_tab) * 1.1)
    ax2.plot(snr_tab, 'ro')
    ax3 = plt.subplot(514)
    ax3.set_title('FWHM [pix] - %.2f' % (fwhm_tab[-1]))
    ax3.set_xlim(-0.2, len(fwhm_tab) - 0.8)
    ax3.set_ylim(min(fwhm_tab) * 0.9, max(fwhm_tab) * 1.1)
    ax3.plot(fwhm_tab, 'ro')
    ax4 = plt.subplot(515)
    if pol_tab:
        ax4.cla()
        ax4.set_title('PD [%]' + ' - %.2f' % avg_pol)
        ax4.set_xlim(-0.2, len(pol_tab) - 0.8)
        ax4.set_ylim(min(pol_tab) * 0.9, max(pol_tab) * 1.1)
        ax4.axhline(y=avg_pol, xmin=0, xmax=len(pol_tab) * 1.1,
                    c="blue", linewidth=0.5, zorder=0)
    ax4.plot(pol_tab, 'ro')

    plt.tight_layout()
    plt.draw()
    plt.pause(0.0001)

def clear():
    print 'cleaning.....'
    for i in files_to_rm:
        files = glob.glob(path + i)
        print path + i
        for j in files:
            os.remove(j)

if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) > 0:
        path = args[0]
    else:
        print "error - write path to files"

    event_handler = WatchObs()
    observer = Observer()
    observer.schedule(event_handler, path=path, recursive=False)
    observer.start()
    try:
        while True:
            time.sleep(sleep_time)
    except KeyboardInterrupt:
        clear()
        observer.stop()
        observer.join()