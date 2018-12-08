import time
import pyds9 as ds9
import os
import glob
from collections import deque
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
import pygame
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import matplotlib.image as mpimg
from pymodbus.client.sync import ModbusTcpClient as ModbusClient
from pymodbus.constants import Endian
from pymodbus.payload import BinaryPayloadBuilder
from pymodbus.payload import BinaryPayloadDecoder
from collections import OrderedDict
import struct


warnings.filterwarnings('ignore')
watchdog_img = mpimg.imread('/opt/ObsAssistant/watchdog.dif')

matplotlib.use('QT4Agg')

# photometry, backgrund estimation, snr etc.
radius = 30  # radius from object center (fits coordinates) for backgroun est.
sigma = 2.0  # The number of standard deviations to use as the clipping limit.
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
filter_key = 'FILTER'
object_key = 'OBJECT'
exp_key = 'EXPTIME'
time_key = 'TIME-OBS'

# watchdog param
file_to_watch = ["*.fits"]
sleep_time = 1.0
max_sleep = 300
sleep_dur = 0
# plot parameters
max_plot_len = 50
x_plot = 10
y_plot = 8
imageSize = 1000
# clear
files_to_rm = ['*.axy', '*.corr', '*.xyls', '*.match',
               '*.new', '*.rdls', '*.solved', '*.wcs']


# MODBUS
modbus_ip = '192.168.2.16'
modbus_port = 502
modbus_UNIT = 0x1
c_registers = 2


##############
flux_tab = deque(maxlen=max_plot_len)
snr_tab = deque(maxlen=max_plot_len)
fwhm_tab = deque(maxlen=max_plot_len)
pol_tab = deque(maxlen=max_plot_len)
pol_flux_tab = deque(maxlen=max_plot_len)
plt.ion()
avg_pol = 0
im_object_name = ''
im_counter = 0
plot_clear = False
##############


class WatchObs(PatternMatchingEventHandler):
    patterns = file_to_watch

    def on_created(self, event):
        print("Got it!", event.src_path)
        self.file_to_open = event.src_path
        global pol_tab
        global avg_pol
        global sleep_dur
        global im_object_name
        global im_counter
        global plot_clear
        sleep_dur = 0
        fits_coo = open_file(self.file_to_open)
        if solve_field(fits_coo):
            im_counter += 1
            solve_coo, solve_file_hdr, solve_file_data = open_solve_file(
                str(self.file_to_open).split(".")[0]+".new")
            if modbus_client:
                if write_modbus(solve_coo):
                    print('MODBUS UPDATE')
                else:
                    print('MODBUS FAIL!')
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
                    signal_to_noise, bkg_median = star.snr(x, y, radius,
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
                    solve_file_object_name = solve_file_hdr[object_key]

                    plot(flux_tab, snr_tab, fwhm_tab, pol_tab, avg_pol,
                         solve_file_hdr, plot_clear, im_counter, bkg_median)
"""
    def on_modified(self, event):
        print("Got it!", event.src_path)
        global sleep_dur
        sleep_dur = 0
        self.file_to_open = event.src_path
        fits_coo = open_file(self.file_to_open)
        if solve_field(fits_coo):
            solve_coo, solve_file_hdr, solve_file_data = open_solve_file(
                str(self.file_to_open).split(".")[0]+".new")
            show_ds9(fits_coo, solve_coo)
"""

class Star:

    # code copy from GINGA https://ginga.readthedocs.org/en/latest/
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

        def errfunc(p, x, y): return gauss_fn(x, p) - y
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

        print('fwhm:', fwhm_x, fwhm_y)

        return fwhm_x, fwhm_y

    def centroid(self, x, y, data, radius):
        dist = 1024
        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)
        daofind = DAOStarFinder(threshold=5.*std, fwhm=fwhm_daofind)
        sources = daofind.find_stars(arr - median)
        cx, cy = x, y
        for i in sources:
            dist_temp = math.sqrt(abs(i[1]+x0-x)**2 + abs(i[2]+y0-y)**2)
            if dist_temp < dist:
                dist = dist_temp
                cx = i[1] + x0
                cy = i[2] + y0

        print('center: ', cx + 1, cy + 1)
        return (cx, cy)

    def cut_region(self, x, y, radius, data):
        n = radius
        ht, wd = data.shape
        x0, x1 = max(0, x-n), min(wd-1, x+n)
        y0, y1 = max(0, y-n), min(ht-1, y+n)
        arr = data[int(y0):int(y1)+1, int(x0):int(x1)+1]

        return (x0, y0, arr)

    def cut_cross(self, x, y, radius, data):
        n = radius
        ht, wd = data.shape
        x0, x1 = max(0, x-n), min(wd-1, x+n)
        y0, y1 = max(0, y-n), min(ht-1, y+n)
        xarr = data[int(y), int(x0):int(x1)+1]
        yarr = data[int(y0):int(y1)+1, int(x)]

        return (x0, y0, xarr, yarr)

    def snr(self, x, y, radius, data, fwhm_x, fwhm_y, flux, aperture_area):
        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)
        try:
            signal_to_noise = flux / math.sqrt(flux +
                                               aperture_area * math.pow(std, 2))
        except ValueError:
            signal_to_noise = 0
        print('snr:', signal_to_noise)

        return signal_to_noise, median

    def phot(self, x, y, radius, data, fwhm_x, fwhm_y):

        x0, y0, arr = self.cut_region(x, y, radius, data)
        mean, median, std = sigma_clipped_stats(arr, sigma=sigma)

        aperture_r = ((fwhm_x + fwhm_y)/2.0) * fwhm_multiplier
        print('aper:', aperture_r)
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
        print('flux:', final_sum[0])
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
                           str(event_handler.file_to_open)]
    sub.Popen(solve_field_command, stdout=sub.PIPE,
              stderr=sub.PIPE).communicate()
    if os.path.exists(str(event_handler.file_to_open).split(".")[0]+'.new'):
        return True
    else:
        print('solve error')
        return False


def open_file(file_to_open):
    hdr = fits.getheader(file_to_open)
    print('fits coo:', str(hdr[ra_key])+" "+str(hdr[dec_key]))
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
        pol_tab.append(0)
        avg_pol = reduce(lambda a, b: a + b, pol_tab) / len(pol_tab)
        del pol_flux_tab[:]

    return pol_tab, avg_pol


def filter_check(solve_file_hdr):
    filter_name = solve_file_hdr['FILTER']
    if 'P' in filter_name:
        return True
    else:
        return False


def plot(flux_tab, snr_tab, fwhm_tab, pol_tab, avg_pol, hdr, clear, im_counter, bkg_median):
    plt.figure(1, figsize=(x_plot, y_plot))
    if clear:
        plt.clf()
    plt.clf()
    ax01 = plt.subplot(521)
    ax01.cla()
    time.sleep(0.1)
    ax01.text(0.1, 0.8, "Object: "+hdr[object_key], fontsize=14)
    ax01.text(0.1, 0.55, "Filter: "+hdr[filter_key], fontsize=14)
    ax01.text(0.1, 0.3, "Exptime: "+str(hdr[exp_key]), fontsize=14)
    ax01.text(0.1, 0.05, "BKG: "+str(bkg_median), fontsize=14)
    ax01.text(0.5, 0.8, "NUM: "+str(im_counter),
              fontsize=15, fontweight='bold')
    ax01.text(0.5, 0.4, str(hdr[time_key][:8]),
              fontsize=15, fontweight='bold')

    ax01.axes.get_xaxis().set_visible(False)
    ax01.axes.get_yaxis().set_visible(False)
    ax02 = plt.subplot(522)
    ax02.imshow(watchdog_img)
    ax02.axes.get_xaxis().set_visible(False)
    ax02.axes.get_yaxis().set_visible(False)
    ax1 = plt.subplot(512)
    ax1.set_title('FLUX [counts]: %.2f' % (flux_tab[-1]))
    ax1.set_xlim(-0.2, len(flux_tab) - 0.8)
    ax1.set_ylim(min(flux_tab) * 0.9, max(flux_tab) * 1.1)
    ax1.plot(flux_tab, 'ro')
    ax2 = plt.subplot(513)
    ax2.set_title('SNR: %.2f' % (snr_tab[-1]))
    ax2.set_xlim(-0.2, len(snr_tab) - 0.8)
    ax2.set_ylim(min(snr_tab) * 0.9, max(snr_tab) * 1.1)
    ax2.plot(snr_tab, 'ro')
    ax3 = plt.subplot(514)
    ax3.set_title('FWHM [pix]: %.2f' % (fwhm_tab[-1]))
    ax3.set_xlim(-0.2, len(fwhm_tab) - 0.8)
    ax3.set_ylim(min(fwhm_tab) * 0.9, max(fwhm_tab) * 1.1)
    ax3.plot(fwhm_tab, 'ro')

    plt.tight_layout()
    plt.draw()
    plt.pause(0.01)


def clear():
    print('cleaning.....')
    for i in files_to_rm:
        files = glob.glob(path + i)
        print(path + i)
        for j in files:
            os.remove(j)


def connect_modbus(modbus_ip, modbus_port):
    try:
        client = ModbusClient(modbus_ip, port=modbus_port)
    except:
        return None
    if client.connect():
        return client
    return None


def read_modbus():
    pass


def write_modbus(solve_coo):

    if solve_coo is None:
        return False
        
    ra = solve_coo.ra.deg
    if ra > 180:
        ra -= 360
    ra = ra * np.pi/180.
    dec = solve_coo.dec.deg * np.pi/180.

    val_dict = {
        24592: ra,
        24590: dec,
        24594: time.time()
    }

    for address, value in val_dict.items():
        builder = BinaryPayloadBuilder(byteorder=Endian.Big,
                                       wordorder=Endian.Big)
        builder.add_32bit_float(value)
        payload = builder.build()
        registers = builder.to_registers()
        rr = modbus_client.write_registers(address, registers, unit=modbus_UNIT)
        time.sleep(0.1)

        if rr.isError():
            return False

    return True


if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) > 0:
        path = args[0]
    else:
        print("error - write path to files")

    modbus_client = connect_modbus(modbus_ip, modbus_port)
    event_handler = WatchObs()
    observer = Observer()
    observer.schedule(event_handler, path=path, recursive=False)
    observer.start()
    pygame.init()
    try:
        while True:
            time.sleep(sleep_time)
            sleep_dur += 1
            if sleep_dur >= max_sleep:
                pygame.mixer.music.load('/opt/ObsAssistant/Angry-dog.mp3')
                pygame.mixer.music.play()
                time.sleep(7)
            else:
                pygame.mixer.music.stop()

    except KeyboardInterrupt:
        clear()
        observer.stop()
        observer.join()
