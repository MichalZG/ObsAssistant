import pylab as plt
import numpy as np
import json
import base64
import time

x1 = y1 = 2048
x2 = y2 = 256
start_time = 0

t_sleep = 2  # sleep time


while True:
    t = start_time + 0.1 
    d1 = np.random.random((x1, y1))
    d2 = np.random.random((x2, y2))

    flux = np.mean(d1)
    err = flux * 0.1
    fwhm = np.std(d1)
    bkg = flux * 0.2

    d1_b64 = base64.b64encode(d1)
    d2_b64 = base64.b64encode(d2)
    
    data_pack = json.dumps({'flux': flux,
		            'err': err,
		            'fwhm': fwhm,
		            'bkg': bkg,
		            'im1': d1_b64,
			    'im2': d2_b64,
			    'time': t})

    time.sleep(t_sleep) 
