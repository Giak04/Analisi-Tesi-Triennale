import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])

for folder in './Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/' :
    with fits.open(folder, f'*_0201_*.fits.gz') as hdul:
        d0114_000cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
        histra = np.append(histra, d0114_000cc00['RA'])
        histdec = np.append(histdec, d0114_000cc00['DEC'])
        gammaness = np.append(gammaness, d0114_000cc00['GAMMANESS'])
        hdul.close()

plt.hist2d(histra, histdec, bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.savefig('./Analisi-Preliminare/Imm_cc00Tot.png')
plt.close()
