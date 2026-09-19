import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0502_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_00cc00 = Table(hdul[3].data) #mostra la prima tabella come tabella
    hdul.close()

print(d0114_00cc00.columns)
print(d0114_00cc00['ENERG_LO'])
print(d0114_00cc00['ENERG_HI'])
print(d0114_00cc00['THETA_LO'])
print(d0114_00cc00['THETA_HI'])
print(d0114_00cc00['MIGRA_LO'])
print(d0114_00cc00['MIGRA_HI'])
print(d0114_00cc00['MATRIX'])

plt.hist(d0114_00cc00['THETA_LO'], bins=50, edgecolor='black')
plt.savefig('ProvaIRF2')
plt.close()


