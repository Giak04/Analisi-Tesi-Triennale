import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
#import gammapy.utils.regions
#from gammapy.data import EventList
#D0114cc00_00 = '$Analisi-Gammmapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc00.lv3.fits'
#events = EventList.read(D0114cc00_00)
#events.peek()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0502_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    #hdul[1].data #mostra la prima tabella (per ASTRI sono gli eventi), così restituisce una lista di liste
    data = Table(hdul[1].data) #mostra la prima tabella come tabella
    hdul.close()

#guarda le funzioni disponibili per Table, in genere i dati più interessanti sono nomi delle colonne, tipi delle colonne e lunghezza
print(data.columns)
