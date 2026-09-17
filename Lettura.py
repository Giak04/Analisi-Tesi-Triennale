import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import gammapy.utils.regions
from gammapy.data import EventList
D0114cc00_00 = '$Analisi-Gammmapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc00.lv3.fits'
events = EventList.read(D0114cc00_00)
