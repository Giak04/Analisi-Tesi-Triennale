import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
from decimal import Decimal

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy.table import Table
from gammapy.data import Observation, Observations, observatory_locations, EventList
from gammapy.datasets import MapDataset, MapDatasetEventSampler, FluxPointsDataset
from gammapy.estimators import ASmoothMapEstimator, TSMapEstimator, ExcessMapEstimator
from gammapy.estimators import FluxPointsEstimator
from gammapy.estimators.utils import resample_energy_edges
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker, RingBackgroundMaker, SafeMaskMaker
from gammapy.datasets import MapDataset, MapDatasetOnOff
from gammapy.makers import MapDatasetMaker, WobbleRegionsFinder
from gammapy.maps import MapAxis, RegionNDMap, WcsGeom, Map, TimeMapAxis
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ConstantSpectralModel,
    FoVBackgroundModel,
    LightCurveTemplateTemporalModel,
    Models,
    PointSpatialModel,
    PowerLawSpectralModel,
    LogParabolaSpectralModel,
    SkyModel,
)
from regions import PointSkyRegion, CircleSkyRegion

from IPython.display import display
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
    SpectrumDatasetOnOff,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
    FoVBackgroundMaker
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    SkyModel,
    create_crab_spectral_model,
    PiecewiseNormSpectralModel,
    PowerLawNormSpectralModel,
    LogParabolaSpectralModel,
)

from gammapy.visualization import plot_spectrum_datasets_off_regions

import glob, os, pathlib, itertools, warnings

#Caricare i file di dati e irf

dir0114_cc00 = './Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/'

d0114_cc00 = []
for folder in dir0114_cc00 :
    d0114_cc00.append(glob.glob(os.path.join(dir0114_cc00, '*Mrk421*_0201_*.fits.gz')))

d0114_cc00 = tuple(itertools.chain.from_iterable(d0114_cc00))


i0114_cc00 = []
for folder in dir0114_cc00 :
    i0114_cc00.append(glob.glob(os.path.join(dir0114_cc00, f'*Mrk421*_0502_*.fits.gz')))

i0114_cc00 = tuple(itertools.chain.from_iterable(i0114_cc00))

#Controllo corretta lettura

print(f"Il numero di run effettuate è {len(d0114_cc00)}; il numero di IRF è {len(i0114_cc00)}.")
if len(i0114_cc00) != len(d0114_cc00):
    warnings.warn("Il numero di run effettuate non coincide con quello delle IRF, controllare l'input e riprovare.")
else:
    print("Le run presenti corrispondono alla IRF.")

dst0114_cc00 = DataStore.from_events_files(np.sort(d0114_cc00), irfs_paths=np.sort(i0114_cc00))

dst0114_cc00.hdu_table.write("hdu-index.fits.gz", overwrite=True)
dst0114_cc00.obs_table.write("obs-index.fits.gz", overwrite=True)

obs_date = dst0114_cc00.obs_table[0]["EVENTS_FILENAME"].split("/")[-1][0:8]

#Creazione colonne AZ/ZD

ev_header = Table.read(dst0114_cc00.obs_table["EVENTS_FILENAME"][0], hdu=1)
meta = ev_header.meta
loc = EarthLocation(
            lat=28.2973583333,
            lon=-16.5106222222,
        )
#Latitudine e Longitudine prese dall'header dei file, lette direttamente attraverso fv

mjdtime = ( ((dst0114_cc00.obs_table["TSTART"] + dst0114_cc00.obs_table["TSTOP"])/2).to("d") 
           + (dst0114_cc00.obs_table["MJDREFI"] + dst0114_cc00.obs_table["MJDREFF"]) * u.d )
observing_time = Time(mjdtime, format="mjd", scale="tt")  

altazframe = AltAz(location=loc, obstime=observing_time)
dst0114_cc00.obs_table["DATE-OBS"] = observing_time.isot 
coord = SkyCoord(dst0114_cc00.obs_table["RA_PNT"], dst0114_cc00.obs_table["DEC_PNT"], unit="deg", frame="icrs")
coord_altaz = coord.transform_to(altazframe)

dst0114_cc00.obs_table["AZ_PNT"] = coord_altaz.az.deg * u.deg
dst0114_cc00.obs_table["ZD_PNT"] = coord_altaz.zen.deg * u.deg

target_position = SkyCoord(ra=166.1138, 
                           dec=38.2088, 
                           unit="deg", 
                           frame='icrs')
#Coordinate e frame presi dai file aperti con fv


selection_offset = dict(type='sky_circle',
                lon=target_position.ra.deg,
                lat=target_position.dec.deg,
                frame="icrs",
                radius=Angle(3, 'deg')
                )
#L'angolo '3' è stato preso in quanto presente nell'esempio fornito, non ho trovato questo valore negli header dei file aperti in fv

#Ho scelto di non fare selezioni particolari, non avendo ancora mai aperto questi file

print(f"The mean exposure of the runs is {np.mean(dst0114_cc00.obs_table['LIVETIME']) * u.s} \n")
print(f"The total exposure time is {(np.sum(dst0114_cc00.obs_table['LIVETIME']) * u.s).to('hr')} \n")
print(f"The total runs are {len(dst0114_cc00.obs_table['OBS_ID'])}. These are the selected ObsID: \n {dst0114_cc00.obs_table['OBS_ID']} \n")

print("Ampiezza in gradi RA:", np.max(dst0114_cc00.obs_table['RA_PNT']) - np.min(dst0114_cc00.obs_table['RA_PNT']), '\n', "DEC:", np.max(dst0114_cc00.obs_table['DEC_PNT']) - np.min(dst0114_cc00.obs_table['DEC_PNT']))


#Presa delle osservazioni e reflected regions

observations = dst0114_cc00.get_observations(obs_id=dst0114_cc00.obs_table['OBS_ID'])

width = 3
binsz = 0.02
radius= np.sqrt(0.0225)
#dati copiati dall'esempio, non ho trovato le corrispondenze nei file aperti con fv

on_region_radius = Angle(radius * u.deg)
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)

npix = int(width/binsz)

exclusion_region = [CircleSkyRegion(
    center=target_position,
    radius=radius * u.deg,
    )]

for i, src in enumerate([]):
    src_position = SkyCoord(ra= [166.60, 166.90][i], 
                            dec= [38.200, 38.218][i], 
                            unit="deg", 
                            frame= ['icrs', 'icrs'][i])
    
    src_exclusion_region = CircleSkyRegion(
                                    center=src_position,
                                    radius=[0.05, 0.05][i] * u.deg,
                                    )
    exclusion_region.append(src_exclusion_region)

geom_excl = WcsGeom.create(
    npix=(npix, npix), binsz=binsz, skydir=target_position.galactic, proj="TAN", frame="galactic"
)

exclusion_mask = ~geom_excl.region_mask(exclusion_region)
exclusion_mask.plot()

plt.savefig('ProvaMaschera1D.png')

print(dst0114_cc00.obs_table.columns)

'''
e_min, e_max = (np.min(dst0114_cc00.hdu_table['ENERGY']) * u.TeV, 
                np.max(dst0114_cc00.hdu_table['ENERGY']) * u.TeV) 
print(e_min, e_max)

energy_axis = MapAxis.from_energy_bounds(
                    e_min, e_max, 
                    nbin=16, 
                    per_decade=False, unit="TeV", name="energy") 
#numero di bin non trovato

# Reduced IRFs are defined in true energy (i.e. not measured energy).
energy_axis_true = MapAxis.from_energy_bounds(
                    0.01, 
                    10000, 
                    nbin=300, 
                    per_decade=False, unit="TeV", name="energy_true" 
                    )
'''
'''
# geometry defining the ON region and SpectrumDataset based on it
geom = RegionGeom.create(region=on_region, axes=[energy_axis])

dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

dataset_maker = SpectrumDatasetMaker(
    containment_correction=True, selection=["counts", "exposure", "edisp"]
)

if global_parameters['DATASET_SETTINGS']['full_reflected_regions'] == False:
    # tell the background maker to use the WobbleRegionsFinder, let us use 3 off by default
    region_finder = WobbleRegionsFinder(n_off_regions=3)
    maker_fov = ReflectedRegionsBackgroundMaker(region_finder=region_finder, exclusion_mask=exclusion_mask)
else:
    maker_fov = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)

offset_max = global_parameters['SAFEMASK_DATSET']['offset_max'] * u.deg
safe_mask_masker = SafeMaskMaker(
#    methods=["offset-max", "aeff-max"], aeff_percent=5, offset_max=offset_max
    methods=["aeff-default","offset-max"], offset_max=offset_max
 )
'''
#
