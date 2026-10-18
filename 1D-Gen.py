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
from gammapy.estimators import FluxPointsEstimator, LightCurveEstimator
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
    create_crab_spectral_model,
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

from gammapy.irf import EffectiveAreaTable2D


from gammapy.maps.utils import edges_from_lo_hi

from gammapy.stats import WStatCountsStatistic, CashCountsStatistic
from matplotlib.patches import Patch

import glob, os, pathlib, itertools, warnings

#Caricare i file di dati e irf

dir0114_cc10 = './Analisi-Gammapy/Mrk421/data/2026-01-14/cc10/'
dir0115_cc10 = './Analisi-Gammapy/Mrk421/data/2026-01-15/cc10/'
dir0116_cc10 = './Analisi-Gammapy/Mrk421/data/2026-01-16/cc10/'

d0114_cc10 = []
for folder in dir0114_cc10 :
    d0114_cc10.append(glob.glob(os.path.join(dir0114_cc10, '*Mrk421*_0201_*.fits.gz')))

d0114_cc10 = tuple(itertools.chain.from_iterable(d0114_cc10))


i0114_cc10 = []
for folder in dir0114_cc10 :
    i0114_cc10.append(glob.glob(os.path.join(dir0114_cc10, f'*Mrk421*_0502_*.fits.gz')))

i0114_cc10 = tuple(itertools.chain.from_iterable(i0114_cc10))

d0115_cc10 = []
for folder in dir0115_cc10 :
    d0115_cc10.append(glob.glob(os.path.join(dir0115_cc10, '*Mrk421*_0201_*.fits.gz')))

d0115_cc10 = tuple(itertools.chain.from_iterable(d0115_cc10))


i0115_cc10 = []
for folder in dir0115_cc10 :
    i0115_cc10.append(glob.glob(os.path.join(dir0115_cc10, f'*Mrk421*_0502_*.fits.gz')))

i0115_cc10 = tuple(itertools.chain.from_iterable(i0115_cc10))

d0116_cc10 = []
for folder in dir0116_cc10 :
    d0116_cc10.append(glob.glob(os.path.join(dir0116_cc10, '*Mrk421*_0201_*.fits.gz')))

d0116_cc10 = tuple(itertools.chain.from_iterable(d0116_cc10))


i0116_cc10 = []
for folder in dir0116_cc10 :
    i0116_cc10.append(glob.glob(os.path.join(dir0116_cc10, f'*Mrk421*_0502_*.fits.gz')))

i0116_cc10 = tuple(itertools.chain.from_iterable(i0116_cc10))

print(len(d0114_cc10))

#Controllo corretta lettura

#print(f"Il numero di run effettuate è ", len(d0114_cc10)+len(d0115_cc10)+len(d0116_cc10)"; il numero di IRF è ", len(i0114_cc10)+len(i0115_cc10)+len(i0116_cc10),".")
#if (len(i0114_cc10)+len(i0115_cc10)+len(i0116_cc10)) != (len(d0114_cc10)+len(d0115_cc10)+len(d0116_cc10)):
    #warnings.warn("Il numero di run effettuate non coincide con quello delle IRF, controllare l'input e riprovare.")
#else:
    #print("Le run presenti corrispondono alla IRF.")

print(np.sort(d0114_cc10))

