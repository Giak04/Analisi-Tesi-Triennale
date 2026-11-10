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

#Leggo il parfile
import yaml

def read_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)

parfile = "./Analisi_parfile.yaml"
parameters = read_yaml(parfile)

#Importo il modulo
import modulo_analisi as modan

dir_gen = parameters['INPUT_DATA']['path_gen']

d_gen = []
i_gen= []

d_gen,i_gen = modan.readfile(d_gen, i_gen, dir_gen)
n_obs = 16

dst_gen = modan.store(d_gen, i_gen)

dst_gen.hdu_table.write("hdu-index.fits.gz", overwrite=True)
dst_gen.obs_table.write("obs-index.fits.gz", overwrite=True)

obs_date = dst_gen.obs_table[0]["EVENTS_FILENAME"].split("/")[-1][0:8]

observations = dst_gen.get_observations(obs_id=dst_gen.obs_table['OBS_ID'])

for i in range(n_obs):
    print('\n', observations[0].tstart.isot, '\n')


