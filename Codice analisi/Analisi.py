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

dir_gen = parameters['INPUT_DATA']['path_feb']

d_gen = []
i_gen= []

d_gen,i_gen = modan.readfile(d_gen, i_gen, dir_gen)
#n_obs = 16
n_obs =28

dst_gen = modan.store(d_gen, i_gen)

dst_gen.hdu_table.write("hdu-index.fits.gz", overwrite=True)
dst_gen.obs_table.write("obs-index.fits.gz", overwrite=True)

obs_date = dst_gen.obs_table[0]["EVENTS_FILENAME"].split("/")[-1][0:8]

observations = dst_gen.get_observations(obs_id=dst_gen.obs_table['OBS_ID'])


e_min, e_max = (0.199526 * u.TeV, 316.228 * u.TeV) 


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

modan.Theta2_full(energy_axis, energy_axis_true, dst_gen, observations, e_min, e_max)

#e_max_cut = 3.1601 * u.TeV
#e_min_cut = 0.499 * u.TeV
e_max_cut = 1.2601 * u.TeV
e_min_cut = 0.499 * u.TeV


modan.Theta2_cut(energy_axis, energy_axis_true, dst_gen, observations, e_min_cut, e_max_cut)

print('\n', 'Fatta Theta2','\n')

ex_reg, on_reg, geom = modan.mask(dst_gen, observations)

datasets = modan.CountMap(e_min_cut, e_max_cut, observations, ex_reg, on_reg, geom)


print('\n', 'Fatta Countmap','\n')

#print(datasets.energy_ranges)

modan.excess(datasets)

print('\n', 'Fatta Excess','\n')

time_range_rbr = modan.run_by_run(datasets, n_obs)

time_range_fm = modan.forty_min(datasets, n_obs)

#time_range_nbn = modan.night_by_night_gen(datasets, n_obs)
time_range_nbn = modan.night_by_night_feb(datasets, n_obs)

lind, ampl = modan.lightcurve_fit(datasets, e_min, e_max)

modan.lightcurve(datasets, e_min_cut, e_max, lind, ampl, time_range_rbr, '_rbr')

modan.lightcurve(datasets, e_min_cut, e_max, lind, ampl, time_range_fm, '_fm')

modan.lightcurve(datasets, e_min_cut, e_max, lind, ampl, time_range_nbn, '_nbn')

print('\n', 'Fatta Lightcurve','\n')

modan.SED(datasets, e_min_cut, e_max_cut)

print('\n', 'Fatta SED','\n')



