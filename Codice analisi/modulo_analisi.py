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


def gaus(x, a, sigma):
    return a * np.exp(-(x - 0)**2 / (2 * sigma**2))

def gaus2(x, a1, x1, sigma1, a2, x2, sigma2):
    return gaus(x, a1, sigma1) + gaus(x, a2, sigma2)


def gaus3(x, a1, sigma1, a2, sigma2, a3, sigma3):
    return gaus(x, a1, sigma1) + gaus(x, a2, sigma2)+ gaus(x, a3, sigma3)
    
def gaus_cost(x, a, sigma, C):
    return a * np.exp(-(x - 0)**2 / (2 * sigma**2)) + C

def king(r, a, sigma, gamma):
    return a * 1 / (2 * np.pi * sigma**2) * (1 - 1/gamma) * (1 + r**2 / (2 * gamma * sigma**2))**(-gamma)


def read_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)

parfile = "./Analisi_parfile.yaml"
parameters = read_yaml(parfile)



#Lettura dei file
#data = nome della variabile dei dati
#ird = nome della variabile delle irf
#dirdata = path della directori contenente i file .fits.gz
def readfile(data, irf, dirdata):

    for folder in dirdata :
        data.append(glob.glob(os.path.join(dirdata, f'*Mrk421*_0201_*.fits.gz')))

    for folder in dirdata :
        irf.append(glob.glob(os.path.join(dirdata, f'*Mrk421*_0502_*.fits.gz')))

    data = tuple(itertools.chain.from_iterable(data))
    irf = tuple(itertools.chain.from_iterable(irf))
    print(f"Il numero di run effettuate è {len(data)}; il numero di IRF è , {len(irf)}.")
    if len(irf) != len(data):
        warnings.warn("Il numero di run effettuate non coincide con quello delle IRF, controllare l'input e riprovare.")
    else:
        print("Le run presenti corrispondono alla IRF.")
    return data, irf

def store(data, irf):
    dst = DataStore.from_events_files(np.sort(data), irfs_paths=np.sort(irf))
    return dst

##########################################################################################

def readtime(data, irf, dirdata):
    data.append(glob.glob(os.path.join(dirdata, f'*Mrk421*_0201_*.fits.gz')))

    irf.append(glob.glob(os.path.join(dirdata, f'*Mrk421*_0502_*.fits.gz')))

    data = tuple(itertools.chain.from_iterable(data))
    irf = tuple(itertools.chain.from_iterable(irf))
    print(f"Il numero di run effettuate è {len(data)}; il numero di IRF è , {len(irf)}.")
    if len(irf) != len(data):
        warnings.warn("Il numero di run effettuate non coincide con quello delle IRF, controllare l'input e riprovare.")

    dst = DataStore.from_events_files(np.sort(data), irfs_paths=np.sort(irf))
    obs_date = dst_gen.obs_table[0]["EVENTS_FILENAME"].split("/")[-1][0:8]

    observations = dst_gen.get_observations(obs_id=dst_gen.obs_table['OBS_ID'])
    return dst




##########################################################################################

def Theta2_full(energy_axis, energy_axis_true, dst, observations, e_min, e_max):

    min_rad, max_rad, nbin = (0, 
                            0.5, 
                            50)
#ho cambiato il 50 da grafico a grafico ad un 300 per un solo grafico
    rad_axis = MapAxis.from_bounds(min_rad, max_rad, nbin=nbin, unit="deg2", interp="lin")

    energy_axis = MapAxis.from_energy_bounds(
                        e_min, e_max, 
                        nbin=16, 
                        per_decade=False, unit="TeV", name="energy")

    energy_axis = energy_axis.downsample(2) 
    erange = energy_axis.edges

    exp_tot = (np.sum(dst.obs_table['LIVETIME']) * u.s).to('hr')

    def make_theta_squared_table_modified(
        observations, 
        theta_squared_axis, 
        position, 
        position_off=None, 
        angle_pos_off = Angle(180 * u.deg),
        energy_edges=None
    ):
        if not theta_squared_axis.edges.unit.is_equivalent("deg2"):
            raise ValueError("The theta2 axis should be equivalent to deg2")

        table = Table()

        table["theta2_min"] = theta_squared_axis.edges[:-1]
        table["theta2_max"] = theta_squared_axis.edges[1:]
        table["counts"] = 0
        table["counts_off"] = 0
        table["acceptance"] = 0.0
        table["acceptance_off"] = 0.0

        alpha_tot = np.zeros(len(table))
        livetime_tot = 0

        create_off = position_off is None

        if energy_edges is not None:
            if len(energy_edges) == 2:
                table.meta["Energy_filter"] = energy_edges
            else:
                raise ValueError(
                    f"Only  supports one energy interval but {len(energy_edges) - 1} passed."
                )

        for observation in observations:
            events = observation.events
            if energy_edges is not None:
                events = events.select_energy(energy_range=energy_edges)

            event_position = events.radec
            pointing = observation.get_pointing_icrs(observation.tmid)

            separation = position.separation(event_position)
            counts, _ = np.histogram(separation**2, theta_squared_axis.edges)
            table["counts"] += counts

            if create_off:
                # Estimate the position of the mirror position
                pos_angle = pointing.position_angle(position)
                sep_angle = pointing.separation(position)
                position_off = pointing.directional_offset_by(
                    pos_angle + angle_pos_off, sep_angle
                )

            # Angular distance of the events from the mirror position
            separation_off = position_off.separation(event_position)

            # Extract the ON and OFF theta2 distribution from the two positions.
            counts_off, _ = np.histogram(separation_off**2, theta_squared_axis.edges)
            table["counts_off"] += counts_off

            # Normalisation between ON and OFF is one
            acceptance = np.ones(theta_squared_axis.nbin)
            acceptance_off = np.ones(theta_squared_axis.nbin)

            table["acceptance"] += acceptance
            table["acceptance_off"] += acceptance_off
            alpha = acceptance / acceptance_off
            alpha_tot += alpha * observation.observation_live_time_duration.to_value("s")
            livetime_tot += observation.observation_live_time_duration.to_value("s")

        alpha_tot /= livetime_tot
        table["alpha"] = alpha_tot

        stat = WStatCountsStatistic(table["counts"], table["counts_off"], table["alpha"])
        table["excess"] = stat.n_sig
        table["sqrt_ts"] = stat.sqrt_ts
        table["excess_errn"] = stat.compute_errn()
        table["excess_errp"] = stat.compute_errp()

        table.meta["ON_RA"] = position.icrs.ra
        table.meta["ON_DEC"] = position.icrs.dec
        return table

    def make_theta_squared_table_reflected(
        observations, 
        theta_squared_axis, 
        position, 
        position_off=None,
        off_regions_number=1,
        energy_edges=None
    ):
        if not theta_squared_axis.edges.unit.is_equivalent("deg2"):
            raise ValueError("The theta2 axis should be equivalent to deg2")
        
        table = Table()
        
        table["theta2_min"] = theta_squared_axis.edges[:-1]
        table["theta2_max"] = theta_squared_axis.edges[1:]
        table["counts"] = 0
        table["counts_off"] = 0.0
        table["acceptance"] = 0.0
        table["acceptance_off"] = 0.0

        alpha_tot = np.zeros(len(table))
        livetime_tot = 0

        create_off = position_off is None

        if energy_edges is not None:
            if len(energy_edges) == 2:
                table.meta["Energy_filter"] = energy_edges
            else:
                raise ValueError(
                    f"Only  supports one energy interval but {len(energy_edges) - 1} passed."
                )

        for observation in observations:
            events = observation.events
            if energy_edges is not None:
                events = events.select_energy(energy_range=energy_edges)

            event_position = events.radec
            pointing = observation.get_pointing_icrs(observation.tmid)

            separation = position.separation(event_position)
            counts, _ = np.histogram(separation**2, theta_squared_axis.edges)
            table["counts"] += counts

            pos_angle = pointing.position_angle(position)
            sep_angle = pointing.separation(position)

            if off_regions_number > 1 and create_off is False:
                raise ValueError("If `off_regions_number` is larger than 1, `position_off` has to be set to None.")

            if off_regions_number > 1:
                off_region_angles = Angle(360 * u.deg) / (off_regions_number + 1)
                off_region_angle = pos_angle + off_region_angles * np.arange(1, off_regions_number + 1)
                off_regions = pointing.directional_offset_by(off_region_angle, sep_angle)

                counts_off = np.zeros_like(counts)
                for off_region in off_regions:
                    # Angular distance of the events from the mirror position
                    separation_off = off_region.separation(event_position)
     
                    # Extract the ON and OFF theta2 distribution from the two positions.
                    counts_off_regions, _ = np.histogram(separation_off**2, theta_squared_axis.edges)
                    counts_off += counts_off_regions

            if off_regions_number == 1: 
                if create_off:
                    # Estimate the position of the mirror position
                    position_off = pointing.directional_offset_by(
                        pos_angle + Angle(np.pi, "rad"), sep_angle
                )
    
                separation_off = position_off.separation(event_position)
                counts_off, _ = np.histogram(separation_off**2, theta_squared_axis.edges)

            table["counts_off"] += np.array(counts_off)
    
            # Normalisation between ON and OFF is one
            acceptance = np.ones(theta_squared_axis.nbin)
            acceptance_off = np.ones(theta_squared_axis.nbin) * off_regions_number

            table["acceptance"] += acceptance
            table["acceptance_off"] += acceptance_off
            alpha = acceptance / acceptance_off
            alpha_tot += alpha * observation.observation_live_time_duration.to_value("s")
            livetime_tot += observation.observation_live_time_duration.to_value("s")

        alpha_tot /= livetime_tot
        table["alpha"] = alpha_tot

        stat = WStatCountsStatistic(table["counts"], table["counts_off"], table["alpha"])
        table["excess"] = stat.n_sig
        table["sqrt_ts"] = stat.sqrt_ts
        table["excess_errn"] = stat.compute_errn()
        table["excess_errp"] = stat.compute_errp()

        table.meta["ON_RA"] = position.icrs.ra
        table.meta["ON_DEC"] = position.icrs.dec
        return table


    n, sigma, c = 50, 0.02, 0
    col = 4
    row = int((len(erange) - 1) / col)
    ratio = int(col / row)

    fig, axs = plt.subplots(row, col, figsize=(18, 18/ratio))

    j, jj = 0, 0
    best, best_err = [], []
    
    '''
    print('erange', erange)
    #qui devo ancora capire come sistemare erange, per ora lo ridefinisco e basta
    #erange=np.array([1.99526000e-01, 3.16228000e+02])*u.TeV
    #Ho modificato il range guardando quali range di energia aveva senso osservare e quali no
    erange=np.array([0.5, 3.16228000e+02])*u.TeV
    '''
    target_position = SkyCoord(ra=parameters['TARGET']['ra'], 
                            dec=parameters['TARGET']['dec'], 
                            unit="deg", 
                            frame='icrs')
    
    for i in np.arange(len(erange)-1):
        energy = np.sqrt(erange[1:] * erange[:-1])[i]


        #print("Default OFF position(s)")
        theta2_dist = make_theta_squared_table_reflected(observations=observations, 
                                               theta_squared_axis=rad_axis, 
                                               position=target_position, 
                                               position_off=None,
                                               off_regions_number = 1,
                                               energy_edges=[erange[i], erange[i+1]])
        #print(theta2_dist, " energy: ", energy)
    #off_regions_number messo ad 1 perchè era il valore che avevo visto prima

        theta2_dist_mod = theta2_dist[theta2_dist["theta2_max"] <= 0.04]

        alpha = np.sum(theta2_dist_mod["alpha"]) / len(theta2_dist_mod["alpha"])
        n_on = theta2_dist_mod["counts"].sum()
        n_off = theta2_dist_mod["counts_off"].sum() 

        stat_ws = WStatCountsStatistic(n_on=n_on, n_off=n_off, alpha=alpha)
        excs = stat_ws.n_sig
        excess_errn = stat_ws.compute_errn()
        excess_errp = stat_ws.compute_errp()


        theta2_edges = edges_from_lo_hi(
        theta2_dist["theta2_min"].quantity, theta2_dist["theta2_max"].quantity
                )
        theta2_axis = MapAxis.from_edges(theta2_edges, interp="lin", name="theta_squared")

        x = theta2_axis.center.value
        x_edges = theta2_axis.edges.value
        xerr = (x - x_edges[:-1], x_edges[1:] - x)

        axs[j,jj].errorbar(
            x,
            theta2_dist["counts"],
            xerr=xerr,
            yerr=np.sqrt(theta2_dist["counts"]),
            linestyle="None",
            label=r"N$_{on}$ = "+ f"{n_on}" + r"$\pm$" + f"{int(n_on**0.5)}",
        )

        axs[j,jj].errorbar(
            x,
            theta2_dist["counts_off"] * theta2_dist["alpha"],
            xerr=xerr,
            yerr=np.sqrt(theta2_dist["counts_off"] * theta2_dist["alpha"]),
            linestyle="None",
            label=r"N$_{off}$ = "+ f"{np.round(n_off * alpha,2)}" + r"$\pm$" + f"{int((n_off * alpha)**0.5)}; " + r"$\alpha$ = " + f"{np.round(alpha,2)}",
        )

    
        axs[j,jj].errorbar(
            x,
            theta2_dist["excess"],
            xerr=xerr,
            yerr=(theta2_dist["excess_errn"], theta2_dist["excess_errp"]),
            fmt="+",
            linestyle="None",
            label=r"N$_{exc}$ = "+ f"{np.round(excs,2)}" + r"$\pm$" + f"{int(excess_errn)}",
        )

        axs[j,jj].set_ylabel("Counts", fontsize=14)
        axs[j,jj].set_xlabel("Theta2 (deg2)", fontsize=14)
        axs[j,jj].legend()

        if False:
            popt, pcov = curve_fit(gaus_cost, x, theta2_dist["excess"],
                                    sigma=theta2_dist["excess_errn"],
                                    p0=[n, sigma, c],
                                    maxfev = 100000000)
    
            best.append(popt)
            best_err.append(pcov)
        
            sig_irf = irf_psf.containment_radius(fraction=0.68, energy_true=[energy]*u.TeV, offset=[0.5]*u.deg) * 0.6624305
            axs[j,jj].plot(x, gaus_cost(x, *popt), label=f"{np.round(popt[1]**0.5,3)} deg, {np.round(sig_irf[0],3)}", color="green")

        y_min, y_max = axs[j,jj].get_ylim()
        axs[j,jj].fill_between(x_edges[(x_edges <= 0.04)], y_min, y_max, alpha=0.3)

        axs[j,jj].set_title(f"{np.round(erange[i],2)}-{np.round(erange[i+1],2)}; " + r"$\sigma$ = " + f"{np.round(stat_ws.sqrt_ts,2)}", fontsize="16")

        legend = axs[j,jj].legend(fontsize="6", loc ="upper right")
        handles, labels = axs[j,jj].get_legend_handles_labels()
        handles.append(Patch(facecolor='w', edgecolor='w'))
        labels.append(f"Time = {np.round(exp_tot,2)} \n" + 
                  r"$N_{exc}/Time$ = " + f"{np.round((excs) / exp_tot.value,2)} " + r"evt/h" + "\n" +
                  r"$N_{off}/Time$ = " + f"{np.round((n_off * alpha) / exp_tot.value,2)} " + r"evt/h" + "\n" +
                  f"Significance (Li&Ma) = {np.round(stat_ws.sqrt_ts,2)}" + r"$\sigma$"
                 )
    
        legend._legend_box = None
        legend._init_legend_box(handles, labels)
        legend._set_loc(legend._loc)
        legend.set_title(legend.get_title().get_text())

        axs[j,jj].set_ylim(y_min,y_max)
        axs[j,jj].set_xlim(0)

        axs[j,jj].hlines(0, x_edges[0], x_edges[-1], color="green", linestyles="dashed")
        jj += 1
        if jj == 4:
            jj = 0
            j += 1

        
    fig.tight_layout(pad=2.0)

    plt.savefig(f"./Graphs/Theta2_plot_full.png")
    plt.close()

########################################################################################

def Theta2_cut(energy_axis, energy_axis_true, dst, observations, e_min, e_max):
    energy_axis = energy_axis.downsample(2) 
    erange = energy_axis.edges

    min_rad, max_rad, nbin = (0, 
                            0.5, 
                            50)

    rad_axis = MapAxis.from_bounds(min_rad, max_rad, nbin=nbin, unit="deg2", interp="lin")

    energy_axis = MapAxis.from_energy_bounds(
                        e_min, e_max, 
                        nbin=2, 
                        per_decade=False, unit="TeV", name="energy")

    energy_axis = energy_axis.downsample(2) 
    erange = energy_axis.edges    
    #erange = np.array([e_min, e_max]) * u.TeV

    #ho cambiato il 50 da grafico a grafico ad un 300 per un solo grafico

    exp_tot = (np.sum(dst.obs_table['LIVETIME']) * u.s).to('hr')

    def make_theta_squared_table_modified(
        observations, 
        theta_squared_axis, 
        position, 
        position_off=None, 
        angle_pos_off = Angle(180 * u.deg),
        energy_edges=None
    ):
        if not theta_squared_axis.edges.unit.is_equivalent("deg2"):
            raise ValueError("The theta2 axis should be equivalent to deg2")

        table = Table()

        table["theta2_min"] = theta_squared_axis.edges[:-1]
        table["theta2_max"] = theta_squared_axis.edges[1:]
        table["counts"] = 0
        table["counts_off"] = 0
        table["acceptance"] = 0.0
        table["acceptance_off"] = 0.0

        alpha_tot = np.zeros(len(table))
        livetime_tot = 0
    
        create_off = position_off is None

        if energy_edges is not None:
            if len(energy_edges) == 2:
                table.meta["Energy_filter"] = energy_edges
            else:
                raise ValueError(
                    f"Only  supports one energy interval but {len(energy_edges) - 1} passed."
                )

        for observation in observations:
            events = observation.events
            if energy_edges is not None:
                events = events.select_energy(energy_range=energy_edges)

            event_position = events.radec
            pointing = observation.get_pointing_icrs(observation.tmid)

            separation = position.separation(event_position)
            counts, _ = np.histogram(separation**2, theta_squared_axis.edges)
            table["counts"] += counts

            if create_off:
                # Estimate the position of the mirror position
                pos_angle = pointing.position_angle(position)
                sep_angle = pointing.separation(position)
                position_off = pointing.directional_offset_by(
                    pos_angle + angle_pos_off, sep_angle
                )
    
            # Angular distance of the events from the mirror position
            separation_off = position_off.separation(event_position)

            # Extract the ON and OFF theta2 distribution from the two positions.
            counts_off, _ = np.histogram(separation_off**2, theta_squared_axis.edges)
            table["counts_off"] += counts_off

            # Normalisation between ON and OFF is one
            acceptance = np.ones(theta_squared_axis.nbin)
            acceptance_off = np.ones(theta_squared_axis.nbin)

            table["acceptance"] += acceptance
            table["acceptance_off"] += acceptance_off
            alpha = acceptance / acceptance_off
            alpha_tot += alpha * observation.observation_live_time_duration.to_value("s")
            livetime_tot += observation.observation_live_time_duration.to_value("s")

        alpha_tot /= livetime_tot
        table["alpha"] = alpha_tot

        stat = WStatCountsStatistic(table["counts"], table["counts_off"], table["alpha"])
        table["excess"] = stat.n_sig
        table["sqrt_ts"] = stat.sqrt_ts
        table["excess_errn"] = stat.compute_errn()
        table["excess_errp"] = stat.compute_errp()

        table.meta["ON_RA"] = position.icrs.ra
        table.meta["ON_DEC"] = position.icrs.dec
        return table

    def make_theta_squared_table_reflected(
        observations, 
        theta_squared_axis, 
        position, 
        position_off=None,
        off_regions_number=1,
        energy_edges=None
    ):
        if not theta_squared_axis.edges.unit.is_equivalent("deg2"):
            raise ValueError("The theta2 axis should be equivalent to deg2")

        table = Table()

        table["theta2_min"] = theta_squared_axis.edges[:-1]
        table["theta2_max"] = theta_squared_axis.edges[1:]
        table["counts"] = 0
        table["counts_off"] = 0.0
        table["acceptance"] = 0.0
        table["acceptance_off"] = 0.0

        alpha_tot = np.zeros(len(table))
        livetime_tot = 0

        create_off = position_off is None
    
        if energy_edges is not None:
            if len(energy_edges) == 2:
                table.meta["Energy_filter"] = energy_edges
            else:
                raise ValueError(
                    f"Only  supports one energy interval but {len(energy_edges) - 1} passed."
                )

        for observation in observations:
            events = observation.events
            if energy_edges is not None:
                events = events.select_energy(energy_range=energy_edges)

            event_position = events.radec
            pointing = observation.get_pointing_icrs(observation.tmid)

            separation = position.separation(event_position)
            counts, _ = np.histogram(separation**2, theta_squared_axis.edges)
            table["counts"] += counts

            pos_angle = pointing.position_angle(position)
            sep_angle = pointing.separation(position)

            if off_regions_number > 1 and create_off is False:
                raise ValueError("If `off_regions_number` is larger than 1, `position_off` has to be set to None.")

            if off_regions_number > 1:
                off_region_angles = Angle(360 * u.deg) / (off_regions_number + 1)
                off_region_angle = pos_angle + off_region_angles * np.arange(1, off_regions_number + 1)
                off_regions = pointing.directional_offset_by(off_region_angle, sep_angle)
    
                counts_off = np.zeros_like(counts)
                for off_region in off_regions:
                    # Angular distance of the events from the mirror position
                    separation_off = off_region.separation(event_position)
     
                    # Extract the ON and OFF theta2 distribution from the two positions.
                    counts_off_regions, _ = np.histogram(separation_off**2, theta_squared_axis.edges)
                    counts_off += counts_off_regions

            if off_regions_number == 1: 
                if create_off:
                    # Estimate the position of the mirror position
                    position_off = pointing.directional_offset_by(
                        pos_angle + Angle(np.pi, "rad"), sep_angle
                    )
    
                separation_off = position_off.separation(event_position)
                counts_off, _ = np.histogram(separation_off**2, theta_squared_axis.edges)

            table["counts_off"] += np.array(counts_off)

            # Normalisation between ON and OFF is one
            acceptance = np.ones(theta_squared_axis.nbin)
            acceptance_off = np.ones(theta_squared_axis.nbin) * off_regions_number

            table["acceptance"] += acceptance
            table["acceptance_off"] += acceptance_off
            alpha = acceptance / acceptance_off
            alpha_tot += alpha * observation.observation_live_time_duration.to_value("s")
            livetime_tot += observation.observation_live_time_duration.to_value("s")

        alpha_tot /= livetime_tot
        table["alpha"] = alpha_tot

        stat = WStatCountsStatistic(table["counts"], table["counts_off"], table["alpha"])
        table["excess"] = stat.n_sig
        table["sqrt_ts"] = stat.sqrt_ts
        table["excess_errn"] = stat.compute_errn()
        table["excess_errp"] = stat.compute_errp()

        table.meta["ON_RA"] = position.icrs.ra
        table.meta["ON_DEC"] = position.icrs.dec
        return table


    fig, axs = plt.subplots(1, 1, figsize=(18, 18))

    best, best_err = [], []

    #qui devo ancora capire come sistemare erange, per ora lo ridefinisco e basta
    #erange=np.array([1.99526000e-01, 3.16228000e+02])*u.TeV
    #Ho modificato il range guardando quali range di energia aveva senso osservare e quali no

    target_position = SkyCoord(ra=parameters['TARGET']['ra'], 
                            dec=parameters['TARGET']['dec'], 
                            unit="deg", 
                            frame='icrs')

    for i in np.arange(len(erange)-1):
        energy = np.sqrt(erange[1:] * erange[:-1])[i]


        theta2_dist = make_theta_squared_table_reflected(observations=observations, 
                                               theta_squared_axis=rad_axis, 
                                               position=target_position, 
                                               position_off=None,
                                               off_regions_number = 1,
                                               energy_edges=[erange[i], erange[i+1]])


        theta2_dist_mod = theta2_dist[theta2_dist["theta2_max"] <= 0.04]

        alpha = np.sum(theta2_dist_mod["alpha"]) / len(theta2_dist_mod["alpha"])
        n_on = theta2_dist_mod["counts"].sum()
        n_off = theta2_dist_mod["counts_off"].sum() 

        stat_ws = WStatCountsStatistic(n_on=n_on, n_off=n_off, alpha=alpha)
        excs = stat_ws.n_sig
        excess_errn = stat_ws.compute_errn()
        excess_errp = stat_ws.compute_errp()
    

        theta2_edges = edges_from_lo_hi(
        theta2_dist["theta2_min"].quantity, theta2_dist["theta2_max"].quantity
                )
        theta2_axis = MapAxis.from_edges(theta2_edges, interp="lin", name="theta_squared")

        x = theta2_axis.center.value
        x_edges = theta2_axis.edges.value
        xerr = (x - x_edges[:-1], x_edges[1:] - x)

        axs.errorbar(
            x,
            theta2_dist["counts"],
            xerr=xerr,
            yerr=np.sqrt(theta2_dist["counts"]),
            linestyle="None",
            label=r"N$_{on}$ = "+ f"{n_on}" + r"$\pm$" + f"{int(n_on**0.5)}",
        )

        axs.errorbar(
            x,
            theta2_dist["counts_off"] * theta2_dist["alpha"],
            xerr=xerr,
            yerr=np.sqrt(theta2_dist["counts_off"] * theta2_dist["alpha"]),
            linestyle="None",
            label=r"N$_{off}$ = "+ f"{np.round(n_off * alpha,2)}" + r"$\pm$" + f"{int((n_off * alpha)**0.5)}; " + r"$\alpha$ = " + f"{np.round(alpha,2)}",
        )

    
        axs.errorbar(
            x,
            theta2_dist["excess"],
            xerr=xerr,
            yerr=(theta2_dist["excess_errn"], theta2_dist["excess_errp"]),
            fmt="+",
            linestyle="None",
            label=r"N$_{exc}$ = "+ f"{np.round(excs,2)}" + r"$\pm$" + f"{int(excess_errn)}",
        )

        axs.set_ylabel("Counts", fontsize=14)
        axs.set_xlabel("Theta2 (deg2)", fontsize=14)
        axs.legend()

        if False:
            popt, pcov = curve_fit(gaus_cost, x, theta2_dist["excess"],
                                    sigma=theta2_dist["excess_errn"],
                                    p0=[n, sigma, c],
                                    maxfev = 100000000)
    
            best.append(popt)
            best_err.append(pcov)
        
            sig_irf = irf_psf.containment_radius(fraction=0.68, energy_true=[energy]*u.TeV, offset=[0.5]*u.deg) * 0.6624305
            axs.plot(x, gaus_cost(x, *popt), label=f"{np.round(popt[1]**0.5,3)} deg, {np.round(sig_irf[0],3)}", color="green")

        y_min, y_max = axs.get_ylim()
        axs.fill_between(x_edges[(x_edges <= 0.04)], y_min, y_max, alpha=0.3)

        axs.set_title(f"{np.round(erange[i],2)}-{np.round(erange[i+1],2)}; " + r"$\sigma$ = " + f"{np.round(stat_ws.sqrt_ts,2)}", fontsize="22")

        legend = axs.legend(fontsize="18", loc ="upper right")
        handles, labels = axs.get_legend_handles_labels()
        handles.append(Patch(facecolor='w', edgecolor='w'))
        labels.append(f"Time = {np.round(exp_tot,2)} \n" + 
                  r"$N_{exc}/Time$ = " + f"{np.round((excs) / exp_tot.value,2)} " + r"evt/h" + "\n" +
                      r"$N_{off}/Time$ = " + f"{np.round((n_off * alpha) / exp_tot.value,2)} " + r"evt/h" + "\n" +
                      f"Significance (Li&Ma) = {np.round(stat_ws.sqrt_ts,2)}" + r"$\sigma$"
                     )
    
        legend._legend_box = None
        legend._init_legend_box(handles, labels)
        legend._set_loc(legend._loc)
        legend.set_title(legend.get_title().get_text())

        axs.set_ylim(y_min,y_max)
        axs.set_xlim(0)

        axs.hlines(0, x_edges[0], x_edges[-1], color="green", linestyles="dashed")

    plt.savefig(f"./Graphs/Theta2_cut_plot.png")
    plt.close()

############################################################################################

def mask(dst, obs):
    ev_header = Table.read(dst.obs_table["EVENTS_FILENAME"][0], hdu=1)
    meta = ev_header.meta
    loc = EarthLocation(
                lat=parameters['DATA_ORIGIN']['lat'],
                lon=parameters['DATA_ORIGIN']['lon'],
            )
    mjdtime = ( ((dst.obs_table["TSTART"] + dst.obs_table["TSTOP"])/2).to("d") 
           + (dst.obs_table["MJDREFI"] + dst.obs_table["MJDREFF"]) * u.d )
    observing_time = Time(mjdtime, format="mjd", scale="tt") 
    
    altazframe = AltAz(location=loc, obstime=observing_time)
    dst.obs_table["DATE-OBS"] = observing_time.isot 
    coord = SkyCoord(dst.obs_table["RA_PNT"], dst.obs_table["DEC_PNT"], unit="deg", frame="icrs")
    coord_altaz = coord.transform_to(altazframe)

    dst.obs_table["AZ_PNT"] = coord_altaz.az.deg * u.deg
    dst.obs_table["ZD_PNT"] = coord_altaz.zen.deg * u.deg

    target_position = SkyCoord(ra=parameters['TARGET']['ra'], 
                            dec=parameters['TARGET']['dec'], 
                            unit="deg", 
                            frame='icrs')

    selection_offset = dict(type='sky_circle',
                    lon=target_position.ra.deg,
                    lat=target_position.dec.deg,
                    frame="icrs",
                    radius=Angle(5, 'deg')
                    )

    print(f"The mean exposure of the runs is {np.mean(dst.obs_table['LIVETIME']) * u.s} \n")
    print(f"The total exposure time is {(np.sum(dst.obs_table['LIVETIME']) * u.s).to('hr')} \n")
    print(f"The total runs are {len(dst.obs_table['OBS_ID'])}. These are the selected ObsID: \n {dst.obs_table['OBS_ID']} \n")

    print("Ampiezza in gradi RA:", np.max(dst.obs_table['RA_PNT']) - np.min(dst.obs_table['RA_PNT']), '\n', "DEC:", np.max(dst.obs_table['DEC_PNT']) - np.min(dst.obs_table['DEC_PNT']))

    width = 3
    binsz = 0.02
    radius= np.sqrt(0.0225)

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
        npix=(npix, npix), binsz=binsz, skydir=target_position.galactic, proj="TAN", frame="galactic")

    exclusion_mask = ~geom_excl.region_mask(exclusion_region)

    return exclusion_mask, on_region, geom_excl

########################################################################################

def CountMap(e_min, e_max, observations, exclusion_mask, on_region, geom_excl):
    energy_axis = MapAxis.from_energy_bounds(
                    e_min, e_max, 
                    nbin=16, 
                    per_decade=False, unit="TeV", name="energy") 

    energy_axis_true = MapAxis.from_energy_bounds(
                        0.01, 
                        10000, 
                        nbin=300, 
                        per_decade=False, unit="TeV", name="energy_true" 
                        )

    geom = RegionGeom.create(region=on_region, axes=[energy_axis])

    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )

    maker_fov = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)

    offset_max = 5 * u.deg
    safe_mask_masker = SafeMaskMaker(
    #    methods=["offset-max", "aeff-max"], aeff_percent=5, offset_max=offset_max
        methods=["aeff-default","offset-max"], offset_max=offset_max
    )

    datasets = Datasets()

    counts= Map.from_geom(geom_excl)

    for observation in observations:
        dataset = dataset_maker.run(
            dataset_empty.copy(name=str(observation.obs_id)), observation
        )
        counts.fill_events(observation.events)
        dataset_on_off = maker_fov.run(dataset, observation)
        dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)

    datasets=datasets.slice_by_energy(e_min, e_max)

    ax = counts.smooth(0.05 * u.deg).plot(stretch="linear", cmap="viridis", add_cbar=True)
    geom.plot_region(ax=ax, kwargs_point={"color": "k", "marker": "*"}, color="r", linewidth=2.5)
    plot_spectrum_datasets_off_regions(ax=ax, datasets=datasets, linewidth=2.5)

    plt.savefig(f"./Graphs/1D_sp_an_countmap.png", bbox_inches='tight')
    plt.close()

    return datasets


##########################################################################################

def excess(datasets):
    
    
        
    info_table= datasets.info_table(cumulative=True)

    fig, (ax_excess, ax_sqrt_ts) = plt.subplots(figsize=(10, 4), ncols=2, nrows=1)
    ax_excess.plot(
        info_table["livetime"].to("h"),
        info_table["excess"],
        marker="o",
        ls="none",
    )

    ax_excess.set_title("Excess")
    ax_excess.set_xlabel("Livetime [h]")
    ax_excess.set_ylabel("Excess events")

    ax_sqrt_ts.plot(
        info_table["livetime"].to("h"),
        info_table["sqrt_ts"],
        marker="o",
        ls="none",
    )

    ax_sqrt_ts.set_title("Sqrt(TS)")
    ax_sqrt_ts.set_xlabel("Livetime [h]")
    ax_sqrt_ts.set_ylabel("Sqrt(TS)")
    plt.savefig(f"./Graphs/Excess-SqrtTS.png")
    plt.close()

###########################################################################################

def lightcurve_fit(datasets, emin, emax):

    spectral_model = PowerLawSpectralModel(index=parameters['LIGHTCURVE']['index'], amplitude=parameters['LIGHTCURVE']['amplitude'], reference="1 TeV")



    flux = spectral_model.energy_flux_error(emin,emax).to("erg cm-2 s-1")
    flux_int = spectral_model.integral_error(emin,emax)

    dataset_stacked = Datasets(datasets).stack_reduce()

    model = SkyModel(spectral_model=spectral_model, name=parameters['TARGET']['name'])
    dataset_stacked.models = model
    dataset_stacked.models[0].spectral_model.amplitude.min = 0

    dataset_stacked.mask_fit = dataset_stacked.counts.geom.energy_mask(
        energy_min=emin, energy_max=emax
    )


    sky_model = SkyModel(spectral_model=dataset_stacked.models[0].spectral_model, name=parameters['TARGET']['name'])
    
    models = Models([sky_model])
    datasets.models = models

    fit = Fit()
    result_stacked = fit.run(dataset_stacked)

    #print(result_stacked)

    dataset_stacked.models.parameters.to_table()

    print(dataset_stacked.models[0].parameters[1].value)
    
    index = dataset_stacked.models[0].parameters[0].value
    amplitude = dataset_stacked.models[0].parameters[1].value

    return index, amplitude

##########################################################################################

def run_by_run(datasets, n_obs):
    
    time_range = []
    for i in range(n_obs-1):
            time_range.append(Time([str((datasets[i].gti.time_start - 5*u.s).isot[0]), str((datasets[i].gti.time_stop + 5*u.s).isot[-1])], 
                         format = "isot",
                         scale = datasets[i].gti.time_start.scale
                        ))
        #time_range.append([observations[i].tstart.isot, observations[i].tstop.isot])
        #i += 1

    return time_range

###########################################################################################

def forty_min(datasets, n_obs):
    
    time_range = []
    for i in range(0, n_obs-1, 2):
        time_range.append(Time([str((datasets[i].gti.time_start - 5*u.s).isot[0]), str((datasets[i].gti.time_stop + 5*u.s).isot[-1])], format = "isot", scale = datasets[i].gti.time_start.scale))
        #time_range.append([observations[i].tstart.isot, observations[i+1].tstop.isot])
        #print('Controllo \n', observations[i].tstart.isot, observations[i+1].tstop.isot, '\n', i)


    return time_range

###########################################################################################

def night_by_night_gen(datasets, n_obs):

    time_range = []
    time_range.append(Time([str((datasets[0].gti.time_start - 5*u.s).isot[0]), str((datasets[3].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[0].gti.time_start.scale))
    time_range.append(Time([str((datasets[4].gti.time_start - 5*u.s).isot[0]), str((datasets[9].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[4].gti.time_start.scale))
    time_range.append(Time([str((datasets[10].gti.time_start - 5*u.s).isot[0]), str((datasets[-1].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[10].gti.time_start.scale))
    #time_range.append([observations[0].tstart.isot, observations[3].tstop.isot])
    #time_range.append([observations[4].tstart.isot, observations[9].tstop.isot])
    #time_range.append([observations[10].tstart.isot, observations[-1].tstop.isot])
    print('Controllo \n', time_range)
    return time_range

###########################################################################################

def night_by_night_feb(datasets, n_obs):

    time_range = []
    time_range.append(Time([str((datasets[0].gti.time_start - 5*u.s).isot[0]), str((datasets[6].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[0].gti.time_start.scale))
    time_range.append(Time([str((datasets[7].gti.time_start - 5*u.s).isot[0]), str((datasets[13].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[4].gti.time_start.scale))
    time_range.append(Time([str((datasets[14].gti.time_start - 5*u.s).isot[0]), str((datasets[19].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[10].gti.time_start.scale))
    time_range.append(Time([str((datasets[20].gti.time_start - 5*u.s).isot[0]), str((datasets[-1].gti.time_stop + 5*u.s).isot[-1])], format = "isot",scale = datasets[10].gti.time_start.scale))
    #print('Controllo \n', time_range)
    return time_range

###########################################################################################

def lightcurve(datasets, emin, emax, index, amplitude, time_range, tipo):

    spectral_model = PowerLawSpectralModel(index=parameters['LIGHTCURVE']['index'], amplitude=parameters['LIGHTCURVE']['amplitude'], reference="1 TeV")



    flux = spectral_model.energy_flux_error(emin,emax).to("erg cm-2 s-1")
    flux_int = spectral_model.integral_error(emin,emax)

    dataset_stacked = Datasets(datasets).stack_reduce()

    model = SkyModel(spectral_model=spectral_model, name=parameters['TARGET']['name'])
    dataset_stacked.models = model
    dataset_stacked.models[0].spectral_model.amplitude.min = 0

    dataset_stacked.mask_fit = dataset_stacked.counts.geom.energy_mask(
        energy_min=emin, energy_max= emax
    )

    sky_model = SkyModel(spectral_model=dataset_stacked.models[0].spectral_model, name=parameters['TARGET']['name'])
    
    models = Models([sky_model])
    datasets.models = models


    '''time_range = [['2026-01-15T02:01:04', '2026-01-15T02:39:24'],['2026-01-15T02:41:06','2026-01-15T03:19:26'],['2026-01-15T03:24:16','2026-01-15T04:02:36'],['2026-01-15T04:04:42','2026-01-15T04:43:02'], ['2026-01-16T02:16:24','2026-01-16T02:54:44'], ['2026-01-16T05:01:04','2026-01-16T05:39:25'], ['2026-01-16T02:58:34','2026-01-16T03:36:54'], ['2026-01-16T05:41:43','2026-01-16T06:00:16'], ['2026-01-16T03:39:31','2026-01-16T04:17:51'], ['2026-01-16T04:19:38','2026-01-16T04:57:58'], ['2026-01-17T02:03:21','2026-01-17T02:40:01'], ['2026-01-17T04:39:30','2026-01-17T05:16:40'], ['2026-01-17T02:42:03','2026-01-17T03:18:43'], ['2026-01-17T05:18:24','2026-01-17T05:55:34'], ['2026-01-17T03:21:18','2026-01-17T03:58:28'], ['2026-01-17T04:00:24','2026-01-17T04:37:34']]
'''

    #if time_range is not None:
    #Time_intervals = []

    '''
        for i in range (len(time_range)):
            Time_intervals.append(Time([str((datasets[i].gti.time_start - 5*u.s).isot[0]), str((datasets[i].gti.time_stop + 5*u.s).isot[-1])], 
                         format = "isot",
                         scale = datasets[i].gti.time_start.scale
                        ))


        Time_intervals.append(Time(time_range, 
                         format = "isot",
                         scale = datasets[0].gti.time_start.scale))
    else:
        raise ValueError("For user defined time intervals, one amongst either `time_usr_defined` or `bin_time` keywords has to be defined!")
    '''



    lc_maker_all = LightCurveEstimator(
        energy_edges=[emin,emax],
        time_intervals=time_range,
        source="Mrk421",
        reoptimize=False,
        selection_optional="all",
        n_sigma=1, 
        n_sigma_ul=3, 
        n_jobs=7
    )

    lc_maker_all.norm.scan_max = 3

    time_lc = lc_maker_all.run(datasets)
    time_lc.sqrt_ts_threshold_ul = 1.5 


    fig, ax = plt.subplots(
        figsize=(8, 6),
        gridspec_kw={"left": 0.16, "bottom": 0.2, "top": 0.98, "right": 0.98},
    )

    sed_type = "eflux"

    b = time_lc.plot(ax=ax, sed_type=sed_type, color="tab:orange", label = f"{emin} - {emax}", axis_name="time")

#ax.set_yscale("linear")
    ax.set_yscale("log")

    if sed_type == "eflux":
        ax.plot(ax.xaxis.get_view_interval(), 
                np.ones_like(ax.xaxis.get_view_interval()) * flux[0], label = parameters['TARGET']['name'], ls="dashed")

    if sed_type == "flux":
        ax.plot(ax.xaxis.get_view_interval(), 
                np.ones_like(ax.xaxis.get_view_interval()) * flux_int[0], label = parameters['TARGET']['name'], ls="dashed")

#time_lc.plot_ts_profiles(ax=ax, sed_type=sed_type)
    plt.legend()
    plt.savefig(f"./Graphs/LC_Gen{tipo}.png", bbox_inches='tight')
    plt.close()

##########################################################################################

def SED(datasets, e_min, e_max):

    energy_axis = MapAxis.from_energy_bounds(
                    e_min, e_max, 
                    nbin=16, 
                    per_decade=False, unit="TeV", name="energy") 

    energy_edges = energy_axis.downsample(1).edges 

    spectral_model = PowerLawSpectralModel(index=parameters['LIGHTCURVE']['index'], amplitude=parameters['LIGHTCURVE']['amplitude'], reference="1 TeV")

    model = SkyModel(spectral_model=spectral_model, name="Mrk421")

    dataset_stacked = Datasets(datasets).stack_reduce()
    dataset_stacked.models = model
    
    dataset_stacked.mask_fit = dataset_stacked.counts.geom.energy_mask(
        energy_min=e_min, energy_max=e_max
    )

    sky_model = SkyModel(spectral_model=dataset_stacked.models[0].spectral_model, name='Mrk421')


    models = Models([sky_model])


    fpe = FluxPointsEstimator(
        energy_edges=energy_edges, 
        source="Mrk421", 
        selection_optional="all", 
        n_sigma=1, 
        n_sigma_ul=3, 
        n_jobs=7
    )


    flux_points = fpe.run(datasets=dataset_stacked)
    flux_points.sqrt_ts_threshold_ul = 1.5

    flux_points.to_table()

    fig, ax = plt.subplots()
    sed_type = "e2dnde" 
    flux_points.plot(ax=ax, sed_type=sed_type, color="darkorange")
    flux_points.plot_ts_profiles(ax=ax, sed_type=sed_type)
    model.spectral_model.plot(ax=ax, energy_bounds=(e_min, e_max), sed_type=sed_type)
    model.spectral_model.plot_error(ax=ax, energy_bounds=(e_min, e_max), sed_type=sed_type)
    plt.savefig(f"./Graphs/SED.png")
    plt.close()






#End
