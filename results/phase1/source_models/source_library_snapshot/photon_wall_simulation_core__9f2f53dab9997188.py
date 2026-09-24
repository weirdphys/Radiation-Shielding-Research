# 2026-09-01-v7.3k.12e14b-clipped-grid-stage-hotfix: installed
# 2026-08-17-v7.3k.8-support-aware-sparse-tail-provenance: immutable coarse/fine ledger + boundary provenance
# 2026-08-17-v7.3k.7-inline-vr-single-native-run: cumulative v7.3k.6 + full-spectrum/fine-N + inline VR
# 2026-08-11-v7.2.3-final-correctness: installed
# 2026-08-11-v7.2.1-runtime-hotfix-spectrum-trim-import
# 2026-08-11-v7.2-correctness-report-neutron-cleanup: installed
# 2026-08-11-photon-core-v7-all-photon-air-kerma: installed
try:
    get_ipython().run_line_magic('matplotlib', 'inline')
except NameError:
    pass

import hashlib
import os
import math
import time
from copy import deepcopy
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap
from geant4_pybind import *
from mpl_toolkits.mplot3d import Axes3D


PLOT_BLACK = '#000000'
PLOT_RED = '#ff2020'
PLOT_ORANGE = '#ff7a00'
PLOT_YELLOW = '#ffd400'
PLOT_DARK_RED = '#5a0000'

PLOT_LINE_WIDTH = 0.35
PLOT_TRACK_LINE_WIDTH = 0.10
PLOT_EDGE_WIDTH = 0.08
PLOT_ERROR_LINE_WIDTH = 0.25
PLOT_MARKER_EDGE_WIDTH = 0.20
PLOT_UNCERTAINTY_ALPHA = 0.10

RED_ORANGE_YELLOW_CMAP = LinearSegmentedColormap.from_list(
    'black_red_orange_yellow',
    [
        (0.00, PLOT_BLACK),
        (0.18, PLOT_DARK_RED),
        (0.52, PLOT_RED),
        (0.78, PLOT_ORANGE),
        (1.00, PLOT_YELLOW),
    ],
)

mpl.rcParams.update(
    {
        'figure.facecolor': PLOT_BLACK,
        'figure.edgecolor': PLOT_BLACK,
        'axes.facecolor': PLOT_BLACK,
        'axes.edgecolor': PLOT_RED,
        'axes.linewidth': 0.5,
        'axes.labelcolor': PLOT_RED,
        'axes.titlecolor': PLOT_RED,
        'axes.prop_cycle': cycler(color=[PLOT_RED, PLOT_ORANGE]),
        'text.color': PLOT_RED,
        'xtick.color': PLOT_RED,
        'ytick.color': PLOT_RED,
        'xtick.labelcolor': PLOT_RED,
        'ytick.labelcolor': PLOT_RED,
        'lines.color': PLOT_RED,
        'lines.linewidth': PLOT_LINE_WIDTH,
        'lines.markerfacecolor': PLOT_RED,
        'lines.markeredgecolor': PLOT_RED,
        'lines.markeredgewidth': PLOT_MARKER_EDGE_WIDTH,
        'patch.facecolor': PLOT_RED,
        'patch.edgecolor': PLOT_RED,
        'patch.linewidth': PLOT_EDGE_WIDTH,
        'patch.force_edgecolor': True,
        'legend.facecolor': PLOT_BLACK,
        'legend.edgecolor': PLOT_RED,
        'grid.color': PLOT_RED,
        'grid.alpha': 0.22,
        'savefig.facecolor': PLOT_BLACK,
        'savefig.edgecolor': PLOT_BLACK,
        'savefig.transparent': False,
    }
)


def alternating_red_orange_colors(count: int) -> List[str]:
    return [
        PLOT_RED if i % 2 == 0 else PLOT_ORANGE
        for i in range(max(0, int(count)))
    ]


def style_black_red_figure(fig) -> None:
    fig.patch.set_facecolor(PLOT_BLACK)
    fig.patch.set_edgecolor(PLOT_BLACK)

    for ax in fig.axes:
        ax.set_facecolor(PLOT_BLACK)
        ax.tick_params(axis='both', which='both', colors=PLOT_RED)

        ax.xaxis.label.set_color(PLOT_RED)
        ax.yaxis.label.set_color(PLOT_RED)
        ax.title.set_color(PLOT_RED)

        ax.xaxis.get_offset_text().set_color(PLOT_RED)
        ax.yaxis.get_offset_text().set_color(PLOT_RED)

        for spine in ax.spines.values():
            spine.set_color(PLOT_RED)

        for text_artist in ax.texts:
            text_artist.set_color(PLOT_RED)

        legend = ax.get_legend()

        if legend is not None:
            legend.get_frame().set_facecolor(PLOT_BLACK)
            legend.get_frame().set_edgecolor(PLOT_RED)

            for legend_text in legend.get_texts():
                legend_text.set_color(PLOT_RED)

        if hasattr(ax, 'zaxis'):
            ax.zaxis.label.set_color(PLOT_RED)
            ax.zaxis.get_offset_text().set_color(PLOT_RED)
            ax.tick_params(axis='z', which='both', colors=PLOT_RED)

            for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                try:
                    axis.pane.set_facecolor(PLOT_BLACK)
                    axis.pane.set_edgecolor(PLOT_RED)
                    axis.pane.set_alpha(1.0)
                except AttributeError:
                    pass

                try:
                    axis._axinfo['grid']['color'] = (
                        1.0,
                        0.125,
                        0.125,
                        0.22,
                    )
                    axis._axinfo['axisline']['color'] = PLOT_RED
                except (AttributeError, KeyError, TypeError):
                    pass


CFG = SimpleNamespace(
    beam_half_angle_deg=11.3,
    source_mode='digitized_spectrum',
    spectrum_case='15MV_40x40_interpolated',
    mono_energy_divisor=1.6,
    production_cut_mm=0.01,
    plot_dir='plots',
    plot_dpi=100,
    mev=0.0,
    wall_thickness_cm=50.0,
    num_photons=10_000,
    transport_threads=None,
    n_visual_tracks=1000,
    concrete_density_g_cm3=None,
    wall_size_y_cm=1000.0,
    wall_size_z_cm=1000.0,
    source_to_wall_gap_cm=100.0,
    depth_nbins=600,
    max_scatter_bin=120,
    energy_hist_bins=120,
    angle_hist_bins=120,
    corr_energy_bins=80,
    corr_angle_bins=80,
    secondary_energy_hist_bins=120,
    max_raw_exit_samples=300_000,
    export_obj=True,
    obj_filepath='photon_tracks_and_wall.obj',
    obj_scale=0.01,
    obj_decimate=1,
    obj_track_radius_cm=0.1,
    include_wall_in_obj=True,
    min_pts_after_decimate=3,
    plot_track_decimate=5,
    run_thickness_sweep=False,
    sweep_thicknesses_cm=[50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300],
    sweep_num_photons_per_point=1_000_000,
    geant4_control_verbose=0,
    geant4_run_verbose=0,
    geant4_event_verbose=0,
    geant4_tracking_verbose=0,
    verbose=True,
    # User-selected design target for the transmitted fraction. Required
    # concrete thickness is derived from this run's Geant4 transmission.
    target_transmission=1.0e-5,
    mu_ratio_E_tab=None,
    mu_ratio_values=None,
)

_LIVE_OBJECTS: List[Any] = []
MEV_TO_J = 1.602176634e-13

PRIMARY_GAMMA_PROCESSES = ('compt', 'phot', 'conv', 'Rayl', 'photonNuclear')
SECONDARY_PARTICLES_TO_TALLY = ('gamma', 'e-', 'e+')


SIMULATION_THICKNESS_FEATURE_BUILD = '2026-08-10-photon-concrete-sweep-v3'


def _suffix_maximum(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float).ravel()
    if values.size == 0:
        return values.copy()
    output = values.copy()
    for index in range(output.size - 2, -1, -1):
        output[index] = max(output[index], output[index + 1])
    return output


def _log_interpolate_target_thickness(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    target: float,
) -> float:
    """Interpolate thickness linearly in log-transmission between two simulated points."""
    if not (x1 > x0):
        return math.nan
    if not (y0 > 0.0 and y1 > 0.0 and target > 0.0):
        return math.nan
    log0 = math.log(y0)
    log1 = math.log(y1)
    logt = math.log(target)
    if abs(log1 - log0) <= 1.0e-15:
        return x1
    fraction = (logt - log0) / (log1 - log0)
    fraction = min(1.0, max(0.0, fraction))
    return x0 + fraction * (x1 - x0)


def calculate_sweep_derived_concrete_requirement(
    sweep_results: Dict[str, np.ndarray],
    cfg,
) -> Dict[str, Any]:
    """Derive required concrete thickness only from the Geant4 thickness sweep.

    Each thickness is a separate Geant4 transport run. The design curve uses
    the two-sided 95% Wilson upper bound on the transmitted-photon fraction at
    every simulated thickness. A suffix maximum makes the design curve
    conservatively non-increasing, so a compliant point is accepted only when
    all thicker sampled points are also compliant. The reported crossing is
    log-linearly interpolated only between two simulated points. No TVL and no
    extrapolation beyond the simulated thickness range are used.
    """
    target = float(getattr(cfg, 'target_transmission', 1.0e-5))
    if not (0.0 < target <= 1.0):
        raise ValueError('TARGET_TRANSMISSION must be greater than 0 and no greater than 1.')

    thicknesses = np.asarray(sweep_results.get('thicknesses_cm', []), dtype=float).ravel()
    counts = np.asarray(sweep_results.get('transmitted_count', []), dtype=np.int64).ravel()
    totals = np.asarray(sweep_results.get('num_photons', []), dtype=np.int64).ravel()
    observed = np.asarray(sweep_results.get('transmitted_fraction', []), dtype=float).ravel()
    wilson_low = np.asarray(sweep_results.get('transmitted_fraction_wilson95_low', []), dtype=float).ravel()
    wilson_high = np.asarray(sweep_results.get('transmitted_fraction_wilson95_high', []), dtype=float).ravel()

    n = thicknesses.size
    if n < 2:
        raise ValueError('Concrete thickness sweep requires at least two simulated thicknesses.')
    if any(array.size != n for array in (counts, totals, observed, wilson_low, wilson_high)):
        raise ValueError('Concrete sweep result arrays have inconsistent lengths.')
    if np.any(~np.isfinite(thicknesses)) or np.any(np.diff(thicknesses) <= 0.0):
        raise ValueError('Concrete sweep thicknesses must be finite and strictly increasing.')
    if np.any(totals <= 0) or np.any(counts < 0) or np.any(counts > totals):
        raise ValueError('Concrete sweep contains invalid photon counts.')

    design95 = _suffix_maximum(wilson_high)
    design_central = _suffix_maximum(observed)

    compliant = np.flatnonzero(design95 <= target)
    first_index = int(compliant[0]) if compliant.size else -1

    required_95 = math.nan
    required_central = math.nan
    relation = f'> {thicknesses[-1]:.6g} cm'
    first_compliant = math.nan
    previous_noncompliant = math.nan
    bracket = None

    if first_index == 0:
        required_95 = float(thicknesses[0])
        relation = f'<= {thicknesses[0]:.6g} cm (target already satisfied at minimum simulated thickness)'
        first_compliant = float(thicknesses[0])
    elif first_index > 0:
        previous_noncompliant = float(thicknesses[first_index - 1])
        first_compliant = float(thicknesses[first_index])
        required_95 = _log_interpolate_target_thickness(
            float(thicknesses[first_index - 1]),
            float(design95[first_index - 1]),
            float(thicknesses[first_index]),
            float(design95[first_index]),
            target,
        )
        relation = 'log-interpolated between adjacent simulated thicknesses'
        bracket = [previous_noncompliant, first_compliant]

    central_compliant = np.flatnonzero(design_central <= target)
    if central_compliant.size:
        ci = int(central_compliant[0])
        if ci == 0:
            required_central = float(thicknesses[0])
        else:
            required_central = _log_interpolate_target_thickness(
                float(thicknesses[ci - 1]),
                float(max(design_central[ci - 1], np.finfo(float).tiny)),
                float(thicknesses[ci]),
                float(max(design_central[ci], np.finfo(float).tiny)),
                target,
            )

    # Exact two-sided Wilson upper limit for zero transmitted events at the
    # configured per-point history count. This is a useful statistical-floor
    # diagnostic for very small target transmissions.
    n_per_point = int(np.min(totals))
    zero_low, zero_high = wilson_interval(0, n_per_point)
    z2 = 1.959963984540054 ** 2
    min_histories_zero = int(math.ceil(z2 * (1.0 - target) / target))

    if first_index < 0:
        status = 'TARGET NOT REACHED WITHIN SIMULATED CONCRETE RANGE; INCREASE MAX_CONCRETE_CM'
    elif first_index == 0:
        status = 'TARGET MET AT MINIMUM SIMULATED CONCRETE THICKNESS'
    else:
        status = 'REQUIRED CONCRETE THICKNESS BRACKETED BY GEANT4 SWEEP'

    return {
        'feature_build': SIMULATION_THICKNESS_FEATURE_BUILD,
        'method': 'multi-thickness Geant4 transmission sweep with no TVL and no out-of-range extrapolation',
        'statistical_basis': '95% Wilson upper bound on transmitted fraction; suffix-maximum conservative design curve',
        'target_transmission': target,
        'sweep_point_count': int(n),
        'sweep_min_thickness_cm': float(thicknesses[0]),
        'sweep_max_thickness_cm': float(thicknesses[-1]),
        'sweep_step_nominal_cm': float(np.median(np.diff(thicknesses))),
        'sweep_num_photons_per_point': n_per_point,
        'zero_count_wilson95_upper_at_sweep_n': float(zero_high),
        'minimum_histories_for_zero_count_to_bound_target_95': min_histories_zero,
        'thicknesses_cm': thicknesses.tolist(),
        'num_primary_photons_by_thickness': totals.astype(int).tolist(),
        'transmitted_photons_by_thickness': counts.astype(int).tolist(),
        'observed_transmission_by_thickness': observed.tolist(),
        'wilson95_low_by_thickness': wilson_low.tolist(),
        'wilson95_high_by_thickness': wilson_high.tolist(),
        'conservative_design_transmission95_by_thickness': design95.tolist(),
        'central_design_transmission_by_thickness': design_central.tolist(),
        'required_thickness_central_cm': float(required_central),
        'required_thickness_95_cm': float(required_95),
        'required_thickness_relation': relation,
        'previous_noncompliant_simulated_thickness_cm': float(previous_noncompliant),
        'first_compliant_simulated_thickness_cm': float(first_compliant),
        'interpolation_bracket_cm': bracket,
        'selected_validation_thickness_cm': float(first_compliant if math.isfinite(first_compliant) else thicknesses[-1]),
        'status': status,
        'no_extrapolation': True,
        'note': (
            'Required concrete thickness is derived directly from separate Geant4 simulations at multiple concrete '
            'thicknesses. The 95% upper transmitted-fraction curve is conservatively forced non-increasing with '
            'thickness by a suffix maximum. The crossing is interpolated only inside a simulated bracket. '
            'If the target is not reached, the maximum sweep thickness must be increased.'
        ),
    }


def attach_final_concrete_validation(
    requirement: Dict[str, Any],
    results: Dict[str, Any],
) -> Dict[str, Any]:
    """Attach the final high-statistics validation run to a sweep requirement."""
    output = dict(requirement)
    summary = results.get('summary', {})
    total = int(summary.get('num_photons', 0))
    transmitted = int(summary.get('transmitted_count', 0))
    thickness = float(summary.get('wall_thickness_cm', math.nan))
    target = float(output.get('target_transmission', math.nan))
    if total <= 0 or transmitted < 0 or transmitted > total:
        raise ValueError('Invalid final transmitted-photon counts.')
    observed = transmitted / total
    low, high = wilson_interval(transmitted, total)
    meets = bool(math.isfinite(target) and high <= target)
    output.update({
        'final_validation_thickness_cm': thickness,
        'final_validation_primary_photons': total,
        'final_validation_transmitted_photons': transmitted,
        'final_validation_observed_transmission': observed,
        'final_validation_wilson95_low': float(low),
        'final_validation_wilson95_high': float(high),
        'final_validation_meets_target_95': meets,
    })
    if math.isfinite(float(output.get('required_thickness_95_cm', math.nan))):
        if meets:
            output['status'] = 'GEANT4 SWEEP REQUIREMENT FOUND; FINAL HIGH-STATISTICS RUN MEETS TARGET AT 95% UPPER BOUND'
        else:
            output['status'] = 'SWEEP FOUND A COMPLIANT GRID POINT, BUT FINAL HIGH-STATISTICS VALIDATION EXCEEDS TARGET'
    else:
        output['status'] = 'TARGET NOT REACHED WITHIN SIMULATED CONCRETE RANGE; FINAL RUN USED MAX_CONCRETE_CM'
    return output


def _energy_centers(max_mev: float) -> List[float]:
    n = int(
        math.floor(
            (float(max_mev) - 0.05) / 0.1 + 1e-12
        )
    ) + 1

    return [
        round(0.05 + 0.1 * i, 2)
        for i in range(n)
    ]


def _coerce_to_grid(
    arr: Sequence[float],
    target_grid: Sequence[float],
) -> List[float]:
    arr = list(arr)
    n_target = len(target_grid)

    if len(arr) < n_target:
        arr = arr + [0.0] * (n_target - len(arr))
    elif len(arr) > n_target:
        arr = arr[:n_target]

    return arr


def _make_spectrum(
    case: str,
    beam_mv: int,
    field_cm: str,
    peak_mev: float,
    energy_mev_center: Sequence[float],
    probability_mass_bin: Sequence[float],
    origin: str,
    interpolation_method: Optional[str] = None,
) -> Dict[str, Any]:
    energy_mev_center = list(energy_mev_center)
    probability_mass_bin = list(probability_mass_bin)

    if len(energy_mev_center) != len(probability_mass_bin):
        raise ValueError(
            f'{case}: energy grid length '
            f'{len(energy_mev_center)} != probability length '
            f'{len(probability_mass_bin)}'
        )

    spec = {
        'case': case,
        'beam_mv': int(beam_mv),
        'field_cm': str(field_cm),
        'bin_width_mev': 0.1,
        'peak_mev_paper_text': float(peak_mev),
        'energy_mev_center': energy_mev_center,
        'probability_mass_bin': probability_mass_bin,
        'origin': str(origin),
    }

    if interpolation_method is not None:
        spec['interpolation_method'] = str(interpolation_method)

    return spec


E_7_05 = _energy_centers(7.05)
E_18_05 = _energy_centers(18.05)

E_40x40_COMMON = np.round(
    np.arange(0.1, 19.0 + 0.0001, 0.1),
    1,
).tolist()

E_7_plot = np.round(
    np.arange(0.1, 7.0 + 0.0001, 0.1),
    1,
).tolist()

E_19_plot = np.round(
    np.arange(0.1, 19.0 + 0.0001, 0.1),
    1,
).tolist()


def _pick_grid_for_array(arr):
    n = len(arr)

    if n == len(E_7_plot):
        return E_7_plot

    if n == len(E_19_plot):
        return E_19_plot

    if n == len(E_7_05):
        return E_7_05

    if n == len(E_18_05):
        return E_18_05

    raise ValueError(f'Unsupported spectrum length {n}')


P_3MV_40x40 = [
    0.029362826661446607,
    0.03621415288245081,
    0.041107957326025256,
    0.04208671821474013,
    0.041107957326025256,
    0.03963981599295292,
    0.038661055104238035,
    0.03817167465988059,
    0.03797592248213761,
    0.03787804639326612,
    0.037829108348830374,
    0.037780170304394636,
    0.03773123225995889,
    0.037682294215523146,
    0.037584418126651656,
    0.03748654203778017,
    0.03729078986003719,
    0.03699716159342273,
    0.036605657237936774,
    0.03611627679357933,
    0.03543114417147891,
    0.034452383282764024,
    0.033082118038563174,
    0.031320348438876384,
    0.02906919839483214,
    0.02623079181755897,
    0.022511500440442398,
    0.017617695996867962,
    0.011255750220221199,
    0.00371929137711657,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
] 

P_6MV_40x40 = [
    0.00019968051118210862,
    0.00039936102236421724,
    0.000998402555910543,
    0.001996805111821086,
    0.003993610223642172,
    0.007987220447284345,
    0.014976038338658145,
    0.02496006389776358,
    0.03394568690095847,
    0.03194888178913738,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.031449680511182104,
    0.030950479233226837,
    0.02995207667731629,
    0.028953674121405752,
    0.027955271565495207,
    0.026956869009584664,
    0.02496006389776358,
    0.022963258785942492,
    0.020966453674121407,
    0.01896964856230032,
    0.016972843450479235,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]

P_10MV_40x40 = [
    0.015208379515996531,
    0.03223584531393477,
    0.04500445448211188,
    0.04944051807032153,
    0.048711045551301284,
    0.04549980023371701,
    0.04153152137010203,
    0.03653502270027422,
    0.029653454385722873,
    0.025903016985611906,
    0.023567506181062042,
    0.021733033468112974,
    0.020074895522926402,
    0.01871939876175759,
    0.017617144723790553,
    0.016657307096760804,
    0.016491210349127212,
    0.01640657421195152,
    0.016453159680711253,
    0.01656877406830686,
    0.016782075190503953,
    0.0174887115778214,
    0.01811374481420654,
    0.01861466612837256,
    0.01882318486078525,
    0.018925685174813023,
    0.01902211994997879,
    0.01908877237858356,
    0.01911303980494982,
    0.019193150533736034,
    0.019286993898087718,
    0.019335519536216712,
    0.01936093373961756,
    0.01935556687846909,
    0.019308522638542883,
    0.019242636636575973,
    0.01906523114647181,
    0.018930059873052928,
    0.01876721060556861,
    0.018570809957836647,
    0.018360372522258733,
    0.01811873977460741,
    0.01786432074425621,
    0.017594845465545354,
    0.017309087753596398,
    0.017017023071570523,
    0.016674071506235066,
    0.01636337228136356,
    0.0160459947642318,
    0.01571790230058058,
    0.015385182617229234,
    0.01504367841310873,
    0.014686863252496019,
    0.014315832244668453,
    0.013930687814840319,
    0.013529408733710321,
    0.01269716281208644,
    0.011032972915266588,
    0.007872031389754279,
    0.0037272402854149167,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]

P_15MV_40x40 = [
    0.00951344954280765,
    0.016710314917218406,
    0.022971150382536654,
    0.028227084172875987,
    0.030102040365155402,
    0.030754361651027323,
    0.030362125109743482,
    0.029004853246207614,
    0.027455599725059884,
    0.025726207680378976,
    0.02361437678922386,
    0.02099262743610313,
    0.018295794634098278,
    0.01671062668344399,
    0.015532123160354057,
    0.014587345100413127,
    0.013811982911770639,
    0.013090282346743877,
    0.012420149784848452,
    0.011847232603786603,
    0.01135566287006342,
    0.010906636131491381,
    0.010500949092386393,
    0.010363649437239492,
    0.010304464329553183,
    0.010268691837492893,
    0.010317661762264159,
    0.010448148923327056,
    0.010622089477816154,
    0.010810941255700496,
    0.011030010989475732,
    0.011340288380930279,
    0.011617601090407834,
    0.011868068899603257,
    0.012088018737274995,
    0.012184713000408893,
    0.012253926803228806,
    0.012298933727059239,
    0.012341277329824348,
    0.012376820815494636,
    0.012405342835767122,
    0.012415998434439053,
    0.012440385383674,
    0.012476043640897704,
    0.012517249379547669,
    0.012546515965786528,
    0.012566605355450703,
    0.012577764493132595,
    0.012580273598216368,
    0.012574989036759247,
    0.012554332353478786,
    0.012528049866660184,
    0.012489326511211912,
    0.012411429429081456,
    0.012346884493018297,
    0.012284615302954539,
    0.01221310972126898,
    0.012130407747002124,
    0.012042444311788926,
    0.01195004321308048,
    0.011846684059871057,
    0.011738788805365006,
    0.011627075802057263,
    0.011509809510832317,
    0.011388911643695298,
    0.011263438120051645,
    0.011135527584915882,
    0.01099834714929184,
    0.010847760342558752,
    0.01071020229922831,
    0.010572486946705442,
    0.01039585413146584,
    0.010156251705221367,
    0.009895803785780206,
    0.009659768102945212,
    0.009479706050088844,
    0.00922506483438884,
    0.0088003147157128,
    0.008315362559644934,
    0.007901861546167075,
    0.007656557790366106,
    0.007507304967210684,
    0.007353980314706763,
    0.007197750991751175,
    0.007040127473409324,
    0.006880617416778187,
    0.0067178458819528485,
    0.00654352276999353,
    0.006352502575468079,
    0.00614350647316775,
    0.005894036542245355,
    0.005602325945632791,
    0.005284021593159442,
    0.004933522278023205,
    0.00455220014864509,
    0.004162163443420539,
    0.0037630679529008794,
    0.0033585832651606,
    0.0029726946899836177,
    0.002623187933673389,
    0.0023197762867113924,
    0.002047506218478392,
    0.0017990909534174742,
    0.0015798790740448215,
    0.0013919751483565694,
    0.0012352933275207375,
    0.0011106478313576663,
    0.0010086796779015896,
    0.0009167359664743734,
    0.0008260243529778099,
    0.0007310422788810781,
    0.0006319321972841428,
    0.0005344470049494491,
    0.00044290877627310655,
    0.000360220726181998,
    0.00028796332590214,
    0.00021393334719623725,
    0.00014668023981953394,
    0.00010006132006052461,
    7.32170130994748e-05,
    6.195246945685248e-05,
    5.856687014996861e-05,
    5.644005134000426e-05,
    5.485779521489391e-05,
    5.3208183840862214e-05,
    5.095434611663834e-05,
    4.783023273493817e-05,
    4.4181366639063134e-05,
    4.040420076398639e-05,
    3.6804417782030605e-05,
    3.3596557037457774e-05,
    3.069135554064832e-05,
    2.7941888763568026e-05,
    2.539974874138158e-05,
    2.309684359776491e-05,
    2.1049460928438712e-05,
    1.9242360888725777e-05,
    1.7621804952796665e-05,
    1.6145639297213084e-05,
    1.4781247463008413e-05,
    1.350380259186705e-05,
    1.2354603250182816e-05,
    1.1328815210752619e-05,
    1.0355863084949528e-05,
    9.38636363401059e-06,
    8.390365282154897e-06,
    7.349030506887416e-06,
    6.305896806735839e-06,
    5.321786557565569e-06,
    4.434942406757977e-06,
    3.6640771691949466e-06,
    2.977152603028851e-06,
    2.36722386373954e-06,
    1.8577944674429228e-06,
    1.45142451023347e-06,
    1.1385393480092267e-06,
    9.030919689976586e-07,
    7.2020193930921e-07,
    5.713273896474897e-07,
    4.460610852938863e-07,
    3.391190943689592e-07,
    2.398137217950027e-07,
    1.5809376374075601e-07,
    1.0330090545350468e-07,
    7.11331355450687e-08,
    5.4884279593230566e-08,
    4.701005435960954e-08,
    4.187851512015701e-08,
    3.793830003168341e-08,
    3.417266626352104e-08,
    2.992420387580164e-08,
    2.5390746547944727e-08,
    2.088978661007888e-08,
    1.6391616214873074e-08,
    1.178442616589583e-08,
    7.671265780858619e-09,
    4.640622094961494e-09,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]

P_16MV_40x40 = [
    0.010717704249746699,
    0.018443832866893474,
    0.02533259195297145,
    0.03129927189856643,
    0.034254209413129656,
    0.03552492318541492,
    0.03568080642248269,
    0.03480668722175438,
    0.033373088182649946,
    0.03172109091277494,
    0.02976372893468661,
    0.027329267037289016,
    0.024458652009797983,
    0.02174675897141974,
    0.02007868197778266,
    0.01881102831210025,
    0.017777401257282586,
    0.01692404402289769,
    0.016117036578762255,
    0.015355348384823202,
    0.014681707300804663,
    0.014093005241915018,
    0.01355923984825422,
    0.013070204813160059,
    0.012764007916318914,
    0.012674719538958244,
    0.012666043997290244,
    0.012648737864947343,
    0.012623268632257765,
    0.012590103789549733,
    0.012549710827151471,
    0.012502557235391204,
    0.012449110504597153,
    0.01238946053675552,
    0.012324524332125508,
    0.01225476938103534,
    0.012180663173813246,
    0.012102673200787447,
    0.012021266952286166,
    0.011936911918637626,
    0.011850075590170052,
    0.011761225457211666,
    0.011670829010090696,
    0.011579353739135363,
    0.01148726713467389,
    0.011395036687034505,
    0.011303129886545426,
    0.011212014223534881,
    0.011122157188331091,
    0.011034188094837437,
    0.010948484531395943,
    0.010865513988334836,
    0.010785743955982335,
    0.010709641924666672,
    0.010637675384716064,
    0.010570311826458734,
    0.010508018740222912,
    0.010451263616336815,
    0.010400513945128672,
    0.010356237216926705,
    0.010318900922059138,
    0.010288972550854193,
    0.010266919593640097,
    0.010253209540745073,
    0.010248309882497341,
    0.009965823999564305,
    0.009678165186581017,
    0.009385878080581991,
    0.00909166386754827,
    0.008796673537749015,
    0.008503138413665205,
    0.008209820484539106,
    0.007914896369701973,
    0.007626633994675811,
    0.007350472513217651,
    0.00707415811599635,
    0.006780418321877046,
    0.006466770867035058,
    0.006157342984693832,
    0.005874766607515145,
    0.005633310425377189,
    0.005362561870426925,
    0.005019953872639046,
    0.004660088070471712,
    0.004341544043595854,
    0.004105033992548348,
    0.003930922220140474,
    0.0037774263705132867,
    0.003627806799985297,
    0.003482431001802897,
    0.003341342975222799,
    0.003204036476578603,
    0.0030690545358617732,
    0.00293354082427308,
    0.002795946293824692,
    0.0026554436106295187,
    0.002504471903879964,
    0.002342377541301879,
    0.002176223598743008,
    0.0020045356943617417,
    0.0018280132317908138,
    0.0016529584098966626,
    0.0014809518595902157,
    0.001312703498583544,
    0.001152905396340636,
    0.0010072430334465707,
    0.0008790540975137985,
    0.0007678195923990851,
    0.0006695894147726198,
    0.0005821635441655876,
    0.0005062210345875118,
    0.00044164856981854684,
    0.0003878317636471774,
    0.00034514675448375147,
    0.00030981688551203195,
    0.0002784282487861484,
    0.00024863729727348454,
    0.00021897848479900442,
    0.00018920094306554218,
    0.00016024670172063083,
    0.00013319358564776703,
    0.0001087596417778645,
    8.733850607709295e-05,
    6.80487560242064e-05,
    4.9601695496373467e-05,
    3.3952324598245705e-05,
    2.3194281358906295e-05,
    1.6498990008233436e-05,
    1.3515895418747572e-05,
    1.2276668672322647e-05,
    1.1315915765599805e-05,
    1.0497146857324773e-05,
    9.718111279486984e-06,
    8.9038211083398e-06,
    8.022561294198244e-06,
    7.110565508364304e-06,
    6.2255859332821355e-06,
    5.411382315915137e-06,
    4.728058658465195e-06,
    4.188486611823056e-06,
    3.7092759085701346e-06,
    3.2764284795298913e-06,
    2.89158805188676e-06,
    2.554151543690333e-06,
    2.261965471970279e-06,
    2.0126896610198776e-06,
    1.7953237687336743e-06,
    1.604030936143155e-06,
    1.4342063712030246e-06,
    1.228983931322315e-06,
    1.0423546398449129e-06,
    8.884078424577591e-07,
    7.552519027936282e-07,
    6.373399121815499e-07,
    5.312881214305309e-07,
    4.351837441513834e-07,
    3.4859877931733227e-07,
    2.729390128447621e-07,
    2.0933651144993076e-07,
    1.583809368058914e-07,
    1.2922625245293637e-07,
    1.0408622221088629e-07,
    8.226638381435719e-08,
    6.421156386625434e-08,
    4.981493563712492e-08,
    3.8664263898893435e-08,
    3.0201021366060366e-08,
    2.3723366545194895e-08,
    1.861198157172918e-08,
    1.4456632487117212e-08,
    1.102046927118915e-08,
    8.104815955602544e-09,
    5.586159596732693e-09,
    3.6379415428270186e-09,
    2.3586639400956762e-09,
    1.6065567125818368e-09,
    1.2120660932488522e-09,
    1.0112160673471842e-09,
    8.715960335628771e-10,
    7.627853593390796e-10,
    6.658264530509658e-10,
    5.694362378784134e-10,
    4.740046204543355e-10,
    3.837838157493178e-10,
    3.001338289528116e-10,
    2.2189955563366337e-10,
    1.5143482358985935e-10,
    9.512068547083345e-11,
    4.741345434950861e-11,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]

P_18MV_40x40 = [
    0.008179741124264233,
    0.014498387433257168,
    0.02010899680781405,
    0.0249271924263894,
    0.026811313457140693,
    0.02762147263246124,
    0.027490422043684175,
    0.026467214970175044,
    0.025242017578371162,
    0.023822081897352996,
    0.02201580263609074,
    0.01969746026431407,
    0.017270282661619403,
    0.015861819407102063,
    0.014818306441500591,
    0.013980892615180755,
    0.013291592321640832,
    0.01264130176169754,
    0.012029260136354731,
    0.011500952725999122,
    0.011466449867821126,
    0.011461811756585848,
    0.01145717364535057,
    0.011452535534115292,
    0.011447897422880014,
    0.011443259311644735,
    0.011438621200409458,
    0.011433983089174178,
    0.011429344977938901,
    0.011424706866703624,
    0.011420068755468344,
    0.011415430644233067,
    0.01141079253299779,
    0.011406154421762512,
    0.011401516310527235,
    0.011396878199291955,
    0.011392240088056678,
    0.0113876019768214,
    0.011382963865586121,
    0.011378325754350844,
    0.011373687643115566,
    0.011369049531880289,
    0.011364411420645011,
    0.011359773309409732,
    0.011355135198174455,
    0.011350497086939177,
    0.011345858975703898,
    0.01134122086446862,
    0.011336582753233343,
    0.011325257495737607,
    0.01112782476233401,
    0.010921388055911025,
    0.010702214079751192,
    0.010450094641449727,
    0.010211844405074336,
    0.009979351954918376,
    0.009744680916984664,
    0.009507908133875373,
    0.009275153013848954,
    0.009052339016737024,
    0.00881871641535196,
    0.008576566685804539,
    0.008328304419301831,
    0.008074476173103688,
    0.007818127187450648,
    0.007560140226584159,
    0.007303348815084027,
    0.0070445778610827865,
    0.006782624840642196,
    0.006535056381036616,
    0.006294141263107,
    0.0060379680202159135,
    0.005755048084115355,
    0.005471607539161846,
    0.005213055482193659,
    0.004995116786781104,
    0.004748523651963094,
    0.004427810081744682,
    0.004092481269869861,
    0.003807261521899033,
    0.0036136856982531954,
    0.003471719759199435,
    0.0033329579843401247,
    0.0031977658456304276,
    0.0030666289604069277,
    0.0029391388352792295,
    0.002814547163875543,
    0.0026893197964230823,
    0.002561458121412205,
    0.002430649075927786,
    0.002288373770275667,
    0.0021346568815259954,
    0.001976059147950849,
    0.0018108850944824836,
    0.001640095346528458,
    0.0014719363782504137,
    0.0013062763169802917,
    0.0011443683270304145,
    0.000994175548512676,
    0.0008610419093733643,
    0.0007472059679234679,
    0.0006478122995357113,
    0.0005593727683755132,
    0.00048289564581939173,
    0.0004183710867629954,
    0.0003651598389134773,
    0.0003229360278848907,
    0.00028848789306159313,
    0.0002578827304637225,
    0.00022851066698475258,
    0.0001988319725071463,
    0.00016892785478387308,
    0.00014035954118032747,
    0.00011421944168849987,
    9.116483275882686e-05,
    7.147195358481961e-05,
    5.203388077770891e-05,
    3.493176465217351e-05,
    2.3310454586076205e-05,
    1.6168166397516464e-05,
    1.3113339008168005e-05,
    1.185455219482598e-05,
    1.0900687128028838e-05,
    1.0089625395325778e-05,
    9.30259125552525e-06,
    8.454647774804493e-06,
    7.52119003314466e-06,
    6.575891433621866e-06,
    5.686110520411272e-06,
    4.89309511300197e-06,
    4.3082665220437435e-06,
    3.7960540958291135e-06,
    3.333898952255864e-06,
    2.924594292629934e-06,
    2.5678811883184543e-06,
    2.261430638762962e-06,
    1.99959450257748e-06,
    1.7733068970988392e-06,
    1.5755591432772128e-06,
    1.4009411323203946e-06,
    1.1670855968943561e-06,
    9.858488336898767e-07,
    8.329736107853534e-07,
    6.998621078700888e-07,
    5.812614853077172e-07,
    4.743385664030164e-07,
    3.7759378749874884e-07,
    2.9287851766490877e-07,
    2.2198426962860873e-07,
    1.6483812410639142e-07,
    1.3313029214567614e-07,
    1.057153431673223e-07,
    8.211489927538205e-08,
    6.291961411288048e-08,
    4.796093325657464e-08,
    3.6676174703032866e-08,
    2.83326896665341e-08,
    2.1980924760802775e-08,
    1.694211289087208e-08,
    1.2833990485506551e-08,
    9.454556988210996e-09,
    6.471606180090994e-09,
    4.125158277636932e-09,
    2.6035607297712525e-09,
    1.7299529438569526e-09,
    1.286711764541813e-09,
    1.061398337419818e-09,
    9.097620486816128e-10,
    7.922623397996074e-10,
    6.853958858420519e-10,
    5.759499705971969e-10,
    4.6857203369096e-10,
    3.6933801792672846e-10,
    2.774357659623978e-10,
    1.907974126381206e-10,
    1.1872394475780025e-10,
    6.860405431069818e-11,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]

P_3MV_40x40_common = _coerce_to_grid(
    P_3MV_40x40,
    E_40x40_COMMON,
)

P_6MV_40x40_common = _coerce_to_grid(
    P_6MV_40x40,
    E_40x40_COMMON,
)

P_10MV_40x40_common = _coerce_to_grid(
    P_10MV_40x40,
    E_40x40_COMMON,
)

P_15MV_40x40_common = _coerce_to_grid(
    P_15MV_40x40,
    E_40x40_COMMON,
)

P_16MV_40x40_common = _coerce_to_grid(
    P_16MV_40x40,
    E_40x40_COMMON,
)

P_18MV_40x40_common = _coerce_to_grid(
    P_18MV_40x40,
    E_40x40_COMMON,
)

SPECTRUM_LIBRARY: Dict[str, Dict[str, Any]] = {
    '3MV_40x40_extrapolated': _make_spectrum(
        case='3MV_40x40_extrapolated',
        beam_mv=3,
        field_cm='40x40',
        peak_mev=0.0,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_3MV_40x40_common,
        origin=(
            'linear_extrapolation_from_6MV_and_18MV_same_field'
        ),
        interpolation_method=(
            'linear_extrapolation_from_6MV_and_18MV_same_field'
        ),
    ),
    '6MV_40x40': _make_spectrum(
        case='6MV_40x40',
        beam_mv=6,
        field_cm='40x40',
        peak_mev=0.3,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_6MV_40x40_common,
        origin='project_processed_from_Ding2002_Figure3_embedded',
    ),
    '10MV_40x40_interpolated': _make_spectrum(
        case='10MV_40x40_interpolated',
        beam_mv=10,
        field_cm='40x40',
        peak_mev=0.0,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_10MV_40x40_common,
        origin='Version_9_6MV_shape_scaled_to_10MV',
        interpolation_method='energy_axis_scaled_by_10_over_6',
    ),
    '15MV_40x40_interpolated': _make_spectrum(
        case='15MV_40x40_interpolated',
        beam_mv=15,
        field_cm='40x40',
        peak_mev=0.0,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_15MV_40x40_common,
        origin='linear_between_6MV_and_18MV_same_field',
        interpolation_method='linear_between_6MV_and_18MV_same_field',
    ),
    '16MV_40x40_interpolated': _make_spectrum(
        case='16MV_40x40_interpolated',
        beam_mv=16,
        field_cm='40x40',
        peak_mev=0.0,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_16MV_40x40_common,
        origin='linear_between_6MV_and_18MV_same_field',
        interpolation_method='linear_between_6MV_and_18MV_same_field',
    ),
    '18MV_40x40': _make_spectrum(
        case='18MV_40x40',
        beam_mv=18,
        field_cm='40x40',
        peak_mev=0.5,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_18MV_40x40_common,
        origin='project_processed_from_Ding2002_Figure3_embedded',
    ),
}

SPECTRUM_FILES_USED: List[str] = []


def get_available_spectrum_cases() -> List[str]:
    return sorted(SPECTRUM_LIBRARY.keys())


def get_spectrum_origin(case: str) -> str:
    return str(
        SPECTRUM_LIBRARY.get(
            case,
            {},
        ).get(
            'origin',
            'unknown',
        )
    )


def get_spectrum_provenance_note(case: str) -> str:
    origin = get_spectrum_origin(case)
    if origin.startswith("project_processed_from_Ding2002"):
        return "Project-processed/smoothed spectrum based on embedded Ding 2002 digitization; not asserted as a raw digitization or commissioned machine spectrum"
    return origin


def get_spectrum_revision(case: str) -> str:
    return "2026-08-11-v7.2-provenance-metadata-only"


def get_spectrum_array_sha256(case: str) -> str:
    entry = SPECTRUM_LIBRARY.get(case, {})
    energy = np.asarray(entry.get("energy_mev_center", []), dtype=np.float64)
    probability = np.asarray(entry.get("probability_mass_bin", []), dtype=np.float64)
    digest = hashlib.sha256()
    digest.update(energy.tobytes(order="C"))
    digest.update(probability.tobytes(order="C"))
    return digest.hexdigest()


def get_source_mode(cfg) -> str:
    return str(
        getattr(
            cfg,
            'source_mode',
            'mono',
        )
    ).strip().lower()


def get_spectrum_case(cfg) -> str:
    return str(
        getattr(
            cfg,
            'spectrum_case',
            '18MV_40x40',
        )
    ).strip()


def get_digitized_spectrum(cfg) -> Dict[str, Any]:
    case = get_spectrum_case(cfg)

    if case not in SPECTRUM_LIBRARY:
        raise ValueError(
            f'Unknown spectrum_case={case}. '
            f'Available: {get_available_spectrum_cases()}'
        )

    spec = deepcopy(SPECTRUM_LIBRARY[case])

    energies = np.asarray(
        spec['energy_mev_center'],
        dtype=float,
    )

    probs = np.asarray(
        spec['probability_mass_bin'],
        dtype=float,
    )

    if energies.ndim != 1 or probs.ndim != 1:
        raise ValueError(
            f'{case}: energy and probability arrays '
            f'must be one-dimensional'
        )

    if len(energies) != len(probs):
        raise ValueError(
            f'{case}: energy grid length {len(energies)} '
            f'!= probability length {len(probs)}'
        )

    if not np.all(np.isfinite(energies)):
        raise ValueError(
            f'{case}: energy grid contains NaN or infinite values'
        )

    if not np.all(np.isfinite(probs)):
        raise ValueError(
            f'{case}: probability array contains NaN or infinite values'
        )

    probs = np.maximum(probs, 0.0)
    positive_indices = np.flatnonzero(probs > 0.0)

    if positive_indices.size == 0:
        raise ValueError(
            f'Spectrum case {case} has non-positive total probability'
        )

    last_positive_index = int(positive_indices[-1])

    energies = energies[:last_positive_index + 1]
    probs = probs[:last_positive_index + 1]

    prob_sum = float(probs.sum())

    if prob_sum <= 0.0:
        raise ValueError(
            f'Spectrum case {case} has non-positive total probability'
        )

    probs = probs / prob_sum
    bin_width_mev = float(
        spec.get(
            'bin_width_mev',
            0.1,
        )
    )

    spec['energy_mev_center'] = energies
    spec['probability_mass_bin'] = probs
    spec['mean_energy_mev'] = float(
        np.sum(energies * probs)
    )
    spec['max_energy_mev'] = float(
        energies[-1] + 0.5 * bin_width_mev
    )
    spec['last_positive_bin_index'] = last_positive_index
    spec['case'] = case
    spec['origin'] = get_spectrum_origin(case)
    spec['files_used'] = list(SPECTRUM_FILES_USED)

    return spec


def source_energy_hist_max_mev(cfg) -> float:
    if get_source_mode(cfg) == 'digitized_spectrum':
        return get_digitized_spectrum(cfg)['max_energy_mev']

    return max(
        float(derived_source_energy_mev(cfg)),
        0.1,
    )


def sample_source_energy_mev(cfg) -> float:
    if get_source_mode(cfg) == 'digitized_spectrum':
        spec = get_digitized_spectrum(cfg)

        return float(
            np.random.choice(
                spec['energy_mev_center'],
                p=spec['probability_mass_bin'],
            )
        )

    return float(derived_source_energy_mev(cfg))


def derived_source_energy_mev(cfg) -> float:
    if get_source_mode(cfg) == 'digitized_spectrum':
        return get_digitized_spectrum(cfg)['mean_energy_mev']

    return float(cfg.mev) / float(
        getattr(
            cfg,
            'mono_energy_divisor',
            1.6,
        )
    )


def clone_cfg(cfg, **updates):
    d = deepcopy(vars(cfg))
    d.update(updates)

    return SimpleNamespace(**d)


def vec_cm(v) -> Tuple[float, float, float]:
    return (
        float(v.x / cm),
        float(v.y / cm),
        float(v.z / cm),
    )


def vec_x_cm(v) -> float:
    return float(v.x / cm)


def vec_x_unit(v) -> float:
    return float(v.x)


def safe_sem_from_sum_sumsq(
    sum_arr: np.ndarray,
    sumsq_arr: np.ndarray,
    n_events: int,
) -> np.ndarray:
    if n_events <= 1:
        return np.zeros_like(
            sum_arr,
            dtype=float,
        )

    mean = sum_arr / n_events
    ex2 = sumsq_arr / n_events
    var = np.maximum(
        ex2 - mean * mean,
        0.0,
    )
    sem = np.sqrt(var / n_events)

    return sem


def binomial_sigma(
    count: int,
    total: int,
) -> float:
    if total <= 0:
        return 0.0

    p = count / total

    return math.sqrt(
        max(
            p * (1.0 - p) / total,
            0.0,
        )
    )


def propagate_ratio_sigma(
    a: float,
    sa: float,
    b: float,
    sb: float,
) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0

    r = a / b

    return abs(r) * math.sqrt(
        (sa / a) ** 2 + (sb / b) ** 2
    )


class EventState:
    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self.primary_track_id: Optional[int] = None
        self.collided_any: bool = False
        self.compton_count: int = 0

        self.had_interaction: bool = False
        self.first_interaction_x_cm: Optional[float] = None
        self.last_interaction_x_cm: Optional[float] = None

        self.terminated: Optional[str] = None
        self.exit_energy_mev: Optional[float] = None
        self.exit_cosine: Optional[float] = None

        self.track_points_cm: List[
            Tuple[float, float, float]
        ] = []

        self.record_this_track: bool = False
        self.primary_processes_seen: set = set()

        self.fluence_event_cm: Optional[np.ndarray] = None
        self.edep_event_mev: Optional[np.ndarray] = None
        self.primary_deltaE_event_mev: Optional[np.ndarray] = None


class SimulationTallies:
    def __init__(self, cfg) -> None:
        self.cfg = cfg

        self.source_energy_mev = derived_source_energy_mev(cfg)
        self.source_energy_hist_max_mev = source_energy_hist_max_mev(cfg)

        self.depth_edges_cm = np.linspace(
            0.0,
            cfg.wall_thickness_cm,
            cfg.depth_nbins + 1,
            dtype=float,
        )

        self.depth_mid_cm = 0.5 * (
            self.depth_edges_cm[:-1]
            + self.depth_edges_cm[1:]
        )

        self.depth_bin_widths_cm = np.diff(
            self.depth_edges_cm
        )

        self.fluence_tracklen_cm = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.fluence_tracklen_cm_sumsq = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.edep_vs_depth_mev = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.edep_vs_depth_mev_sumsq = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.primary_gamma_deltaE_vs_depth_mev = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.primary_gamma_deltaE_vs_depth_mev_sumsq = np.zeros(
            cfg.depth_nbins,
            dtype=np.float64,
        )

        self.first_interaction_hist = np.zeros(
            cfg.depth_nbins,
            dtype=np.int64,
        )

        self.last_interaction_hist = np.zeros(
            cfg.depth_nbins,
            dtype=np.int64,
        )

        self.scatter_hist = np.zeros(
            cfg.max_scatter_bin + 2,
            dtype=np.int64,
        )

        self.energy_hist_edges_mev = np.linspace(
            0.0,
            self.source_energy_hist_max_mev,
            cfg.energy_hist_bins + 1,
            dtype=float,
        )

        self.energy_hist_counts = np.zeros(
            cfg.energy_hist_bins,
            dtype=np.int64,
        )

        self.exit_cosine_hist_edges = np.linspace(
            -1.0,
            1.0,
            cfg.angle_hist_bins + 1,
            dtype=float,
        )

        self.exit_cosine_hist_counts = np.zeros(
            cfg.angle_hist_bins,
            dtype=np.int64,
        )

        self.corr_cos_edges = np.linspace(
            -1.0,
            1.0,
            cfg.corr_angle_bins + 1,
            dtype=float,
        )

        self.corr_energy_edges_mev = np.linspace(
            0.0,
            self.source_energy_hist_max_mev,
            cfg.corr_energy_bins + 1,
            dtype=float,
        )

        self.exit_energy_angle_hist = np.zeros(
            (
                cfg.corr_angle_bins,
                cfg.corr_energy_bins,
            ),
            dtype=np.int64,
        )

        self.transmitted_energy_sample_mev: List[float] = []
        self.exit_cosine_sample: List[float] = []
        self.exit_sample_seen: int = 0

        self.tracks_cm: List[
            Tuple[
                List[float],
                List[float],
                List[float],
            ]
        ] = []

        self.tracks_saved: int = 0
        self.tracks_rejected_short: int = 0

        self.transmitted_count: int = 0
        self.uncollided_transmitted_count: int = 0
        self.backscatter_count: int = 0
        self.removed_count: int = 0
        self.other_escape_count: int = 0

        self.process_counts_all_tracks = {
            p: 0
            for p in PRIMARY_GAMMA_PROCESSES
        }

        self.process_counts_primary_total = {
            p: 0
            for p in PRIMARY_GAMMA_PROCESSES
        }

        self.process_counts_primary_eventwise = {
            p: 0
            for p in PRIMARY_GAMMA_PROCESSES
        }

        self.buildup_primary_counts = {
            'uncollided': 0,
            'rayleigh_only': 0,
            'compton_only': 0,
            'mixed_rayleigh_compton': 0,
        }

        self.secondary_energy_hist_edges_mev = np.linspace(
            0.0,
            self.source_energy_hist_max_mev,
            cfg.secondary_energy_hist_bins + 1,
            dtype=float,
        )

        self.secondary_exit_hist_counts = {
            pname: np.zeros(
                cfg.secondary_energy_hist_bins,
                dtype=np.int64,
            )
            for pname in SECONDARY_PARTICLES_TO_TALLY
        }

        self.secondary_exit_counts_forward = {
            pname: 0
            for pname in SECONDARY_PARTICLES_TO_TALLY
        }

        self.material_name_used: Optional[str] = None
        self.material_density_g_cm3: Optional[float] = None

    def depth_bin_index(
        self,
        x_cm: float,
    ) -> int:
        i = int(
            np.searchsorted(
                self.depth_edges_cm,
                x_cm,
                side='right',
            ) - 1
        )

        if i < 0:
            return 0

        if i >= len(self.depth_mid_cm):
            return len(self.depth_mid_cm) - 1

        return i

    def tally_track_length_1d_into(
        self,
        bins_arr: np.ndarray,
        x0_cm: float,
        x1_cm: float,
        step_length_cm: float,
    ) -> None:
        if step_length_cm <= 0.0:
            return

        if x1_cm == x0_cm:
            if (
                x0_cm < self.depth_edges_cm[0]
                or x0_cm > self.depth_edges_cm[-1]
            ):
                return

            i = (
                np.searchsorted(
                    self.depth_edges_cm,
                    x0_cm,
                    side='right',
                )
                - 1
            )

            if 0 <= i < len(bins_arr):
                bins_arr[i] += step_length_cm

            return

        xmin = min(x0_cm, x1_cm)
        xmax = max(x0_cm, x1_cm)

        if (
            xmax <= self.depth_edges_cm[0]
            or xmin >= self.depth_edges_cm[-1]
        ):
            return

        i0 = max(
            0,
            np.searchsorted(
                self.depth_edges_cm,
                xmin,
                side='right',
            ) - 1,
        )

        i1 = min(
            len(bins_arr) - 1,
            np.searchsorted(
                self.depth_edges_cm,
                xmax,
                side='left',
            ) - 1,
        )

        dx = x1_cm - x0_cm
        inv_dx = 1.0 / dx

        for i in range(i0, i1 + 1):
            a = self.depth_edges_cm[i]
            b = self.depth_edges_cm[i + 1]

            t_a = (a - x0_cm) * inv_dx
            t_b = (b - x0_cm) * inv_dx

            t_lo = max(
                0.0,
                min(t_a, t_b),
            )

            t_hi = min(
                1.0,
                max(t_a, t_b),
            )

            if t_hi > t_lo:
                bins_arr[i] += step_length_cm * (
                    t_hi - t_lo
                )

    def _increment_1d_hist(
        self,
        value: float,
        edges: np.ndarray,
        counts: np.ndarray,
    ) -> None:
        if value < edges[0] or value > edges[-1]:
            return

        idx = int(
            np.searchsorted(
                edges,
                value,
                side='right',
            ) - 1
        )

        if idx == len(counts):
            idx -= 1

        if 0 <= idx < len(counts):
            counts[idx] += 1

    def _increment_2d_hist(
        self,
        x: float,
        y: float,
    ) -> None:
        if (
            x < self.corr_cos_edges[0]
            or x > self.corr_cos_edges[-1]
        ):
            return

        if (
            y < self.corr_energy_edges_mev[0]
            or y > self.corr_energy_edges_mev[-1]
        ):
            return

        ix = int(
            np.searchsorted(
                self.corr_cos_edges,
                x,
                side='right',
            ) - 1
        )

        iy = int(
            np.searchsorted(
                self.corr_energy_edges_mev,
                y,
                side='right',
            ) - 1
        )

        if ix == self.exit_energy_angle_hist.shape[0]:
            ix -= 1

        if iy == self.exit_energy_angle_hist.shape[1]:
            iy -= 1

        if (
            0 <= ix < self.exit_energy_angle_hist.shape[0]
            and 0 <= iy < self.exit_energy_angle_hist.shape[1]
        ):
            self.exit_energy_angle_hist[ix, iy] += 1

    def append_exit_sample(
        self,
        energy_mev: float,
        cosine: float,
    ) -> None:
        self.exit_sample_seen += 1

        if (
            len(self.transmitted_energy_sample_mev)
            < self.cfg.max_raw_exit_samples
        ):
            self.transmitted_energy_sample_mev.append(
                float(energy_mev)
            )

            self.exit_cosine_sample.append(
                float(cosine)
            )

            return

        j = np.random.randint(
            0,
            self.exit_sample_seen,
        )

        if j < self.cfg.max_raw_exit_samples:
            self.transmitted_energy_sample_mev[j] = float(
                energy_mev
            )

            self.exit_cosine_sample[j] = float(
                cosine
            )

    def try_store_track(
        self,
        pts_cm: Sequence[
            Tuple[float, float, float]
        ],
    ) -> None:
        if not pts_cm:
            return

        decimated = pts_cm[
            ::max(
                1,
                self.cfg.obj_decimate,
            )
        ]

        if len(decimated) < self.cfg.min_pts_after_decimate:
            self.tracks_rejected_short += 1
            return

        xs = [p[0] for p in pts_cm]
        ys = [p[1] for p in pts_cm]
        zs = [p[2] for p in pts_cm]

        self.tracks_cm.append(
            (
                xs,
                ys,
                zs,
            )
        )

        self.tracks_saved += 1


class MyDetectorConstruction(
    G4VUserDetectorConstruction
):
    def __init__(
        self,
        cfg,
        tallies: SimulationTallies,
    ):
        super().__init__()

        self.cfg = cfg
        self.tallies = tallies

        self.world_half_x = max(
            5.0 * m,
            0.5
            * (
                cfg.wall_thickness_cm
                + 2.0 * cfg.source_to_wall_gap_cm
                + 500.0
            )
            * cm,
        )

        self.world_half_y = max(
            5.0 * m,
            0.6 * cfg.wall_size_y_cm * cm,
        )

        self.world_half_z = max(
            5.0 * m,
            0.6 * cfg.wall_size_z_cm * cm,
        )

    def _build_concrete_material(self):
        nist = G4NistManager.Instance()

        if self.cfg.concrete_density_g_cm3 is None:
            mat = nist.FindOrBuildMaterial(
                'G4_CONCRETE'
            )

            self.tallies.material_name_used = 'G4_CONCRETE'

            self.tallies.material_density_g_cm3 = float(
                mat.GetDensity() / (g / cm3)
            )

            return mat

        if hasattr(
            nist,
            'BuildMaterialWithNewDensity',
        ):
            mat = nist.BuildMaterialWithNewDensity(
                'CustomConcrete',
                'G4_CONCRETE',
                self.cfg.concrete_density_g_cm3
                * g
                / cm3,
            )

            self.tallies.material_name_used = (
                'CustomConcrete_from_G4_CONCRETE'
            )

            self.tallies.material_density_g_cm3 = float(
                mat.GetDensity() / (g / cm3)
            )

            return mat

        raise RuntimeError(
            'geant4_pybind does not expose '
            'BuildMaterialWithNewDensity() here. '
            'Use concrete_density_g_cm3=None, or manually '
            'define a custom G4Material.'
        )

    def Construct(self):
        nist = G4NistManager.Instance()

        air = nist.FindOrBuildMaterial(
            'G4_AIR'
        )

        concrete = self._build_concrete_material()

        solidWorld = G4Box(
            'World',
            self.world_half_x,
            self.world_half_y,
            self.world_half_z,
        )

        logicWorld = G4LogicalVolume(
            solidWorld,
            air,
            'World',
        )

        physWorld = G4PVPlacement(
            None,
            G4ThreeVector(),
            logicWorld,
            'World',
            None,
            False,
            0,
        )

        shield_half_x = (
            0.5
            * self.cfg.wall_thickness_cm
            * cm
        )

        shield_half_y = (
            0.5
            * self.cfg.wall_size_y_cm
            * cm
        )

        shield_half_z = (
            0.5
            * self.cfg.wall_size_z_cm
            * cm
        )

        shield_center_x = (
            0.5
            * self.cfg.wall_thickness_cm
            * cm
        )

        solidShield = G4Box(
            'Shield',
            shield_half_x,
            shield_half_y,
            shield_half_z,
        )

        logicShield = G4LogicalVolume(
            solidShield,
            concrete,
            'Shield',
        )

        G4PVPlacement(
            None,
            G4ThreeVector(
                shield_center_x,
                0.0,
                0.0,
            ),
            logicShield,
            'Shield',
            logicWorld,
            False,
            0,
        )

        return physWorld


class MyPrimaryGeneratorAction(
    G4VUserPrimaryGeneratorAction
):
    def __init__(self, cfg):
        super().__init__()

        self.cfg = cfg
        self.particleGun = G4ParticleGun(1)

        self.particleGun.SetParticleDefinition(
            G4Gamma.Gamma()
        )

    def GeneratePrimaries(self, anEvent):
        x_cm = -float(
            self.cfg.source_to_wall_gap_cm
        )

        y_cm = 0.0
        z_cm = 0.0

        self.particleGun.SetParticlePosition(
            G4ThreeVector(
                x_cm * cm,
                y_cm * cm,
                z_cm * cm,
            )
        )

        source_energy_mev = sample_source_energy_mev(
            self.cfg
        )

        self.particleGun.SetParticleEnergy(
            source_energy_mev * MeV
        )

        alpha = np.deg2rad(
            float(self.cfg.beam_half_angle_deg)
        )

        u = np.random.random()
        v = np.random.random()

        cos_theta = (
            1.0
            - u * (1.0 - np.cos(alpha))
        )

        sin_theta = np.sqrt(
            max(
                0.0,
                1.0 - cos_theta * cos_theta,
            )
        )

        phi = 2.0 * np.pi * v

        ux = cos_theta
        uy = sin_theta * np.cos(phi)
        uz = sin_theta * np.sin(phi)

        self.particleGun.SetParticleMomentumDirection(
            G4ThreeVector(
                float(ux),
                float(uy),
                float(uz),
            )
        )

        self.particleGun.GeneratePrimaryVertex(
            anEvent
        )


class MyPhysicsList(
    G4VModularPhysicsList
):
    def __init__(self, cfg):
        super().__init__()

        self.cfg = cfg

        if 'G4EmStandardPhysics_option4' in globals():
            self.RegisterPhysics(
                G4EmStandardPhysics_option4()
            )
        elif 'G4EmStandardPhysicsOption4' in globals():
            self.RegisterPhysics(
                G4EmStandardPhysicsOption4()
            )
        else:
            raise RuntimeError(
                'Could not find EM Option 4 constructor '
                'in geant4_pybind. Expected '
                'G4EmStandardPhysics_option4 or '
                'G4EmStandardPhysicsOption4.'
            )

    def SetCuts(self):
        cut = (
            self.cfg.production_cut_mm
            * mm
        )

        self.SetCutValue(
            cut,
            'gamma',
        )

        self.SetCutValue(
            cut,
            'e-',
        )

        self.SetCutValue(
            cut,
            'e+',
        )

        try:
            self.SetCutValue(
                cut,
                'proton',
            )
        except Exception:
            pass


class MyEventAction(
    G4UserEventAction
):
    def __init__(
        self,
        cfg,
        tallies: SimulationTallies,
        state: EventState,
    ):
        super().__init__()

        self.cfg = cfg
        self.tallies = tallies
        self.state = state

    def BeginOfEventAction(
        self,
        event,
    ):
        self.state.reset()

        self.state.record_this_track = (
            self.tallies.tracks_saved
            < self.cfg.n_visual_tracks
        )

        self.state.fluence_event_cm = np.zeros(
            self.cfg.depth_nbins,
            dtype=float,
        )

        self.state.edep_event_mev = np.zeros(
            self.cfg.depth_nbins,
            dtype=float,
        )

        self.state.primary_deltaE_event_mev = np.zeros(
            self.cfg.depth_nbins,
            dtype=float,
        )

    def EndOfEventAction(
        self,
        event,
    ):
        self.tallies.fluence_tracklen_cm += (
            self.state.fluence_event_cm
        )

        self.tallies.fluence_tracklen_cm_sumsq += (
            self.state.fluence_event_cm ** 2
        )

        self.tallies.edep_vs_depth_mev += (
            self.state.edep_event_mev
        )

        self.tallies.edep_vs_depth_mev_sumsq += (
            self.state.edep_event_mev ** 2
        )

        self.tallies.primary_gamma_deltaE_vs_depth_mev += (
            self.state.primary_deltaE_event_mev
        )

        self.tallies.primary_gamma_deltaE_vs_depth_mev_sumsq += (
            self.state.primary_deltaE_event_mev ** 2
        )

        for pname in self.state.primary_processes_seen:
            if (
                pname
                in self.tallies.process_counts_primary_eventwise
            ):
                self.tallies.process_counts_primary_eventwise[
                    pname
                ] += 1

        if self.state.terminated == 'transmitted':
            self.tallies.transmitted_count += 1

            if not self.state.collided_any:
                self.tallies.uncollided_transmitted_count += 1

                self.tallies.buildup_primary_counts[
                    'uncollided'
                ] += 1

            rayl_seen = (
                'Rayl'
                in self.state.primary_processes_seen
            )

            compt_seen = (
                'compt'
                in self.state.primary_processes_seen
            )

            if self.state.collided_any:
                if rayl_seen and compt_seen:
                    self.tallies.buildup_primary_counts[
                        'mixed_rayleigh_compton'
                    ] += 1
                elif rayl_seen:
                    self.tallies.buildup_primary_counts[
                        'rayleigh_only'
                    ] += 1
                elif compt_seen:
                    self.tallies.buildup_primary_counts[
                        'compton_only'
                    ] += 1

            idx = self.state.compton_count

            if idx > self.cfg.max_scatter_bin:
                idx = self.cfg.max_scatter_bin + 1

            self.tallies.scatter_hist[idx] += 1

            if (
                self.state.exit_energy_mev is not None
                and self.state.exit_cosine is not None
            ):
                self.tallies._increment_1d_hist(
                    self.state.exit_energy_mev,
                    self.tallies.energy_hist_edges_mev,
                    self.tallies.energy_hist_counts,
                )

                self.tallies._increment_1d_hist(
                    self.state.exit_cosine,
                    self.tallies.exit_cosine_hist_edges,
                    self.tallies.exit_cosine_hist_counts,
                )

                self.tallies._increment_2d_hist(
                    self.state.exit_cosine,
                    self.state.exit_energy_mev,
                )

                self.tallies.append_exit_sample(
                    self.state.exit_energy_mev,
                    self.state.exit_cosine,
                )

        elif self.state.terminated == 'backescaped':
            self.tallies.backscatter_count += 1

        elif self.state.terminated == 'removed':
            self.tallies.removed_count += 1

        else:
            self.tallies.other_escape_count += 1

        if self.state.had_interaction:
            self.tallies.first_interaction_hist[
                self.tallies.depth_bin_index(
                    self.state.first_interaction_x_cm
                )
            ] += 1

            self.tallies.last_interaction_hist[
                self.tallies.depth_bin_index(
                    self.state.last_interaction_x_cm
                )
            ] += 1

        if (
            self.state.record_this_track
            and self.state.track_points_cm
        ):
            self.tallies.try_store_track(
                self.state.track_points_cm
            )


class MyTrackingAction(
    G4UserTrackingAction
):
    def __init__(
        self,
        cfg,
        state: EventState,
    ):
        super().__init__()

        self.cfg = cfg
        self.state = state

    def PreUserTrackingAction(
        self,
        track,
    ):
        if (
            track.GetParentID() == 0
            and track.GetDefinition().GetParticleName()
            == 'gamma'
        ):
            self.state.primary_track_id = (
                track.GetTrackID()
            )

    def PostUserTrackingAction(
        self,
        track,
    ):
        if self.state.primary_track_id is None:
            return

        if (
            track.GetTrackID()
            != self.state.primary_track_id
        ):
            return

        if (
            track.GetDefinition().GetParticleName()
            != 'gamma'
        ):
            return

        if self.state.terminated is not None:
            return

        x_cm = vec_x_cm(
            track.GetPosition()
        )

        if x_cm > self.cfg.wall_thickness_cm:
            self.state.terminated = 'transmitted'

            self.state.exit_energy_mev = float(
                track.GetKineticEnergy() / MeV
            )

            self.state.exit_cosine = vec_x_unit(
                track.GetMomentumDirection()
            )

        elif x_cm < 0.0:
            self.state.terminated = 'backescaped'

        elif 0.0 <= x_cm <= self.cfg.wall_thickness_cm:
            self.state.terminated = 'removed'

        else:
            self.state.terminated = 'other_escape'


class MySteppingAction(
    G4UserSteppingAction
):
    INTERACTION_IGNORE = {
        'Transportation',
        'CoupledTransportation',
        'StepLimiter',
        'UserSpecialCut',
    }

    def __init__(
        self,
        cfg,
        tallies: SimulationTallies,
        state: EventState,
    ):
        super().__init__()

        self.cfg = cfg
        self.tallies = tallies
        self.state = state

    def UserSteppingAction(
        self,
        step,
    ):
        track = step.GetTrack()
        pre = step.GetPreStepPoint()
        post = step.GetPostStepPoint()

        pre_pos = pre.GetPosition()
        post_pos = post.GetPosition()

        x0_cm, y0_cm, z0_cm = vec_cm(
            pre_pos
        )

        x1_cm, y1_cm, z1_cm = vec_cm(
            post_pos
        )

        proc = post.GetProcessDefinedStep()

        proc_name = (
            proc.GetProcessName()
            if proc is not None
            else None
        )

        if (
            proc_name
            in self.tallies.process_counts_all_tracks
        ):
            self.tallies.process_counts_all_tracks[
                proc_name
            ] += 1

        edep_mev = float(
            step.GetTotalEnergyDeposit() / MeV
        )

        if edep_mev > 0.0:
            x_mid_cm = max(
                0.0,
                min(
                    self.cfg.wall_thickness_cm,
                    0.5 * (x0_cm + x1_cm),
                ),
            )

            self.state.edep_event_mev[
                self.tallies.depth_bin_index(
                    x_mid_cm
                )
            ] += edep_mev

        pname = (
            track.GetDefinition().GetParticleName()
        )

        if (
            x0_cm <= self.cfg.wall_thickness_cm
            and x1_cm > self.cfg.wall_thickness_cm
        ):
            if (
                track.GetParentID() > 0
                and pname
                in self.tallies.secondary_exit_hist_counts
            ):
                ekin_mev = float(
                    post.GetKineticEnergy() / MeV
                )

                self.tallies._increment_1d_hist(
                    ekin_mev,
                    self.tallies.secondary_energy_hist_edges_mev,
                    self.tallies.secondary_exit_hist_counts[
                        pname
                    ],
                )

                self.tallies.secondary_exit_counts_forward[
                    pname
                ] += 1

        if self.state.primary_track_id is None:
            return

        if (
            track.GetTrackID()
            != self.state.primary_track_id
        ):
            return

        if pname != 'gamma':
            return

        if self.state.record_this_track:
            if not self.state.track_points_cm:
                self.state.track_points_cm.append(
                    (
                        x0_cm,
                        y0_cm,
                        z0_cm,
                    )
                )

            last = self.state.track_points_cm[-1]

            new_pt = (
                x1_cm,
                y1_cm,
                z1_cm,
            )

            if new_pt != last:
                self.state.track_points_cm.append(
                    new_pt
                )

        self.tallies.tally_track_length_1d_into(
            self.state.fluence_event_cm,
            x0_cm=x0_cm,
            x1_cm=x1_cm,
            step_length_cm=float(
                step.GetStepLength() / cm
            ),
        )

        if (
            proc_name is not None
            and proc_name
            not in self.INTERACTION_IGNORE
        ):
            self.state.collided_any = True

            if (
                proc_name
                in self.tallies.process_counts_primary_total
            ):
                self.tallies.process_counts_primary_total[
                    proc_name
                ] += 1

                self.state.primary_processes_seen.add(
                    proc_name
                )

            x_int_cm = max(
                0.0,
                min(
                    self.cfg.wall_thickness_cm,
                    x1_cm,
                ),
            )

            if not self.state.had_interaction:
                self.state.had_interaction = True
                self.state.first_interaction_x_cm = x_int_cm

            self.state.last_interaction_x_cm = x_int_cm

            if proc_name == 'compt':
                self.state.compton_count += 1

                preE_mev = float(
                    pre.GetKineticEnergy() / MeV
                )

                postE_mev = float(
                    post.GetKineticEnergy() / MeV
                )

                dE_mev = max(
                    0.0,
                    preE_mev - postE_mev,
                )

                self.state.primary_deltaE_event_mev[
                    self.tallies.depth_bin_index(
                        x_int_cm
                    )
                ] += dE_mev

        if (
            self.state.terminated is None
            and x0_cm <= self.cfg.wall_thickness_cm
            and x1_cm > self.cfg.wall_thickness_cm
        ):
            self.state.terminated = 'transmitted'

            self.state.exit_energy_mev = float(
                post.GetKineticEnergy() / MeV
            )

            self.state.exit_cosine = vec_x_unit(
                track.GetMomentumDirection()
            )

        elif (
            self.state.terminated is None
            and x0_cm >= 0.0
            and x1_cm < 0.0
        ):
            self.state.terminated = 'backescaped'

        elif self.state.terminated is None:
            half_y = (
                0.5 * self.cfg.wall_size_y_cm
            )

            half_z = (
                0.5 * self.cfg.wall_size_z_cm
            )

            if (
                abs(y1_cm) > half_y
                or abs(z1_cm) > half_z
            ):
                self.state.terminated = 'other_escape'


class MyActionInitialization(
    G4VUserActionInitialization
):
    def __init__(
        self,
        cfg,
        tallies: SimulationTallies,
        state: EventState,
    ):
        super().__init__()

        self.cfg = cfg
        self.tallies = tallies
        self.state = state

        self.primary = None
        self.event = None
        self.tracking = None
        self.stepping = None

    def Build(self):
        self.primary = MyPrimaryGeneratorAction(
            self.cfg
        )

        self.event = MyEventAction(
            self.cfg,
            self.tallies,
            self.state,
        )

        self.tracking = MyTrackingAction(
            self.cfg,
            self.state,
        )

        self.stepping = MySteppingAction(
            self.cfg,
            self.tallies,
            self.state,
        )

        self.SetUserAction(
            self.primary
        )

        self.SetUserAction(
            self.event
        )

        self.SetUserAction(
            self.tracking
        )

        self.SetUserAction(
            self.stepping
        )


def export_scene_to_obj(
    tracks_cm: Sequence[
        Tuple[
            Sequence[float],
            Sequence[float],
            Sequence[float],
        ]
    ],
    filepath: str,
    scale: float,
    decimate: int,
    wall_thickness_cm: float,
    wall_size_y_cm: float,
    wall_size_z_cm: float,
    include_wall: bool = True,
    track_radius_cm: float = 0.5,
) -> Dict[str, Any]:
    decimate = max(
        1,
        int(decimate),
    )

    track_radius_cm = float(
        track_radius_cm
    )

    out_dir = (
        os.path.dirname(filepath)
        or '.'
    )

    base = os.path.splitext(
        os.path.basename(filepath)
    )[0]

    mtl_filename = base + '.mtl'

    mtl_filepath = os.path.join(
        out_dir,
        mtl_filename,
    )

    vertex_lines: List[str] = []
    object_lines: List[str] = []

    def add_vertex(
        p_cm: np.ndarray,
    ) -> int:
        x, y, z = p_cm * scale

        vertex_lines.append(
            f'v {x:.8g} {y:.8g} {z:.8g}'
        )

        return len(vertex_lines)

    def add_face(
        i0: int,
        i1: int,
        i2: int,
        i3: int,
    ) -> str:
        return f'f {i0} {i1} {i2} {i3}'

    wall_faces_written = 0

    if include_wall:
        object_lines.append('o Wall')
        object_lines.append('usemtl WallMaterial')

        x0 = 0.0
        x1 = wall_thickness_cm

        y0 = -0.5 * wall_size_y_cm
        y1 = 0.5 * wall_size_y_cm

        z0 = -0.5 * wall_size_z_cm
        z1 = 0.5 * wall_size_z_cm

        wall_vertices = [
            np.array(
                [x0, y0, z0],
                dtype=float,
            ),
            np.array(
                [x1, y0, z0],
                dtype=float,
            ),
            np.array(
                [x1, y1, z0],
                dtype=float,
            ),
            np.array(
                [x0, y1, z0],
                dtype=float,
            ),
            np.array(
                [x0, y0, z1],
                dtype=float,
            ),
            np.array(
                [x1, y0, z1],
                dtype=float,
            ),
            np.array(
                [x1, y1, z1],
                dtype=float,
            ),
            np.array(
                [x0, y1, z1],
                dtype=float,
            ),
        ]

        idx = [
            add_vertex(v)
            for v in wall_vertices
        ]

        wall_face_defs = [
            (
                idx[0],
                idx[1],
                idx[2],
                idx[3],
            ),
            (
                idx[4],
                idx[5],
                idx[6],
                idx[7],
            ),
            (
                idx[0],
                idx[1],
                idx[5],
                idx[4],
            ),
            (
                idx[3],
                idx[2],
                idx[6],
                idx[7],
            ),
            (
                idx[0],
                idx[3],
                idx[7],
                idx[4],
            ),
            (
                idx[1],
                idx[2],
                idx[6],
                idx[5],
            ),
        ]

        for f in wall_face_defs:
            object_lines.append(
                add_face(*f)
            )

            wall_faces_written += 1

    object_lines.append('o Tracks')
    object_lines.append('usemtl TrackMaterial')

    tracks_written = 0
    segment_prisms_written = 0

    def add_segment_prism(
        p0_cm: np.ndarray,
        p1_cm: np.ndarray,
        r_cm: float,
    ) -> bool:
        d = p1_cm - p0_cm
        length = np.linalg.norm(d)

        if length <= 1e-12:
            return False

        d = d / length

        ref = np.array(
            [0.0, 0.0, 1.0],
            dtype=float,
        )

        if abs(np.dot(d, ref)) > 0.9:
            ref = np.array(
                [0.0, 1.0, 0.0],
                dtype=float,
            )

        u = np.cross(d, ref)
        u_norm = np.linalg.norm(u)

        if u_norm <= 1e-12:
            ref = np.array(
                [1.0, 0.0, 0.0],
                dtype=float,
            )

            u = np.cross(d, ref)
            u_norm = np.linalg.norm(u)

            if u_norm <= 1e-12:
                return False

        u = u / u_norm
        v = np.cross(d, u)

        a = r_cm / np.sqrt(2.0)

        corners0 = [
            p0_cm + a * (u + v),
            p0_cm + a * (u - v),
            p0_cm + a * (-u - v),
            p0_cm + a * (-u + v),
        ]

        corners1 = [
            p1_cm + a * (u + v),
            p1_cm + a * (u - v),
            p1_cm + a * (-u - v),
            p1_cm + a * (-u + v),
        ]

        ids = [
            add_vertex(p)
            for p in corners0 + corners1
        ]

        face_defs = [
            (
                ids[0],
                ids[1],
                ids[2],
                ids[3],
            ),
            (
                ids[4],
                ids[5],
                ids[6],
                ids[7],
            ),
            (
                ids[0],
                ids[1],
                ids[5],
                ids[4],
            ),
            (
                ids[1],
                ids[2],
                ids[6],
                ids[5],
            ),
            (
                ids[2],
                ids[3],
                ids[7],
                ids[6],
            ),
            (
                ids[3],
                ids[0],
                ids[4],
                ids[7],
            ),
        ]

        for f in face_defs:
            object_lines.append(
                add_face(*f)
            )

        return True

    for xs, ys, zs in tracks_cm:
        xs_d = np.asarray(
            xs[::decimate],
            dtype=float,
        )

        ys_d = np.asarray(
            ys[::decimate],
            dtype=float,
        )

        zs_d = np.asarray(
            zs[::decimate],
            dtype=float,
        )

        if len(xs_d) < 2:
            continue

        wrote_any = False

        for i in range(len(xs_d) - 1):
            p0 = np.array(
                [
                    xs_d[i],
                    ys_d[i],
                    zs_d[i],
                ],
                dtype=float,
            )

            p1 = np.array(
                [
                    xs_d[i + 1],
                    ys_d[i + 1],
                    zs_d[i + 1],
                ],
                dtype=float,
            )

            ok = add_segment_prism(
                p0,
                p1,
                track_radius_cm,
            )

            if ok:
                segment_prisms_written += 1
                wrote_any = True

        if wrote_any:
            tracks_written += 1

    with open(
        mtl_filepath,
        'w',
        encoding='utf-8',
    ) as f:
        f.write('newmtl WallMaterial\n')
        f.write('Ka 0.2 0.2 0.2\n')
        f.write('Kd 0.7 0.7 0.7\n')
        f.write('Ks 0.1 0.1 0.1\n')
        f.write('Ns 10.0\n')
        f.write('illum 2\n\n')

        f.write('newmtl TrackMaterial\n')
        f.write('Ka 0.2 0.05 0.05\n')
        f.write('Kd 0.9 0.2 0.2\n')
        f.write('Ks 0.1 0.1 0.1\n')
        f.write('Ns 20.0\n')
        f.write('illum 2\n')

    with open(
        filepath,
        'w',
        encoding='utf-8',
    ) as f:
        f.write(
            f'mtllib {mtl_filename}\n'
        )

        for line in vertex_lines:
            f.write(line + '\n')

        for line in object_lines:
            f.write(line + '\n')

    return {
        'obj_filepath': filepath,
        'mtl_filepath': mtl_filepath,
        'tracks_written': tracks_written,
        'segment_prisms_written': segment_prisms_written,
        'vertices': len(vertex_lines),
        'faces': (
            wall_faces_written
            + segment_prisms_written * 6
        ),
        'wall_included': bool(include_wall),
        'track_radius_cm': track_radius_cm,
    }


def resultsafe_int(
    arr: np.ndarray,
) -> np.ndarray:
    return np.maximum(
        arr.astype(float),
        0.0,
    )


def build_results(
    cfg,
    tallies: SimulationTallies,
    runtime_s: float,
    obj_info: Optional[
        Dict[str, Any]
    ],
) -> Dict[str, Any]:
    N = max(
        1,
        int(cfg.num_photons),
    )

    source_energy_mev = derived_source_energy_mev(
        cfg
    )

    T_total = (
        tallies.transmitted_count
        / N
    )

    T_unc = (
        tallies.uncollided_transmitted_count
        / N
    )

    backscatter_fraction = (
        tallies.backscatter_count
        / N
    )

    removed_fraction = (
        tallies.removed_count
        / N
    )

    other_escape_fraction = (
        tallies.other_escape_count
        / N
    )

    sigma_T_total = binomial_sigma(
        tallies.transmitted_count,
        N,
    )

    sigma_T_unc = binomial_sigma(
        tallies.uncollided_transmitted_count,
        N,
    )

    sigma_backscatter = binomial_sigma(
        tallies.backscatter_count,
        N,
    )

    sigma_removed = binomial_sigma(
        tallies.removed_count,
        N,
    )

    sigma_other_escape = binomial_sigma(
        tallies.other_escape_count,
        N,
    )

    buildup = (
        T_total / T_unc
        if T_unc > 0
        else math.inf
    )

    sigma_buildup = (
        propagate_ratio_sigma(
            T_total,
            sigma_T_total,
            T_unc,
            sigma_T_unc,
        )
        if T_unc > 0
        else math.inf
    )

    # A zero transmitted count does not establish an infinite attenuation
    # coefficient; it censors this log-transmission diagnostic.  Keep it
    # explicitly unavailable so reports never print the misleading `inf ± inf`.
    mu_eff = (
        -math.log(T_total)
        / cfg.wall_thickness_cm
        if T_total > 0
        else math.nan
    )

    sigma_mu_eff = (
        sigma_T_total
        / (
            cfg.wall_thickness_cm
            * T_total
        )
        if T_total > 0
        else math.nan
    )

    energy_hist_sigma = np.sqrt(
        resultsafe_int(
            tallies.energy_hist_counts
        )
    )

    angle_hist_sigma = np.sqrt(
        resultsafe_int(
            tallies.exit_cosine_hist_counts
        )
    )

    scatter_hist_sigma = np.sqrt(
        resultsafe_int(
            tallies.scatter_hist
        )
    )

    first_interaction_sigma = np.sqrt(
        resultsafe_int(
            tallies.first_interaction_hist
        )
    )

    last_interaction_sigma = np.sqrt(
        resultsafe_int(
            tallies.last_interaction_hist
        )
    )

    exit_energy_angle_sigma = np.sqrt(
        resultsafe_int(
            tallies.exit_energy_angle_hist
        )
    )

    secondary_exit_hist_sigma = {
        pname: np.sqrt(
            resultsafe_int(arr)
        )
        for pname, arr
        in tallies.secondary_exit_hist_counts.items()
    }

    fluence_sem_cm = safe_sem_from_sum_sumsq(
        tallies.fluence_tracklen_cm,
        tallies.fluence_tracklen_cm_sumsq,
        N,
    )

    edep_sem_mev = safe_sem_from_sum_sumsq(
        tallies.edep_vs_depth_mev,
        tallies.edep_vs_depth_mev_sumsq,
        N,
    )

    deltaE_sem_mev = safe_sem_from_sum_sumsq(
        tallies.primary_gamma_deltaE_vs_depth_mev,
        tallies.primary_gamma_deltaE_vs_depth_mev_sumsq,
        N,
    )

    fluence_per_primary_per_cm = (
        tallies.fluence_tracklen_cm
        / (
            N
            * tallies.depth_bin_widths_cm
        )
    )

    fluence_per_primary_per_cm_sigma = (
        fluence_sem_cm
        / tallies.depth_bin_widths_cm
    )

    edep_per_primary_per_cm = (
        tallies.edep_vs_depth_mev
        / (
            N
            * tallies.depth_bin_widths_cm
        )
    )

    edep_per_primary_per_cm_sigma = (
        edep_sem_mev
        / tallies.depth_bin_widths_cm
    )

    deltaE_per_primary_per_cm = (
        tallies.primary_gamma_deltaE_vs_depth_mev
        / (
            N
            * tallies.depth_bin_widths_cm
        )
    )

    deltaE_per_primary_per_cm_sigma = (
        deltaE_sem_mev
        / tallies.depth_bin_widths_cm
    )

    density_g_cm3 = (
        tallies.material_density_g_cm3
    )

    bin_volume_cm3 = (
        cfg.wall_size_y_cm
        * cfg.wall_size_z_cm
        * tallies.depth_bin_widths_cm
    )

    bin_mass_kg = (
        density_g_cm3
        * bin_volume_cm3
        / 1000.0
    )

    edep_joule_vs_depth = (
        tallies.edep_vs_depth_mev
        * MEV_TO_J
    )

    dose_vs_depth_gy = (
        edep_joule_vs_depth
        / bin_mass_kg
    )

    dose_vs_depth_gy_sigma = (
        edep_sem_mev
        * MEV_TO_J
        / bin_mass_kg
    )

    total_wall_mass_kg = (
        density_g_cm3
        * cfg.wall_size_y_cm
        * cfg.wall_size_z_cm
        * cfg.wall_thickness_cm
        / 1000.0
    )

    total_edep_mev = float(
        np.sum(
            tallies.edep_vs_depth_mev
        )
    )

    total_edep_mev_sem = float(
        np.sqrt(
            np.sum(
                edep_sem_mev ** 2
            )
        )
    )

    total_dose_gy = (
        total_edep_mev
        * MEV_TO_J
        / total_wall_mass_kg
    )

    total_dose_gy_sigma = (
        total_edep_mev_sem
        * MEV_TO_J
        / total_wall_mass_kg
    )

    equivalent_dose_proxy_sv_wR1 = total_dose_gy
    equivalent_dose_proxy_sv_wR1_sigma = total_dose_gy_sigma

    equivalent_dose_proxy_vs_depth_sv_wR1 = (
        dose_vs_depth_gy.copy()
    )

    equivalent_dose_proxy_vs_depth_sv_wR1_sigma = (
        dose_vs_depth_gy_sigma.copy()
    )

    buildup_primary_fraction_of_transmitted = {}
    buildup_primary_fraction_of_all = {}
    buildup_primary_sigma_of_transmitted = {}
    buildup_primary_sigma_of_all = {}

    for (
        key,
        count,
    ) in tallies.buildup_primary_counts.items():
        buildup_primary_fraction_of_transmitted[key] = (
            count / tallies.transmitted_count
            if tallies.transmitted_count > 0
            else 0.0
        )

        buildup_primary_fraction_of_all[key] = (
            count / N
        )

        buildup_primary_sigma_of_transmitted[key] = (
            binomial_sigma(
                count,
                tallies.transmitted_count,
            )
            if tallies.transmitted_count > 0
            else 0.0
        )

        buildup_primary_sigma_of_all[key] = (
            binomial_sigma(
                count,
                N,
            )
        )

    return {
        'summary': {
            'beam_half_angle_deg': float(
                cfg.beam_half_angle_deg
            ),
            'source_to_wall_gap_cm': float(
                cfg.source_to_wall_gap_cm
            ),
            'beam_mev_input': float(
                cfg.mev
            ),
            'source_energy_mev_internal': float(
                source_energy_mev
            ),
            'source_mode': get_source_mode(
                cfg
            ),
            'spectrum_case': (
                get_spectrum_case(cfg)
                if get_source_mode(cfg)
                == 'digitized_spectrum'
                else None
            ),
            'spectrum_case_origin': (
                get_spectrum_origin(
                    get_spectrum_case(cfg)
                )
                if get_source_mode(cfg)
                == 'digitized_spectrum'
                else None
            ),
            'spectrum_provenance_note': (
                get_spectrum_provenance_note(get_spectrum_case(cfg))
                if get_source_mode(cfg) == 'digitized_spectrum' else None
            ),
            'spectrum_revision': (
                get_spectrum_revision(get_spectrum_case(cfg))
                if get_source_mode(cfg) == 'digitized_spectrum' else None
            ),
            'spectrum_array_sha256': (
                get_spectrum_array_sha256(get_spectrum_case(cfg))
                if get_source_mode(cfg) == 'digitized_spectrum' else None
            ),
            'spectrum_files_used': (
                list(SPECTRUM_FILES_USED)
                if get_source_mode(cfg)
                == 'digitized_spectrum'
                else []
            ),
            'source_energy_hist_max_mev': float(
                source_energy_hist_max_mev(
                    cfg
                )
            ),
            'wall_thickness_cm': float(
                cfg.wall_thickness_cm
            ),
            'num_photons': int(
                cfg.num_photons
            ),
            'n_visual_tracks_requested': int(
                cfg.n_visual_tracks
            ),
            'material_name_used': (
                tallies.material_name_used
            ),
            'material_density_g_cm3': (
                density_g_cm3
            ),
            'transmitted_count': int(
                tallies.transmitted_count
            ),
            'uncollided_transmitted_count': int(
                tallies.uncollided_transmitted_count
            ),
            'backscatter_count': int(
                tallies.backscatter_count
            ),
            'removed_count': int(
                tallies.removed_count
            ),
            'other_escape_count': int(
                tallies.other_escape_count
            ),
            'transmitted_fraction': float(
                T_total
            ),
            'transmitted_fraction_sigma': float(
                sigma_T_total
            ),
            'uncollided_transmitted_fraction': float(
                T_unc
            ),
            'uncollided_transmitted_fraction_sigma': float(
                sigma_T_unc
            ),
            'backscatter_fraction': float(
                backscatter_fraction
            ),
            'backscatter_fraction_sigma': float(
                sigma_backscatter
            ),
            'removed_fraction': float(
                removed_fraction
            ),
            'removed_fraction_sigma': float(
                sigma_removed
            ),
            'other_escape_fraction': float(
                other_escape_fraction
            ),
            'other_escape_fraction_sigma': float(
                sigma_other_escape
            ),
            'buildup_factor_estimate': float(
                buildup
            ),
            'buildup_factor_sigma': float(
                sigma_buildup
            ),
            'effective_mu_per_cm': float(
                mu_eff
            ),
            'effective_mu_per_cm_sigma': float(
                sigma_mu_eff
            ),
            'total_absorbed_dose_gy': float(
                total_dose_gy
            ),
            'total_absorbed_dose_gy_sigma': float(
                total_dose_gy_sigma
            ),
            'equivalent_dose_proxy_sv_wR1': float(
                equivalent_dose_proxy_sv_wR1
            ),
            'equivalent_dose_proxy_sv_wR1_sigma': float(
                equivalent_dose_proxy_sv_wR1_sigma
            ),
            'tracks_saved': int(
                tallies.tracks_saved
            ),
            'tracks_rejected_short': int(
                tallies.tracks_rejected_short
            ),
            'runtime_s': float(
                runtime_s
            ),
        },
        'process_counts_all_tracks_total': dict(
            tallies.process_counts_all_tracks
        ),
        'process_counts_primary_total': dict(
            tallies.process_counts_primary_total
        ),
        'process_counts_primary_eventwise': dict(
            tallies.process_counts_primary_eventwise
        ),
        'buildup_primary_counts': dict(
            tallies.buildup_primary_counts
        ),
        'buildup_primary_fraction_of_transmitted': (
            buildup_primary_fraction_of_transmitted
        ),
        'buildup_primary_fraction_of_transmitted_sigma': (
            buildup_primary_sigma_of_transmitted
        ),
        'buildup_primary_fraction_of_all': (
            buildup_primary_fraction_of_all
        ),
        'buildup_primary_fraction_of_all_sigma': (
            buildup_primary_sigma_of_all
        ),
        'depth_cm': (
            tallies.depth_mid_cm.copy()
        ),
        'depth_bin_widths_cm': (
            tallies.depth_bin_widths_cm.copy()
        ),
        'fluence_tracklen_cm': (
            tallies.fluence_tracklen_cm.copy()
        ),
        'fluence_tracklen_cm_sem': (
            fluence_sem_cm.copy()
        ),
        'fluence_per_primary_per_cm': (
            fluence_per_primary_per_cm.copy()
        ),
        'fluence_per_primary_per_cm_sigma': (
            fluence_per_primary_per_cm_sigma.copy()
        ),
        'edep_vs_depth_mev': (
            tallies.edep_vs_depth_mev.copy()
        ),
        'edep_vs_depth_mev_sem': (
            edep_sem_mev.copy()
        ),
        'edep_per_primary_per_cm': (
            edep_per_primary_per_cm.copy()
        ),
        'edep_per_primary_per_cm_sigma': (
            edep_per_primary_per_cm_sigma.copy()
        ),
        'primary_gamma_deltaE_vs_depth_mev': (
            tallies.primary_gamma_deltaE_vs_depth_mev.copy()
        ),
        'primary_gamma_deltaE_vs_depth_mev_sem': (
            deltaE_sem_mev.copy()
        ),
        'primary_gamma_deltaE_per_primary_per_cm': (
            deltaE_per_primary_per_cm.copy()
        ),
        'primary_gamma_deltaE_per_primary_per_cm_sigma': (
            deltaE_per_primary_per_cm_sigma.copy()
        ),
        'dose_vs_depth_gy': (
            dose_vs_depth_gy.copy()
        ),
        'dose_vs_depth_gy_sigma': (
            dose_vs_depth_gy_sigma.copy()
        ),
        'equivalent_dose_proxy_vs_depth_sv_wR1': (
            equivalent_dose_proxy_vs_depth_sv_wR1.copy()
        ),
        'equivalent_dose_proxy_vs_depth_sv_wR1_sigma': (
            equivalent_dose_proxy_vs_depth_sv_wR1_sigma.copy()
        ),
        'first_interaction_hist_counts': (
            tallies.first_interaction_hist.copy()
        ),
        'first_interaction_hist_sigma': (
            first_interaction_sigma.copy()
        ),
        'last_interaction_hist_counts': (
            tallies.last_interaction_hist.copy()
        ),
        'last_interaction_hist_sigma': (
            last_interaction_sigma.copy()
        ),
        'scatter_hist_counts': (
            tallies.scatter_hist.copy()
        ),
        'scatter_hist_sigma': (
            scatter_hist_sigma.copy()
        ),
        'scatter_hist_bins': np.arange(
            cfg.max_scatter_bin + 2,
            dtype=int,
        ),
        'energy_hist_edges_mev': (
            tallies.energy_hist_edges_mev.copy()
        ),
        'energy_hist_counts': (
            tallies.energy_hist_counts.copy()
        ),
        'energy_hist_sigma': (
            energy_hist_sigma.copy()
        ),
        'exit_cosine_hist_edges': (
            tallies.exit_cosine_hist_edges.copy()
        ),
        'exit_cosine_hist_counts': (
            tallies.exit_cosine_hist_counts.copy()
        ),
        'exit_cosine_hist_sigma': (
            angle_hist_sigma.copy()
        ),
        'exit_energy_angle_hist_counts': (
            tallies.exit_energy_angle_hist.copy()
        ),
        'exit_energy_angle_hist_sigma': (
            exit_energy_angle_sigma.copy()
        ),
        'exit_energy_angle_cos_edges': (
            tallies.corr_cos_edges.copy()
        ),
        'exit_energy_angle_energy_edges_mev': (
            tallies.corr_energy_edges_mev.copy()
        ),
        'transmitted_energy_sample_mev': np.array(
            tallies.transmitted_energy_sample_mev,
            dtype=float,
        ),
        'exit_cosine_sample': np.array(
            tallies.exit_cosine_sample,
            dtype=float,
        ),
        'secondary_energy_hist_edges_mev': (
            tallies.secondary_energy_hist_edges_mev.copy()
        ),
        'secondary_exit_hist_counts': {
            key: value.copy()
            for key, value
            in tallies.secondary_exit_hist_counts.items()
        },
        'secondary_exit_hist_sigma': {
            key: value.copy()
            for key, value
            in secondary_exit_hist_sigma.items()
        },
        'secondary_exit_counts_forward': dict(
            tallies.secondary_exit_counts_forward
        ),
        'tracks_cm': tallies.tracks_cm,
        'obj_export': {
            'enabled': bool(
                cfg.export_obj
            ),
            'filepath': (
                cfg.obj_filepath
                if cfg.export_obj
                else None
            ),
            'info': obj_info,
        },
    }


def print_summary(
    results: Dict[str, Any],
) -> None:
    s = results['summary']

    print(
        '\n==================== RESULTS ===================='
    )

    print(
        f"beam_half_angle_deg                    = "
        f"{s['beam_half_angle_deg']:.6g}"
    )

    print(
        f"source_to_wall_gap_cm                  = "
        f"{s['source_to_wall_gap_cm']:.6g}"
    )

    print(
        f"beam_mev_input                         = "
        f"{s['beam_mev_input']:.6g}"
    )

    print(
        f"source_energy_mev_internal             = "
        f"{s['source_energy_mev_internal']:.6g}"
    )

    print(
        f"source_mode                            = "
        f"{s['source_mode']}"
    )

    if s['spectrum_case'] is not None:
        print(
            f"spectrum_case                          = "
            f"{s['spectrum_case']}"
        )

    if s.get('spectrum_case_origin') is not None:
        print(
            f"spectrum_case_origin                   = "
            f"{s['spectrum_case_origin']}"
        )

    if s.get('spectrum_files_used'):
        print(
            f"spectrum_files_used                    = "
            f"{' | '.join(s['spectrum_files_used'])}"
        )

    print(
        f"source_energy_hist_max_mev             = "
        f"{s['source_energy_hist_max_mev']:.6g}"
    )

    print(
        f"wall_thickness_cm                      = "
        f"{s['wall_thickness_cm']:.6g}"
    )

    print(
        f"num_photons                            = "
        f"{s['num_photons']}"
    )

    print(
        f"material_name_used                     = "
        f"{s['material_name_used']}"
    )

    print(
        f"material_density_g_cm3                 = "
        f"{s['material_density_g_cm3']:.6g}"
    )

    print('')

    print(
        f"transmitted_count                      = "
        f"{s['transmitted_count']}"
    )

    print(
        f"uncollided_transmitted_count           = "
        f"{s['uncollided_transmitted_count']}"
    )

    print(
        f"backscatter_count                      = "
        f"{s['backscatter_count']}"
    )

    print(
        f"removed_count                          = "
        f"{s['removed_count']}"
    )

    print(
        f"other_escape_count                     = "
        f"{s['other_escape_count']}"
    )

    print('')

    print(
        f"transmitted_fraction                   = "
        f"{s['transmitted_fraction']:.6e} ± "
        f"{s['transmitted_fraction_sigma']:.2e}"
    )

    print(
        f"uncollided_transmitted_fraction        = "
        f"{s['uncollided_transmitted_fraction']:.6e} ± "
        f"{s['uncollided_transmitted_fraction_sigma']:.2e}"
    )

    print(
        f"backscatter_fraction                   = "
        f"{s['backscatter_fraction']:.6e} ± "
        f"{s['backscatter_fraction_sigma']:.2e}"
    )

    print(
        f"removed_fraction                       = "
        f"{s['removed_fraction']:.6e} ± "
        f"{s['removed_fraction_sigma']:.2e}"
    )

    print(
        f"other_escape_fraction                  = "
        f"{s['other_escape_fraction']:.6e} ± "
        f"{s['other_escape_fraction_sigma']:.2e}"
    )

    print(
        f"buildup_factor_estimate                = "
        f"{s['buildup_factor_estimate']:.6e} ± "
        f"{s['buildup_factor_sigma']:.2e}"
    )

    if np.isfinite(float(s.get('effective_mu_per_cm', np.nan))) and np.isfinite(float(s.get('effective_mu_per_cm_sigma', np.nan))):
        print(
            f"effective_mu_per_cm                    = "
            f"{s['effective_mu_per_cm']:.6e} ± "
            f"{s['effective_mu_per_cm_sigma']:.2e}"
        )
    else:
        print(
            "effective_mu_per_cm                    = NOT ESTABLISHED — "
            "zero/sparse transmitted-primary sample censors the log-transmission diagnostic"
        )

    print('')

    print(
        f"total_absorbed_dose_gy                 = "
        f"{s['total_absorbed_dose_gy']:.6e} ± "
        f"{s['total_absorbed_dose_gy_sigma']:.2e}"
    )

    print(
        f"equivalent_dose_proxy_sv_wR1           = "
        f"{s['equivalent_dose_proxy_sv_wR1']:.6e} ± "
        f"{s['equivalent_dose_proxy_sv_wR1_sigma']:.2e}"
    )

    print('')

    print(
        'process_counts_primary_total           =',
        results['process_counts_primary_total'],
    )

    print(
        'process_counts_primary_eventwise       =',
        results['process_counts_primary_eventwise'],
    )

    print(
        'process_counts_all_tracks_total        =',
        results['process_counts_all_tracks_total'],
    )

    print('')

    print(
        'buildup_primary_counts                 =',
        results['buildup_primary_counts'],
    )

    print(
        'secondary_exit_counts_forward          =',
        results['secondary_exit_counts_forward'],
    )

    print('')

    print(
        f"tracks_saved                           = "
        f"{s['tracks_saved']}"
    )

    print(
        f"tracks_rejected_short                  = "
        f"{s['tracks_rejected_short']}"
    )

    print(
        f"runtime_s                              = "
        f"{s['runtime_s']:.3f}"
    )

    print(
        '=================================================\n'
    )


def run_simulation(
    cfg,
) -> Dict[str, Any]:
    _LIVE_OBJECTS.clear()

    tallies = SimulationTallies(cfg)
    state = EventState()

    detector = MyDetectorConstruction(
        cfg,
        tallies,
    )

    physics = MyPhysicsList(cfg)

    actions = MyActionInitialization(
        cfg,
        tallies,
        state,
    )

    run_manager = G4RunManager()

    _LIVE_OBJECTS.extend(
        [
            tallies,
            state,
            detector,
            physics,
            actions,
            run_manager,
        ]
    )

    run_manager.SetUserInitialization(
        detector
    )

    run_manager.SetUserInitialization(
        physics
    )

    run_manager.SetUserInitialization(
        actions
    )

    run_manager.Initialize()

    ui = G4UImanager.GetUIpointer()

    ui.ApplyCommand(
        f'/control/verbose '
        f'{cfg.geant4_control_verbose}'
    )

    ui.ApplyCommand(
        f'/run/verbose '
        f'{cfg.geant4_run_verbose}'
    )

    ui.ApplyCommand(
        f'/event/verbose '
        f'{cfg.geant4_event_verbose}'
    )

    ui.ApplyCommand(
        f'/tracking/verbose '
        f'{cfg.geant4_tracking_verbose}'
    )

    t0 = time.perf_counter()

    run_manager.BeamOn(
        int(cfg.num_photons)
    )

    runtime_s = (
        time.perf_counter()
        - t0
    )

    obj_info = None

    if (
        cfg.export_obj
        and (
            cfg.include_wall_in_obj
            or len(tallies.tracks_cm) > 0
        )
    ):
        obj_info = export_scene_to_obj(
            tracks_cm=tallies.tracks_cm,
            filepath=cfg.obj_filepath,
            scale=cfg.obj_scale,
            decimate=cfg.obj_decimate,
            wall_thickness_cm=cfg.wall_thickness_cm,
            wall_size_y_cm=cfg.wall_size_y_cm,
            wall_size_z_cm=cfg.wall_size_z_cm,
            include_wall=cfg.include_wall_in_obj,
            track_radius_cm=cfg.obj_track_radius_cm,
        )

    results = build_results(
        cfg,
        tallies,
        runtime_s,
        obj_info,
    )

    if cfg.verbose:
        print_summary(results)

        if obj_info is not None:
            print(
                f'[OBJ EXPORT] '
                f'{cfg.obj_filepath}'
            )

            print(
                f"  mtl_filepath           = "
                f"{obj_info['mtl_filepath']}"
            )

            print(
                f"  tracks_written         = "
                f"{obj_info['tracks_written']}"
            )

            print(
                f"  segment_prisms_written = "
                f"{obj_info['segment_prisms_written']}"
            )

            print(
                f"  vertices               = "
                f"{obj_info['vertices']}"
            )

            print(
                f"  faces                  = "
                f"{obj_info['faces']}"
            )

            print(
                f"  wall_included          = "
                f"{obj_info['wall_included']}"
            )

            print(
                f"  track_radius_cm        = "
                f"{obj_info['track_radius_cm']}"
            )

            print(
                f'  scale                  = '
                f'{cfg.obj_scale}'
            )

            print(
                f'  decimate               = '
                f'{cfg.obj_decimate}'
            )

            print('')

    return results



def run_thickness_sweep(cfg) -> Dict[str, np.ndarray]:
    """Run actual independent Geant4 thickness points and return v7 physical + diagnostic response arrays."""
    print('\n==================== CONCRETE THICKNESS SWEEP ====================')
    print('Each point is a separate full Geant4 run.')
    print('Physical shielding response: all forward-crossing photons, dry-air energy-absorption weighted.')
    print('Primary-history transmission remains a diagnostic only.')
    print('=================================================================\n')

    thicknesses = np.asarray(cfg.sweep_thicknesses_cm, dtype=float)
    if thicknesses.size == 0:
        raise ValueError('sweep_thicknesses_cm must contain at least one thickness')
    n_point = int(cfg.sweep_num_photons_per_point)
    if n_point <= 0:
        raise ValueError('sweep_num_photons_per_point must be positive')

    transmitted = []
    transmitted_sigma = []
    mu_eff = []
    mu_eff_sigma = []
    all_response = []
    all_response_sem = []
    all_crossings = []
    source_weight_mean = []
    max_weight = []
    primaries = []
    subrun_ids = []
    subrun_purposes = []
    subrun_seeds = []

    for index, thickness in enumerate(thicknesses, start=1):
        print(
            f'\\n===== CONCRETE SWEEP POINT {index}/{len(thicknesses)}: '
            f'{float(thickness):g} cm, N={n_point:,} =====',
            flush=True,
        )
        point_cfg = deepcopy(cfg)
        point_cfg.wall_thickness_cm = float(thickness)
        point_cfg.num_photons = n_run
        point_cfg.n_visual_tracks = 0
        point_cfg.export_obj = False
        point_cfg.verbose = False
        point_cfg.run_thickness_sweep = False
        parent_run_id = str(getattr(cfg, 'ncrp_run_id', 'run-unavailable'))
        point_cfg.transport_subrun_id = f"{parent_run_id}:photon-sweep:{index:02d}:{float(thickness):g}cm"
        point_cfg.transport_subrun_purpose = "adaptive photon sweep sampled point"
        point_result = run_simulation(point_cfg)
        summary = point_result.get('summary', {})

        transmitted.append(float(summary.get('transmitted_fraction', np.nan)))
        transmitted_sigma.append(float(summary.get('transmitted_fraction_sigma', np.nan)))
        mu_eff.append(float(summary.get('effective_mu_per_cm', np.nan)))
        mu_eff_sigma.append(float(summary.get('effective_mu_per_cm_sigma', np.nan)))
        all_response.append(float(summary.get('all_photon_air_kerma_response', np.nan)))
        all_response_sem.append(float(summary.get('all_photon_air_kerma_response_sem', np.nan)))
        all_crossings.append(int(summary.get('all_photon_forward_crossing_count', 0) or 0))
        source_weight_mean.append(float(summary.get('source_air_kerma_weight_mean', np.nan)))
        max_weight.append(float(summary.get('max_air_kerma_weight_per_photon', np.nan)))
        primaries.append(n_point)
        subrun_ids.append(str(summary.get('transport_subrun_id', point_cfg.transport_subrun_id)))
        subrun_purposes.append(str(summary.get('transport_subrun_purpose', point_cfg.transport_subrun_purpose)))
        subrun_seeds.append(int(summary.get('random_seed', getattr(point_cfg, 'random_seed', 0) or 0)))

        print(
            f'SWEEP RESULT: thickness={float(thickness):g} cm | '
            f'all-photon dry-air energy-absorption response={all_response[-1]:.8e} | '
            f'SEM={all_response_sem[-1]:.3e} | forward gamma crossings={all_crossings[-1]:,} | '
            f'primary-history T={transmitted[-1]:.8e}',
            flush=True,
        )

    return {
        'thicknesses_cm': thicknesses,
        'primaries': np.asarray(primaries, dtype=np.int64),
        'transmitted_fraction': np.asarray(transmitted, dtype=float),
        'transmitted_fraction_sigma': np.asarray(transmitted_sigma, dtype=float),
        'effective_mu_per_cm': np.asarray(mu_eff, dtype=float),
        'effective_mu_per_cm_sigma': np.asarray(mu_eff_sigma, dtype=float),
        'all_photon_air_kerma_response': np.asarray(all_response, dtype=float),
        'all_photon_air_kerma_response_sem': np.asarray(all_response_sem, dtype=float),
        'all_photon_forward_crossing_count': np.asarray(all_crossings, dtype=np.int64),
        'source_air_kerma_weight_mean': np.asarray(source_weight_mean, dtype=float),
        'max_air_kerma_weight_per_photon': np.asarray(max_weight, dtype=float),
        'physical_selection_tally': 'all_photon_downstream_dry_air_kerma_response',
        'primary_history_fraction_role': 'diagnostic only',
        'transport_subrun_ids': np.asarray(subrun_ids, dtype=object),
        'transport_subrun_purposes': np.asarray(subrun_purposes, dtype=object),
        'transport_subrun_seeds': np.asarray(subrun_seeds, dtype=np.int64),
    }



def plot_results(
    results: Dict[str, Any],
    cfg,
    sweep_results: Optional[
        Dict[str, np.ndarray]
    ] = None,
) -> None:
    os.makedirs(
        cfg.plot_dir,
        exist_ok=True,
    )

    plot_counter = {
        'i': 0,
    }

    def show_fig(
        fig,
        name: str,
    ):
        safe_name = (
            name.replace(
                ' ',
                '_',
            )
            .replace(
                '/',
                '_',
            )
            .replace(
                '(',
                '',
            )
            .replace(
                ')',
                '',
            )
            .replace(
                ',',
                '',
            )
        )

        filepath = os.path.join(
            cfg.plot_dir,
            f"{plot_counter['i']:02d}_{safe_name}.png",
        )

        style_black_red_figure(fig)

        # 2026-08-10-coupled-ncrp-shielding-v3: core spectrum save trailing-zero trim
        try:
            from shielding_logic import trim_spectrum_figure_trailing_zero_tail
            trim_spectrum_figure_trailing_zero_tail(fig, name)
        except Exception as spectrum_trim_exc:
            print(f"Spectrum tail trim warning: {type(spectrum_trim_exc).__name__}: {spectrum_trim_exc}", flush=True)

        fig.tight_layout()

        fig.savefig(
            filepath,
            dpi=cfg.plot_dpi,
            bbox_inches='tight',
            facecolor=PLOT_BLACK,
            edgecolor=PLOT_BLACK,
        )

        plot_counter['i'] += 1
        plt.show()

    mu_ratio_E = getattr(
        cfg,
        'mu_ratio_E_tab',
        None,
    )

    mu_ratio_vals = getattr(
        cfg,
        'mu_ratio_values',
        None,
    )

    if (
        mu_ratio_E is None
        and mu_ratio_vals is None
    ):
        if (
            'E_tab' in globals()
            and 'muC_over_muT_tab' in globals()
        ):
            mu_ratio_E = np.asarray(
                globals()['E_tab'],
                dtype=float,
            )

            mu_ratio_vals = np.asarray(
                globals()['muC_over_muT_tab'],
                dtype=float,
            )

    depth_cm = results['depth_cm']

    raw_E_exit = np.asarray(
        results[
            'transmitted_energy_sample_mev'
        ],
        dtype=float,
    )

    raw_exit_cos = np.asarray(
        results['exit_cosine_sample'],
        dtype=float,
    )

    energy_edges = results[
        'energy_hist_edges_mev'
    ]

    energy_counts = results[
        'energy_hist_counts'
    ]

    energy_sigma = results.get(
        'energy_hist_sigma',
        np.sqrt(
            np.maximum(
                energy_counts,
                0,
            )
        ),
    )

    angle_edges = results[
        'exit_cosine_hist_edges'
    ]

    angle_counts = results[
        'exit_cosine_hist_counts'
    ]

    angle_sigma = results.get(
        'exit_cosine_hist_sigma',
        np.sqrt(
            np.maximum(
                angle_counts,
                0,
            )
        ),
    )

    scatter_hist = results[
        'scatter_hist_counts'
    ]

    scatter_sigma = results.get(
        'scatter_hist_sigma',
        np.sqrt(
            np.maximum(
                scatter_hist,
                0,
            )
        ),
    )

    max_scatter_bin = (
        len(scatter_hist) - 2
    )

    if (
        mu_ratio_E is not None
        and mu_ratio_vals is not None
    ):
        fig, ax = plt.subplots()

        ax.plot(
            mu_ratio_E,
            mu_ratio_vals,
            color=PLOT_RED,
            linewidth=PLOT_LINE_WIDTH,
        )

        ax.set_xlabel(
            'Energy (MeV)'
        )

        ax.set_ylabel(
            'muC / muT'
        )

        ax.set_title(
            'muC/muT vs energy (from table)'
        )

        ax.set_xscale(
            'log'
        )

        ax.set_ylim(
            0,
            1.05,
        )

        show_fig(
            fig,
            'muC_over_muT_vs_energy',
        )

    else:
        print(
            '[plot 0 skipped] muC/muT vs energy requires '
            'E_tab and muC_over_muT_tab data from your '
            'custom/XCOM table.'
        )

    if get_source_mode(cfg) == 'digitized_spectrum':
        spec = get_digitized_spectrum(
            cfg
        )

        fig, ax = plt.subplots()

        ax.bar(
            spec['energy_mev_center'],
            spec['probability_mass_bin'],
            width=(
                float(
                    spec.get(
                        'bin_width_mev',
                        0.1,
                    )
                )
                * 0.95
            ),
            align='center',
            color=PLOT_RED,
            edgecolor=PLOT_RED,
            linewidth=PLOT_EDGE_WIDTH,
        )

        ax.set_xlabel(
            'Source photon energy (MeV)'
        )

        ax.set_ylabel(
            'Probability per bin'
        )

        ax.set_title(
            f'Source spectrum used '
            f'({get_spectrum_case(cfg)})'
        )

        show_fig(
            fig,
            'source_spectrum_used',
        )

    if len(raw_E_exit) > 0:
        fig, ax = plt.subplots()

        ax.hist(
            raw_E_exit,
            bins=cfg.energy_hist_bins,
            color=PLOT_RED,
            edgecolor=PLOT_RED,
            linewidth=PLOT_EDGE_WIDTH,
        )

        ax.set_xlabel(
            'Exit photon energy (MeV)'
        )

        ax.set_ylabel(
            'Counts'
        )

        ax.set_title(
            'Transmitted Energy Spectrum'
        )

        show_fig(
            fig,
            'transmitted_energy_spectrum',
        )

    if len(raw_exit_cos) > 0:
        fig, ax = plt.subplots()

        ax.hist(
            raw_exit_cos,
            bins=cfg.angle_hist_bins,
            color=PLOT_RED,
            edgecolor=PLOT_RED,
            linewidth=PLOT_EDGE_WIDTH,
        )

        ax.set_xlabel(
            'cos(theta) relative to slab normal (+x)'
        )

        ax.set_ylabel(
            'Counts'
        )

        ax.set_title(
            'Angular Distribution of Transmitted Photons'
        )

        show_fig(
            fig,
            'angular_distribution_of_transmitted_photons',
        )

    if np.sum(scatter_hist) > 0:
        fig, ax = plt.subplots()

        x_bins = np.arange(
            len(scatter_hist)
        )

        ax.bar(
            x_bins,
            scatter_hist,
            color=PLOT_RED,
            edgecolor=PLOT_RED,
            linewidth=PLOT_EDGE_WIDTH,
        )

        ax.set_xlabel(
            f'# Compton scatters before transmission '
            f'(last bin >= {max_scatter_bin})'
        )

        ax.set_ylabel(
            'Counts'
        )

        ax.set_title(
            'Scatter-order distribution (transmitted)'
        )

        show_fig(
            fig,
            'scatter_order_distribution',
        )

    fig, ax = plt.subplots()

    phi = results[
        'fluence_per_primary_per_cm'
    ]

    ax.plot(
        depth_cm,
        phi,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Track-length fluence proxy '
        '(per primary per cm)'
    )

    ax.set_title(
        f"Fluence vs depth "
        f"(N={results['summary']['num_photons']})"
    )

    show_fig(
        fig,
        'fluence_vs_depth',
    )

    fig, ax = plt.subplots()

    dep_per_cm = results[
        'primary_gamma_deltaE_per_primary_per_cm'
    ]

    ax.plot(
        depth_cm,
        dep_per_cm,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Compton ΔE proxy '
        '(MeV per primary per cm)'
    )

    ax.set_title(
        f"Energy deposition proxy vs depth "
        f"(N={results['summary']['num_photons']})"
    )

    show_fig(
        fig,
        'deposition_proxy_vs_depth',
    )

    if (
        len(raw_E_exit) > 0
        and len(raw_exit_cos) > 0
    ):
        fig, ax = plt.subplots()

        hb = ax.hexbin(
            raw_exit_cos,
            raw_E_exit,
            gridsize=40,
            bins='log',
            cmap=RED_ORANGE_YELLOW_CMAP,
            mincnt=1,
        )

        ax.set_xlabel(
            'Exit cos(angle) relative to +x (u_x)'
        )

        ax.set_ylabel(
            'Exit energy (MeV)'
        )

        ax.set_title(
            'Exit energy vs exit angle '
            '(hexbin density)'
        )

        fig.colorbar(
            hb,
            ax=ax,
            label='log10(count)',
        )

        show_fig(
            fig,
            'exit_energy_vs_exit_angle_hexbin',
        )

    if sweep_results is not None:
        xs = np.asarray(sweep_results['thicknesses_cm'], dtype=float)
        transmissions = np.asarray(sweep_results.get('transmitted_fraction', []), dtype=float)
        upper95 = np.asarray(sweep_results.get('transmitted_fraction_wilson95_high', []), dtype=float)
        target = float(getattr(cfg, 'target_transmission', math.nan))
        if transmissions.size == xs.size and upper95.size == xs.size:
            design95 = _suffix_maximum(upper95)
            fig, ax = plt.subplots()
            ax.semilogy(xs, np.maximum(transmissions, np.finfo(float).tiny), marker='o', label='Observed transmission')
            ax.semilogy(xs, np.maximum(upper95, np.finfo(float).tiny), marker='o', label='95% Wilson upper')
            ax.semilogy(xs, np.maximum(design95, np.finfo(float).tiny), linewidth=max(PLOT_LINE_WIDTH, 0.8), label='Conservative design envelope')
            if math.isfinite(target) and target > 0.0:
                ax.axhline(target, linestyle='--', linewidth=max(PLOT_LINE_WIDTH, 0.8), label='Target transmission')
            ax.set_xlabel('Concrete thickness (cm)')
            ax.set_ylabel('Transmitted photon fraction')
            ax.set_title(f'Photon Transmission vs Concrete Thickness (N={int(cfg.sweep_num_photons_per_point):,} per point)')
            ax.grid(True, which='both', alpha=0.22)
            ax.legend()
            show_fig(fig, 'transmission_vs_concrete_thickness')

        mu_effs = sweep_results[
            'effective_mu_per_cm'
        ]

        fig, ax = plt.subplots()

        ax.plot(
            xs,
            mu_effs,
            marker='o',
            color=PLOT_RED,
            markerfacecolor=PLOT_RED,
            markeredgecolor=PLOT_RED,
            markeredgewidth=PLOT_MARKER_EDGE_WIDTH,
            linewidth=PLOT_LINE_WIDTH,
        )

        ax.set_xlabel(
            'Thickness d (cm)'
        )

        ax.set_ylabel(
            'mu_eff = -ln(T_total)/d  (1/cm)'
        )

        ax.set_title(
            'Effective mu vs thickness '
            '(stage-specific per-point N)'
        )

        ax.set_ylim(
            bottom=0,
        )

        show_fig(
            fig,
            'effective_mu_vs_thickness',
        )

    tracks = results[
        'tracks_cm'
    ]

    if len(tracks) > 0:
        fig = plt.figure()

        ax = fig.add_subplot(
            111,
            projection='3d',
        )

        for (
            track_index,
            (
                xs,
                ys,
                zs,
            ),
        ) in enumerate(tracks):
            xs = np.asarray(
                xs,
                dtype=float,
            )

            ys = np.asarray(
                ys,
                dtype=float,
            )

            zs = np.asarray(
                zs,
                dtype=float,
            )

            ax.plot(
                xs[
                    ::max(
                        1,
                        cfg.plot_track_decimate,
                    )
                ],
                ys[
                    ::max(
                        1,
                        cfg.plot_track_decimate,
                    )
                ],
                zs[
                    ::max(
                        1,
                        cfg.plot_track_decimate,
                    )
                ],
                linewidth=PLOT_TRACK_LINE_WIDTH,
                color=(
                    PLOT_RED
                    if track_index % 2 == 0
                    else PLOT_ORANGE
                ),
            )

        ax.set_xlabel(
            'Depth x (cm)'
        )

        ax.set_ylabel(
            'y (cm)'
        )

        ax.set_zlabel(
            'z (cm)'
        )

        ax.set_title(
            f'3D Photon Tracks '
            f'(accepted={len(tracks)}, '
            f'plot-decimate='
            f'{max(1, cfg.plot_track_decimate)})'
        )

        show_fig(
            fig,
            '3d_photon_tracks',
        )

    fig, ax = plt.subplots()

    y = results[
        'edep_per_primary_per_cm'
    ]

    ysig = results[
        'edep_per_primary_per_cm_sigma'
    ]

    ax.plot(
        depth_cm,
        y,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
        label='Mean',
    )

    ax.fill_between(
        depth_cm,
        np.maximum(
            y - ysig,
            0.0,
        ),
        y + ysig,
        color=PLOT_RED,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
        label='±1σ SEM',
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Deposited energy '
        '(MeV per primary per cm)'
    )

    ax.set_title(
        'Actual Deposited Energy vs Depth'
    )

    ax.legend()

    show_fig(
        fig,
        'actual_deposited_energy_vs_depth',
    )

    fig, ax = plt.subplots()

    y = results[
        'dose_vs_depth_gy'
    ]

    ysig = results[
        'dose_vs_depth_gy_sigma'
    ]

    ax.plot(
        depth_cm,
        y,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
        label='Mean',
    )

    ax.fill_between(
        depth_cm,
        np.maximum(
            y - ysig,
            0.0,
        ),
        y + ysig,
        color=PLOT_RED,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
        label='±1σ SEM',
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Dose (Gy)'
    )

    ax.set_title(
        'Dose vs Depth'
    )

    ax.legend()

    show_fig(
        fig,
        'dose_vs_depth',
    )

    fig, ax = plt.subplots()

    y = results[
        'equivalent_dose_proxy_vs_depth_sv_wR1'
    ]

    ysig = results[
        'equivalent_dose_proxy_vs_depth_sv_wR1_sigma'
    ]

    ax.plot(
        depth_cm,
        y,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
        label='Mean',
    )

    ax.fill_between(
        depth_cm,
        np.maximum(
            y - ysig,
            0.0,
        ),
        y + ysig,
        color=PLOT_RED,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
        label='±1σ SEM',
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Equivalent dose proxy (Sv, wR=1)'
    )

    ax.set_title(
        'Equivalent Dose Proxy vs Depth'
    )

    ax.legend()

    show_fig(
        fig,
        'equivalent_dose_proxy_vs_depth',
    )

    fig, ax = plt.subplots()

    first_counts = results[
        'first_interaction_hist_counts'
    ]

    first_sigma = results[
        'first_interaction_hist_sigma'
    ]

    last_counts = results[
        'last_interaction_hist_counts'
    ]

    last_sigma = results[
        'last_interaction_hist_sigma'
    ]

    ax.plot(
        depth_cm,
        first_counts,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
        label='First interaction',
    )

    ax.fill_between(
        depth_cm,
        np.maximum(
            first_counts - first_sigma,
            0.0,
        ),
        first_counts + first_sigma,
        color=PLOT_RED,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
    )

    ax.plot(
        depth_cm,
        last_counts,
        color=PLOT_ORANGE,
        linewidth=PLOT_LINE_WIDTH,
        label='Last interaction',
    )

    ax.fill_between(
        depth_cm,
        np.maximum(
            last_counts - last_sigma,
            0.0,
        ),
        last_counts + last_sigma,
        color=PLOT_ORANGE,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Counts'
    )

    ax.set_title(
        'Interaction-Depth Distributions'
    )

    ax.legend()

    show_fig(
        fig,
        'interaction_depth_distributions',
    )

    summary = results[
        'summary'
    ]

    labels = [
        'Forward transmit',
        'Backscatter',
        'Removed',
        'Other escape',
    ]

    vals = [
        summary['transmitted_fraction'],
        summary['backscatter_fraction'],
        summary['removed_fraction'],
        summary['other_escape_fraction'],
    ]

    errs = [
        summary['transmitted_fraction_sigma'],
        summary['backscatter_fraction_sigma'],
        summary['removed_fraction_sigma'],
        summary['other_escape_fraction_sigma'],
    ]

    fig, ax = plt.subplots()

    ax.bar(
        labels,
        vals,
        yerr=errs,
        capsize=4,
        color=alternating_red_orange_colors(
            len(labels)
        ),
        edgecolor=alternating_red_orange_colors(
            len(labels)
        ),
        linewidth=PLOT_EDGE_WIDTH,
        error_kw={
            'ecolor': PLOT_RED,
            'elinewidth': PLOT_ERROR_LINE_WIDTH,
            'capthick': PLOT_ERROR_LINE_WIDTH,
        },
    )

    ax.set_ylabel(
        'Fraction of primaries'
    )

    ax.set_title(
        'Forward / Backscatter / Removal Partition'
    )

    ax.tick_params(
        axis='x',
        rotation=20,
    )

    show_fig(
        fig,
        'partition_fractions',
    )

    fig, ax = plt.subplots()

    e_centers = 0.5 * (
        energy_edges[:-1]
        + energy_edges[1:]
    )

    e_widths = np.diff(
        energy_edges
    )

    ax.bar(
        e_centers,
        energy_counts,
        width=e_widths,
        align='center',
        color=PLOT_RED,
        edgecolor=PLOT_RED,
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.errorbar(
        e_centers,
        energy_counts,
        yerr=energy_sigma,
        fmt='none',
        capsize=2,
        ecolor=PLOT_RED,
        elinewidth=PLOT_ERROR_LINE_WIDTH,
        capthick=PLOT_ERROR_LINE_WIDTH,
    )

    ax.set_xlabel(
        'Exit photon energy (MeV)'
    )

    ax.set_ylabel(
        'Counts'
    )

    ax.set_title(
        'Transmitted Energy Spectrum '
        'with Statistical Error Bars'
    )

    show_fig(
        fig,
        'transmitted_energy_spectrum_with_error_bars',
    )

    fig, ax = plt.subplots()

    a_centers = 0.5 * (
        angle_edges[:-1]
        + angle_edges[1:]
    )

    a_widths = np.diff(
        angle_edges
    )

    ax.bar(
        a_centers,
        angle_counts,
        width=a_widths,
        align='center',
        color=PLOT_RED,
        edgecolor=PLOT_RED,
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.errorbar(
        a_centers,
        angle_counts,
        yerr=angle_sigma,
        fmt='none',
        capsize=2,
        ecolor=PLOT_RED,
        elinewidth=PLOT_ERROR_LINE_WIDTH,
        capthick=PLOT_ERROR_LINE_WIDTH,
    )

    ax.set_xlabel(
        'cos(theta) relative to +x'
    )

    ax.set_ylabel(
        'Counts'
    )

    ax.set_title(
        'Angular Distribution '
        'with Statistical Error Bars'
    )

    show_fig(
        fig,
        'angular_distribution_with_error_bars',
    )

    fig, ax = plt.subplots()

    scatter_x = np.arange(
        len(scatter_hist)
    )

    ax.bar(
        scatter_x,
        scatter_hist,
        color=PLOT_RED,
        edgecolor=PLOT_RED,
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.errorbar(
        scatter_x,
        scatter_hist,
        yerr=scatter_sigma,
        fmt='none',
        capsize=2,
        ecolor=PLOT_RED,
        elinewidth=PLOT_ERROR_LINE_WIDTH,
        capthick=PLOT_ERROR_LINE_WIDTH,
    )

    ax.set_xlabel(
        f'# Compton scatters before transmission '
        f'(last bin >= {max_scatter_bin})'
    )

    ax.set_ylabel(
        'Counts'
    )

    ax.set_title(
        'Scatter-Order Distribution '
        'with Statistical Error Bars'
    )

    show_fig(
        fig,
        'scatter_order_distribution_with_error_bars',
    )

    proc_counts = results[
        'process_counts_primary_total'
    ]

    proc_names = list(
        proc_counts.keys()
    )

    proc_vals = [
        proc_counts[key]
        for key in proc_names
    ]

    fig, ax = plt.subplots()

    ax.bar(
        proc_names,
        proc_vals,
        color=alternating_red_orange_colors(
            len(proc_names)
        ),
        edgecolor=alternating_red_orange_colors(
            len(proc_names)
        ),
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.set_xlabel(
        'Process'
    )

    ax.set_ylabel(
        'Count on primary gamma'
    )

    ax.set_title(
        'Primary-Gamma Process Counts'
    )

    show_fig(
        fig,
        'primary_gamma_process_counts',
    )

    proc_evt = results[
        'process_counts_primary_eventwise'
    ]

    proc_names = list(
        proc_evt.keys()
    )

    proc_evt_vals = [
        proc_evt[key]
        for key in proc_names
    ]

    fig, ax = plt.subplots()

    ax.bar(
        proc_names,
        proc_evt_vals,
        color=alternating_red_orange_colors(
            len(proc_names)
        ),
        edgecolor=alternating_red_orange_colors(
            len(proc_names)
        ),
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.set_xlabel(
        'Process'
    )

    ax.set_ylabel(
        'Number of events with ≥1 occurrence'
    )

    ax.set_title(
        'Primary-Gamma Eventwise Process Incidence'
    )

    show_fig(
        fig,
        'primary_gamma_eventwise_process_incidence',
    )

    proc_all = results[
        'process_counts_all_tracks_total'
    ]

    proc_names = list(
        proc_all.keys()
    )

    proc_vals = [
        proc_all[key]
        for key in proc_names
    ]

    fig, ax = plt.subplots()

    ax.bar(
        proc_names,
        proc_vals,
        color=alternating_red_orange_colors(
            len(proc_names)
        ),
        edgecolor=alternating_red_orange_colors(
            len(proc_names)
        ),
        linewidth=PLOT_EDGE_WIDTH,
    )

    ax.set_xlabel(
        'Process'
    )

    ax.set_ylabel(
        'Count across all tracks'
    )

    ax.set_title(
        'All-Track Photon Process Counts'
    )

    show_fig(
        fig,
        'all_track_photon_process_counts',
    )

    buildup_counts = results[
        'buildup_primary_counts'
    ]

    buildup_labels = list(
        buildup_counts.keys()
    )

    buildup_vals = [
        results[
            'buildup_primary_fraction_of_transmitted'
        ][key]
        for key in buildup_labels
    ]

    buildup_errs = [
        results[
            'buildup_primary_fraction_of_transmitted_sigma'
        ][key]
        for key in buildup_labels
    ]

    fig, ax = plt.subplots()

    ax.bar(
        buildup_labels,
        buildup_vals,
        yerr=buildup_errs,
        capsize=4,
        color=alternating_red_orange_colors(
            len(buildup_labels)
        ),
        edgecolor=alternating_red_orange_colors(
            len(buildup_labels)
        ),
        linewidth=PLOT_EDGE_WIDTH,
        error_kw={
            'ecolor': PLOT_RED,
            'elinewidth': PLOT_ERROR_LINE_WIDTH,
            'capthick': PLOT_ERROR_LINE_WIDTH,
        },
    )

    ax.set_ylabel(
        'Fraction of transmitted primaries'
    )

    ax.set_title(
        'Primary-Gamma Buildup Decomposition'
    )

    ax.tick_params(
        axis='x',
        rotation=20,
    )

    show_fig(
        fig,
        'primary_gamma_buildup_decomposition',
    )

    sec_edges = results[
        'secondary_energy_hist_edges_mev'
    ]

    sec_centers = 0.5 * (
        sec_edges[:-1]
        + sec_edges[1:]
    )

    sec_widths = np.diff(
        sec_edges
    )

    for pname in [
        'gamma',
        'e-',
        'e+',
    ]:
        counts = results[
            'secondary_exit_hist_counts'
        ][pname]

        sigma = results[
            'secondary_exit_hist_sigma'
        ][pname]

        if np.sum(counts) == 0:
            print(
                f'[plot skipped] No forward-exit '
                f'secondary {pname} events recorded.'
            )

            continue

        fig, ax = plt.subplots()

        ax.bar(
            sec_centers,
            counts,
            width=sec_widths,
            align='center',
            color=PLOT_RED,
            edgecolor=PLOT_RED,
            linewidth=PLOT_EDGE_WIDTH,
        )

        ax.errorbar(
            sec_centers,
            counts,
            yerr=sigma,
            fmt='none',
            capsize=2,
            ecolor=PLOT_RED,
            elinewidth=PLOT_ERROR_LINE_WIDTH,
            capthick=PLOT_ERROR_LINE_WIDTH,
        )

        ax.set_xlabel(
            f'Forward-exit secondary '
            f'{pname} energy (MeV)'
        )

        ax.set_ylabel(
            'Counts'
        )

        ax.set_title(
            f'Secondary {pname} Exit Spectrum'
        )

        safe_particle_name = (
            pname.replace(
                '+',
                'plus',
            ).replace(
                '-',
                'minus',
            )
        )

        show_fig(
            fig,
            f'secondary_{safe_particle_name}_exit_spectrum',
        )

    corr = results[
        'exit_energy_angle_hist_counts'
    ].T

    xedges = results[
        'exit_energy_angle_cos_edges'
    ]

    yedges = results[
        'exit_energy_angle_energy_edges_mev'
    ]

    fig, ax = plt.subplots()

    mesh = ax.pcolormesh(
        xedges,
        yedges,
        np.log10(
            corr + 1.0
        ),
        shading='auto',
        cmap=RED_ORANGE_YELLOW_CMAP,
    )

    ax.set_xlabel(
        'Exit cos(angle) relative to +x'
    )

    ax.set_ylabel(
        'Exit energy (MeV)'
    )

    ax.set_title(
        'Exit Energy vs Exit Angle (Binned)'
    )

    fig.colorbar(
        mesh,
        ax=ax,
        label='log10(count + 1)',
    )

    show_fig(
        fig,
        'exit_energy_vs_exit_angle_binned',
    )

    fig, ax = plt.subplots()

    cum_dose = np.cumsum(
        results['dose_vs_depth_gy']
        * results['depth_bin_widths_cm']
    )

    ax.plot(
        depth_cm,
        cum_dose,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
    )

    ax.set_xlabel(
        'Depth x (cm)'
    )

    ax.set_ylabel(
        'Cumulative dose metric'
    )

    ax.set_title(
        'Cumulative Deposited Dose vs Depth'
    )

    show_fig(
        fig,
        'cumulative_deposited_dose_vs_depth',
    )

    print(
        f'Saved plots to: '
        f'{os.path.abspath(cfg.plot_dir)}'
    )


_RAYLEIGH_TALLY_PATCH_APPLIED = False

if not globals().get(
    '_RAYLEIGH_TALLY_PATCH_APPLIED',
    False,
):
    _RAYLEIGH_TALLY_PATCH_APPLIED = True

    try:
        if (
            'Rayl'
            not in PRIMARY_GAMMA_PROCESSES
        ):
            PRIMARY_GAMMA_PROCESSES = tuple(
                list(PRIMARY_GAMMA_PROCESSES)
                + ['Rayl']
            )
    except NameError:
        PRIMARY_GAMMA_PROCESSES = (
            'compt',
            'phot',
            'conv',
            'Rayl',
        )

    _orig_EventState_reset = (
        EventState.reset
    )

    def _rayleigh_EventState_reset(
        self,
    ):
        _orig_EventState_reset(self)

        self.rayleigh_count = 0
        self.rayleigh_event_counts_by_depth = None

    EventState.reset = (
        _rayleigh_EventState_reset
    )

    _orig_SimulationTallies_init = (
        SimulationTallies.__init__
    )

    def _rayleigh_SimulationTallies_init(
        self,
        cfg,
    ):
        _orig_SimulationTallies_init(
            self,
            cfg,
        )

        if not hasattr(
            self,
            'rayleigh_vs_depth_counts',
        ):
            self.rayleigh_vs_depth_counts = np.zeros(
                cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            self,
            'rayleigh_vs_depth_counts_sumsq',
        ):
            self.rayleigh_vs_depth_counts_sumsq = np.zeros(
                cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            self,
            'rayleigh_hist',
        ):
            self.rayleigh_hist = np.zeros(
                cfg.max_scatter_bin + 2,
                dtype=np.int64,
            )

        if hasattr(
            self,
            'process_counts_all_tracks',
        ):
            self.process_counts_all_tracks.setdefault(
                'Rayl',
                0,
            )

        if hasattr(
            self,
            'process_counts_primary_total',
        ):
            self.process_counts_primary_total.setdefault(
                'Rayl',
                0,
            )

        if hasattr(
            self,
            'process_counts_primary_eventwise',
        ):
            self.process_counts_primary_eventwise.setdefault(
                'Rayl',
                0,
            )

        if hasattr(
            self,
            'buildup_primary_counts',
        ):
            self.buildup_primary_counts.setdefault(
                'rayleigh_only',
                0,
            )

            self.buildup_primary_counts.setdefault(
                'mixed_rayleigh_compton',
                0,
            )

    SimulationTallies.__init__ = (
        _rayleigh_SimulationTallies_init
    )

    _orig_MyEventAction_BeginOfEventAction = (
        MyEventAction.BeginOfEventAction
    )

    def _rayleigh_BeginOfEventAction(
        self,
        event,
    ):
        _orig_MyEventAction_BeginOfEventAction(
            self,
            event,
        )

        self.state.rayleigh_count = 0

        self.state.rayleigh_event_counts_by_depth = np.zeros(
            self.cfg.depth_nbins,
            dtype=float,
        )

    MyEventAction.BeginOfEventAction = (
        _rayleigh_BeginOfEventAction
    )

    _orig_MyEventAction_EndOfEventAction = (
        MyEventAction.EndOfEventAction
    )

    def _rayleigh_EndOfEventAction(
        self,
        event,
    ):
        before_depth_sum = float(
            np.sum(
                getattr(
                    self.tallies,
                    'rayleigh_vs_depth_counts',
                    np.array([0.0]),
                )
            )
        )

        before_hist_sum = int(
            np.sum(
                getattr(
                    self.tallies,
                    'rayleigh_hist',
                    np.array([0]),
                )
            )
        )

        _orig_MyEventAction_EndOfEventAction(
            self,
            event,
        )

        if not hasattr(
            self.tallies,
            'rayleigh_vs_depth_counts',
        ):
            self.tallies.rayleigh_vs_depth_counts = np.zeros(
                self.cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            self.tallies,
            'rayleigh_vs_depth_counts_sumsq',
        ):
            self.tallies.rayleigh_vs_depth_counts_sumsq = np.zeros(
                self.cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            self.tallies,
            'rayleigh_hist',
        ):
            self.tallies.rayleigh_hist = np.zeros(
                self.cfg.max_scatter_bin + 2,
                dtype=np.int64,
            )

        after_depth_sum = float(
            np.sum(
                self.tallies.rayleigh_vs_depth_counts
            )
        )

        after_hist_sum = int(
            np.sum(
                self.tallies.rayleigh_hist
            )
        )

        rayleigh_event_counts = getattr(
            self.state,
            'rayleigh_event_counts_by_depth',
            None,
        )

        if (
            rayleigh_event_counts is not None
            and after_depth_sum == before_depth_sum
        ):
            self.tallies.rayleigh_vs_depth_counts += (
                rayleigh_event_counts
            )

            self.tallies.rayleigh_vs_depth_counts_sumsq += (
                rayleigh_event_counts ** 2
            )

        if (
            getattr(
                self.state,
                'terminated',
                None,
            )
            == 'transmitted'
            and after_hist_sum == before_hist_sum
        ):
            ridx = int(
                getattr(
                    self.state,
                    'rayleigh_count',
                    0,
                )
            )

            if ridx > self.cfg.max_scatter_bin:
                ridx = (
                    self.cfg.max_scatter_bin
                    + 1
                )

            self.tallies.rayleigh_hist[
                ridx
            ] += 1

    MyEventAction.EndOfEventAction = (
        _rayleigh_EndOfEventAction
    )

    _orig_MySteppingAction_UserSteppingAction = (
        MySteppingAction.UserSteppingAction
    )

    def _rayleigh_UserSteppingAction(
        self,
        step,
    ):
        before_rayleigh_count = int(
            getattr(
                self.state,
                'rayleigh_count',
                0,
            )
        )

        _orig_MySteppingAction_UserSteppingAction(
            self,
            step,
        )

        track = step.GetTrack()
        post = step.GetPostStepPoint()

        proc = post.GetProcessDefinedStep()

        proc_name = (
            proc.GetProcessName()
            if proc is not None
            else None
        )

        if proc_name != 'Rayl':
            return

        if getattr(
            self.state,
            'primary_track_id',
            None,
        ) is None:
            return

        if (
            track.GetTrackID()
            != self.state.primary_track_id
        ):
            return

        if (
            track.GetDefinition().GetParticleName()
            != 'gamma'
        ):
            return

        after_rayleigh_count = int(
            getattr(
                self.state,
                'rayleigh_count',
                0,
            )
        )

        if (
            after_rayleigh_count
            > before_rayleigh_count
        ):
            return

        post_pos = post.GetPosition()

        x_int_cm = max(
            0.0,
            min(
                self.cfg.wall_thickness_cm,
                float(post_pos.x / cm),
            ),
        )

        if getattr(
            self.state,
            'rayleigh_event_counts_by_depth',
            None,
        ) is None:
            self.state.rayleigh_event_counts_by_depth = np.zeros(
                self.cfg.depth_nbins,
                dtype=float,
            )

        self.state.rayleigh_count = (
            before_rayleigh_count + 1
        )

        self.state.rayleigh_event_counts_by_depth[
            self.tallies.depth_bin_index(
                x_int_cm
            )
        ] += 1.0

        self.state.collided_any = True
        self.state.had_interaction = True

        self.state.primary_processes_seen.add(
            'Rayl'
        )

        if hasattr(
            self.tallies,
            'process_counts_primary_total',
        ):
            self.tallies.process_counts_primary_total.setdefault(
                'Rayl',
                0,
            )

        if (
            self.state.first_interaction_x_cm
            is None
        ):
            self.state.first_interaction_x_cm = (
                x_int_cm
            )

        self.state.last_interaction_x_cm = (
            x_int_cm
        )

    MySteppingAction.UserSteppingAction = (
        _rayleigh_UserSteppingAction
    )

    _orig_build_results = build_results

    def _rayleigh_build_results(
        cfg,
        tallies,
        runtime_s,
        obj_info,
    ):
        results = _orig_build_results(
            cfg,
            tallies,
            runtime_s,
            obj_info,
        )

        N = int(
            results['summary']['num_photons']
        )

        if N <= 0:
            N = 1

        if not hasattr(
            tallies,
            'rayleigh_vs_depth_counts',
        ):
            tallies.rayleigh_vs_depth_counts = np.zeros(
                cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            tallies,
            'rayleigh_vs_depth_counts_sumsq',
        ):
            tallies.rayleigh_vs_depth_counts_sumsq = np.zeros(
                cfg.depth_nbins,
                dtype=np.float64,
            )

        if not hasattr(
            tallies,
            'rayleigh_hist',
        ):
            tallies.rayleigh_hist = np.zeros(
                cfg.max_scatter_bin + 2,
                dtype=np.int64,
            )

        rayleigh_sem_counts = safe_sem_from_sum_sumsq(
            tallies.rayleigh_vs_depth_counts,
            tallies.rayleigh_vs_depth_counts_sumsq,
            N,
        )

        rayleigh_per_primary_per_cm = (
            tallies.rayleigh_vs_depth_counts
            / (
                N
                * tallies.depth_bin_widths_cm
            )
        )

        rayleigh_per_primary_per_cm_sigma = (
            rayleigh_sem_counts
            / tallies.depth_bin_widths_cm
        )

        rayleigh_total = int(
            tallies.process_counts_primary_total.get(
                'Rayl',
                0,
            )
        )

        rayleigh_event_total = int(
            tallies.process_counts_primary_eventwise.get(
                'Rayl',
                0,
            )
        )

        results['summary'][
            'primary_rayleigh_total_count'
        ] = rayleigh_total

        results['summary'][
            'primary_rayleigh_event_count'
        ] = rayleigh_event_total

        results['summary'][
            'primary_rayleigh_events_per_primary'
        ] = float(
            rayleigh_total / N
        )

        results['summary'][
            'primary_rayleigh_event_fraction'
        ] = float(
            rayleigh_event_total / N
        )

        results[
            'rayleigh_vs_depth_counts'
        ] = (
            tallies.rayleigh_vs_depth_counts.copy()
        )

        results[
            'rayleigh_vs_depth_sem'
        ] = (
            rayleigh_sem_counts.copy()
        )

        results[
            'rayleigh_per_primary_per_cm'
        ] = (
            rayleigh_per_primary_per_cm.copy()
        )

        results[
            'rayleigh_per_primary_per_cm_sigma'
        ] = (
            rayleigh_per_primary_per_cm_sigma.copy()
        )

        results[
            'rayleigh_hist_counts'
        ] = (
            tallies.rayleigh_hist.copy()
        )

        results[
            'rayleigh_hist_sigma'
        ] = np.sqrt(
            np.maximum(
                tallies.rayleigh_hist,
                0,
            )
        )

        results[
            'rayleigh_hist_bins'
        ] = np.arange(
            cfg.max_scatter_bin + 2,
            dtype=int,
        )

        return results

    build_results = _rayleigh_build_results

    _orig_print_summary = print_summary

    def _rayleigh_print_summary(
        results,
    ):
        _orig_print_summary(
            results
        )

        summary = results.get(
            'summary',
            {},
        )

        print('')

        print(
            f"primary_rayleigh_total_count           = "
            f"{int(summary.get('primary_rayleigh_total_count', 0))}"
        )

        print(
            f"primary_rayleigh_event_count           = "
            f"{int(summary.get('primary_rayleigh_event_count', 0))}"
        )

        print(
            f"primary_rayleigh_events_per_primary    = "
            f"{float(summary.get('primary_rayleigh_events_per_primary', 0.0)):.6e}"
        )

        print(
            f"primary_rayleigh_event_fraction        = "
            f"{float(summary.get('primary_rayleigh_event_fraction', 0.0)):.6e}"
        )

        print('')

    print_summary = _rayleigh_print_summary

    _orig_plot_results = plot_results

    def _rayleigh_plot_results(
        results,
        cfg,
        sweep_results=None,
    ):
        _orig_plot_results(
            results,
            cfg,
            sweep_results=sweep_results,
        )

        os.makedirs(
            cfg.plot_dir,
            exist_ok=True,
        )

        if (
            'rayleigh_hist_counts'
            in results
        ):
            rayleigh_hist = np.asarray(
                results[
                    'rayleigh_hist_counts'
                ],
                dtype=float,
            )

            rayleigh_sigma = np.asarray(
                results.get(
                    'rayleigh_hist_sigma',
                    np.sqrt(
                        np.maximum(
                            rayleigh_hist,
                            0,
                        )
                    ),
                ),
                dtype=float,
            )

            max_scatter_bin = (
                len(rayleigh_hist) - 2
            )

            rayleigh_x = np.arange(
                len(rayleigh_hist)
            )

            fig, ax = plt.subplots()

            ax.bar(
                rayleigh_x,
                rayleigh_hist,
                color=PLOT_RED,
                edgecolor=PLOT_RED,
                linewidth=PLOT_EDGE_WIDTH,
            )

            ax.errorbar(
                rayleigh_x,
                rayleigh_hist,
                yerr=rayleigh_sigma,
                fmt='none',
                capsize=2,
                ecolor=PLOT_RED,
                elinewidth=PLOT_ERROR_LINE_WIDTH,
                capthick=PLOT_ERROR_LINE_WIDTH,
            )

            ax.set_xlabel(
                f'# Rayleigh scatters before transmission '
                f'(last bin >= {max_scatter_bin})'
            )

            ax.set_ylabel(
                'Counts'
            )

            ax.set_title(
                'Rayleigh-Order Distribution '
                'with Statistical Error Bars'
            )

            style_black_red_figure(
                fig
            )

            fig.tight_layout()

            fig.savefig(
                os.path.join(
                    cfg.plot_dir,
                    'rayleigh_order_distribution_with_error_bars.png',
                ),
                dpi=cfg.plot_dpi,
                bbox_inches='tight',
                facecolor=PLOT_BLACK,
                edgecolor=PLOT_BLACK,
            )

            plt.show()

        if (
            'rayleigh_per_primary_per_cm'
            in results
        ):
            depth_cm = np.asarray(
                results['depth_cm'],
                dtype=float,
            )

            rayleigh_y = np.asarray(
                results[
                    'rayleigh_per_primary_per_cm'
                ],
                dtype=float,
            )

            rayleigh_sigma = np.asarray(
                results.get(
                    'rayleigh_per_primary_per_cm_sigma',
                    np.zeros_like(
                        rayleigh_y
                    ),
                ),
                dtype=float,
            )

            fig, ax = plt.subplots()

            ax.plot(
                depth_cm,
                rayleigh_y,
                color=PLOT_RED,
                linewidth=PLOT_LINE_WIDTH,
            )

            ax.fill_between(
                depth_cm,
                rayleigh_y - rayleigh_sigma,
                rayleigh_y + rayleigh_sigma,
                color=PLOT_RED,
                alpha=PLOT_UNCERTAINTY_ALPHA,
                linewidth=0.0,
            )

            ax.set_xlabel(
                'Depth x (cm)'
            )

            ax.set_ylabel(
                'Rayleigh interactions '
                'per primary per cm'
            )

            ax.set_title(
                'Primary-Gamma Rayleigh '
                'Scattering vs Depth'
            )

            style_black_red_figure(
                fig
            )

            fig.tight_layout()

            fig.savefig(
                os.path.join(
                    cfg.plot_dir,
                    'primary_gamma_rayleigh_vs_depth.png',
                ),
                dpi=cfg.plot_dpi,
                bbox_inches='tight',
                facecolor=PLOT_BLACK,
                edgecolor=PLOT_BLACK,
            )

            plt.show()

    plot_results = _rayleigh_plot_results



# -----------------------------------------------------------------------------
# Expert validation corrections
# -----------------------------------------------------------------------------
# These corrections make all depth-resolved material quantities explicitly
# Shield-volume-only, calculate dose and process-count uncertainties on the
# independent primary-history basis, record reproducibility metadata, retire
# the redundant wR=1 equivalent-dose plot, and replace the non-standard Gy*cm
# cumulative metric with a cumulative deposited-energy fraction.
_EXPERT_VALIDATION_PATCH_APPLIED = True
CFG.run_thickness_sweep = False


def _expert_text(value: Any) -> str:
    """Return a stable Python string for G4String/bytes wrapper values."""
    if value is None:
        return ''
    if isinstance(value, bytes):
        try:
            return value.decode('utf-8', errors='replace').strip()
        except Exception:
            return str(value).strip()
    text = str(value).strip()
    if len(text) >= 3 and text[0] == 'b' and text[1] in ("'", '"') and text[-1] == text[1]:
        text = text[2:-1]
    return text.strip()


def _expert_step_volume(step_point):
    try:
        volume = step_point.GetPhysicalVolume()
    except Exception:
        volume = None

    if volume is None:
        try:
            touchable = step_point.GetTouchableHandle()
            volume = touchable.GetVolume() if touchable is not None else None
        except Exception:
            volume = None
    return volume


def _expert_volume_descriptors(step_point) -> Tuple[str, str, str]:
    volume = _expert_step_volume(step_point)
    if volume is None:
        return ('', '', '')

    physical_name = ''
    logical_name = ''
    material_name = ''

    try:
        physical_name = _expert_text(volume.GetName())
    except Exception:
        pass

    logical = None
    try:
        logical = volume.GetLogicalVolume()
    except Exception:
        logical = None

    if logical is not None:
        try:
            logical_name = _expert_text(logical.GetName())
        except Exception:
            pass
        try:
            material = logical.GetMaterial()
            material_name = _expert_text(material.GetName()) if material is not None else ''
        except Exception:
            pass

    return physical_name, logical_name, material_name


def _expert_volume_name(step_point) -> Optional[str]:
    physical_name, logical_name, _material_name = _expert_volume_descriptors(step_point)
    return physical_name or logical_name or None


def _expert_point_inside_shield_geometry(step_point, cfg, tolerance_cm: float = 1.0e-7) -> bool:
    """Coordinate fallback used only when pybind does not expose a usable volume name."""
    try:
        x_cm, y_cm, z_cm = vec_cm(step_point.GetPosition())
    except Exception:
        return False

    half_y = 0.5 * float(cfg.wall_size_y_cm)
    half_z = 0.5 * float(cfg.wall_size_z_cm)
    thickness = float(cfg.wall_thickness_cm)
    tol = max(float(tolerance_cm), 1.0e-9 * max(1.0, thickness))

    return (
        -tol <= x_cm <= thickness + tol
        and abs(y_cm) <= half_y + tol
        and abs(z_cm) <= half_z + tol
    )


def _expert_is_shield(step_point, cfg) -> bool:
    """Identify Shield steps robustly across geant4_pybind name wrappers.

    Earlier code compared ``str(volume.GetName())`` with one exact literal.
    Some bindings expose a bytes-like G4String or a decorated physical-volume
    name, which made every depth/process tally zero and produced empty plots.
    This implementation checks physical/logical/material descriptors and uses
    the known geometry only when no trustworthy descriptor is available.
    """
    physical_name, logical_name, material_name = _expert_volume_descriptors(step_point)
    names = [physical_name.lower(), logical_name.lower()]

    if any(name == 'shield' or 'shield' in name for name in names if name):
        return True

    material_lower = material_name.lower()
    if material_lower in ('g4_concrete', 'concrete') or 'concrete' in material_lower:
        return _expert_point_inside_shield_geometry(step_point, cfg)

    # The Shield is the only placed volume occupying this exact rectangular
    # coordinate region.  Coordinate membership is therefore a safe final
    # fallback even when a boundary-step wrapper reports World or an opaque
    # volume name.  This prevents pybind naming differences from silently
    # zeroing the depth, process and track tallies.
    return _expert_point_inside_shield_geometry(step_point, cfg)


def _expert_scalar_sem(total: float, total_sumsq: float, n: int) -> float:
    if n <= 1:
        return 0.0
    mean = float(total) / float(n)
    ex2 = float(total_sumsq) / float(n)
    variance = max(ex2 - mean * mean, 0.0)
    return math.sqrt(variance / float(n))


def wilson_interval(count: int, total: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if total <= 0:
        return (0.0, 1.0)
    n = float(total)
    p = min(max(float(count) / n, 0.0), 1.0)
    z2 = z * z
    denominator = 1.0 + z2 / n
    centre = (p + z2 / (2.0 * n)) / denominator
    half = z * math.sqrt(max(p * (1.0 - p) / n + z2 / (4.0 * n * n), 0.0)) / denominator
    return (max(0.0, centre - half), min(1.0, centre + half))


def _expert_geant4_version() -> str:
    for name in ('G4Version', 'G4VERSION_NUMBER'):
        value = globals().get(name)
        if value is None:
            continue
        try:
            return str(value() if callable(value) else value)
        except Exception:
            continue
    return 'not available from geant4_pybind'


_expert_previous_event_reset = EventState.reset

def _expert_event_reset(self) -> None:
    _expert_previous_event_reset(self)
    self.process_counts_primary_event = {p: 0 for p in PRIMARY_GAMMA_PROCESSES}
    self.process_counts_all_tracks_event = {p: 0 for p in PRIMARY_GAMMA_PROCESSES}

EventState.reset = _expert_event_reset


_expert_previous_tallies_init = SimulationTallies.__init__

def _expert_tallies_init(self, cfg) -> None:
    _expert_previous_tallies_init(self, cfg)
    self.process_counts_primary_sumsq = {p: 0.0 for p in PRIMARY_GAMMA_PROCESSES}
    self.process_counts_all_tracks_sumsq = {p: 0.0 for p in PRIMARY_GAMMA_PROCESSES}
    self.total_edep_mev_event_sum = 0.0
    self.total_edep_mev_event_sumsq = 0.0

SimulationTallies.__init__ = _expert_tallies_init


_expert_previous_begin_event = MyEventAction.BeginOfEventAction

def _expert_begin_event(self, event) -> None:
    _expert_previous_begin_event(self, event)
    self.state.process_counts_primary_event = {p: 0 for p in PRIMARY_GAMMA_PROCESSES}
    self.state.process_counts_all_tracks_event = {p: 0 for p in PRIMARY_GAMMA_PROCESSES}

MyEventAction.BeginOfEventAction = _expert_begin_event


_expert_previous_end_event = MyEventAction.EndOfEventAction

def _expert_end_event(self, event) -> None:
    primary_counts = dict(getattr(self.state, 'process_counts_primary_event', {}))
    all_counts = dict(getattr(self.state, 'process_counts_all_tracks_event', {}))
    event_edep = float(np.sum(getattr(self.state, 'edep_event_mev', np.zeros(1))))

    _expert_previous_end_event(self, event)

    for process_name in PRIMARY_GAMMA_PROCESSES:
        pcount = float(primary_counts.get(process_name, 0))
        acount = float(all_counts.get(process_name, 0))
        self.tallies.process_counts_primary_sumsq[process_name] = (
            self.tallies.process_counts_primary_sumsq.get(process_name, 0.0) + pcount * pcount
        )
        self.tallies.process_counts_all_tracks_sumsq[process_name] = (
            self.tallies.process_counts_all_tracks_sumsq.get(process_name, 0.0) + acount * acount
        )

    self.tallies.total_edep_mev_event_sum += event_edep
    self.tallies.total_edep_mev_event_sumsq += event_edep * event_edep

MyEventAction.EndOfEventAction = _expert_end_event


def _expert_primary_track_point(self, point: Tuple[float, float, float]) -> None:
    if not self.state.track_points_cm or point != self.state.track_points_cm[-1]:
        self.state.track_points_cm.append(point)


def _expert_user_stepping_action(self, step) -> None:
    track = step.GetTrack()
    pre = step.GetPreStepPoint()
    post = step.GetPostStepPoint()

    x0_cm, y0_cm, z0_cm = vec_cm(pre.GetPosition())
    x1_cm, y1_cm, z1_cm = vec_cm(post.GetPosition())
    pre_is_shield = _expert_is_shield(pre, self.cfg)
    post_is_shield = _expert_is_shield(post, self.cfg)

    proc = post.GetProcessDefinedStep()
    proc_name = proc.GetProcessName() if proc is not None else None
    pname = track.GetDefinition().GetParticleName()

    # All-track photon-process tallies are restricted to interactions whose
    # step originates in the concrete Shield volume.
    if pre_is_shield and proc_name in self.tallies.process_counts_all_tracks:
        self.tallies.process_counts_all_tracks[proc_name] += 1
        self.state.process_counts_all_tracks_event[proc_name] = (
            self.state.process_counts_all_tracks_event.get(proc_name, 0) + 1
        )

    # Local energy deposition is assigned only when the depositing step is in
    # the Shield. This prevents air/world deposition from being clamped into
    # the first or last material bin.
    edep_mev = float(step.GetTotalEnergyDeposit() / MeV)
    if pre_is_shield and edep_mev > 0.0:
        x_mid_cm = max(0.0, min(self.cfg.wall_thickness_cm, 0.5 * (x0_cm + x1_cm)))
        self.state.edep_event_mev[self.tallies.depth_bin_index(x_mid_cm)] += edep_mev

    # Secondary forward exits are boundary crossings from Shield to world.
    if (
        pre_is_shield
        and not post_is_shield
        and x1_cm >= self.cfg.wall_thickness_cm - 1.0e-9
        and track.GetParentID() > 0
        and pname in self.tallies.secondary_exit_hist_counts
    ):
        ekin_mev = float(post.GetKineticEnergy() / MeV)
        self.tallies._increment_1d_hist(
            ekin_mev,
            self.tallies.secondary_energy_hist_edges_mev,
            self.tallies.secondary_exit_hist_counts[pname],
        )
        self.tallies.secondary_exit_counts_forward[pname] += 1

    if self.state.primary_track_id is None:
        return
    if track.GetTrackID() != self.state.primary_track_id or pname != 'gamma':
        return

    # Visual tracks now contain only the segment inside the Shield.
    if self.state.record_this_track and (pre_is_shield or post_is_shield):
        if pre_is_shield:
            _expert_primary_track_point(
                self,
                (max(0.0, min(self.cfg.wall_thickness_cm, x0_cm)), y0_cm, z0_cm),
            )
        _expert_primary_track_point(
            self,
            (max(0.0, min(self.cfg.wall_thickness_cm, x1_cm)), y1_cm, z1_cm),
        )

    if pre_is_shield:
        self.tallies.tally_track_length_1d_into(
            self.state.fluence_event_cm,
            x0_cm=x0_cm,
            x1_cm=x1_cm,
            step_length_cm=float(step.GetStepLength() / cm),
        )

    if pre_is_shield and proc_name is not None and proc_name not in self.INTERACTION_IGNORE:
        self.state.collided_any = True

        if proc_name in self.tallies.process_counts_primary_total:
            self.tallies.process_counts_primary_total[proc_name] += 1
            self.state.process_counts_primary_event[proc_name] = (
                self.state.process_counts_primary_event.get(proc_name, 0) + 1
            )
            self.state.primary_processes_seen.add(proc_name)

        x_int_cm = max(0.0, min(self.cfg.wall_thickness_cm, x1_cm))
        if not self.state.had_interaction:
            self.state.had_interaction = True
            self.state.first_interaction_x_cm = x_int_cm
        self.state.last_interaction_x_cm = x_int_cm

        if proc_name == 'compt':
            self.state.compton_count += 1
            preE_mev = float(pre.GetKineticEnergy() / MeV)
            postE_mev = float(post.GetKineticEnergy() / MeV)
            dE_mev = max(0.0, preE_mev - postE_mev)
            self.state.primary_deltaE_event_mev[
                self.tallies.depth_bin_index(x_int_cm)
            ] += dE_mev

        if proc_name == 'Rayl':
            self.state.rayleigh_count = int(getattr(self.state, 'rayleigh_count', 0)) + 1
            if getattr(self.state, 'rayleigh_event_counts_by_depth', None) is None:
                self.state.rayleigh_event_counts_by_depth = np.zeros(
                    self.cfg.depth_nbins,
                    dtype=float,
                )
            self.state.rayleigh_event_counts_by_depth[
                self.tallies.depth_bin_index(x_int_cm)
            ] += 1.0

    # Classify a primary only when it actually leaves the Shield volume.
    if self.state.terminated is None and pre_is_shield and not post_is_shield:
        if x1_cm >= self.cfg.wall_thickness_cm - 1.0e-9:
            self.state.terminated = 'transmitted'
            self.state.exit_energy_mev = float(post.GetKineticEnergy() / MeV)
            self.state.exit_cosine = vec_x_unit(track.GetMomentumDirection())
        elif x1_cm <= 1.0e-9:
            self.state.terminated = 'backescaped'
        else:
            self.state.terminated = 'other_escape'

MySteppingAction.UserSteppingAction = _expert_user_stepping_action


_expert_previous_run_simulation = run_simulation

def run_simulation(cfg) -> Dict[str, Any]:
    seed = getattr(cfg, 'random_seed', None)
    if seed is None:
        seed = int(time.time_ns() % 2147483647) or 1
        cfg.random_seed = seed
    seed = int(seed)
    np.random.seed(seed % (2 ** 32 - 1))
    try:
        G4Random.setTheSeed(seed)
    except Exception:
        try:
            CLHEP.HepRandom.setTheSeed(seed)
        except Exception:
            pass
    return _expert_previous_run_simulation(cfg)


_expert_previous_build_results = build_results

def build_results(cfg, tallies, runtime_s, obj_info):
    results = _expert_previous_build_results(cfg, tallies, runtime_s, obj_info)
    summary = results.setdefault('summary', {})
    N = max(1, int(cfg.num_photons))

    density_g_cm3 = float(tallies.material_density_g_cm3)
    widths_cm = np.asarray(tallies.depth_bin_widths_cm, dtype=float)
    bin_volume_cm3 = float(cfg.wall_size_y_cm) * float(cfg.wall_size_z_cm) * widths_cm
    bin_mass_kg = density_g_cm3 * bin_volume_cm3 / 1000.0

    # Whole-slice average absorbed dose per independent primary history.
    dose_per_primary_gy = (
        np.asarray(tallies.edep_vs_depth_mev, dtype=float) / N
        * MEV_TO_J
        / bin_mass_kg
    )
    edep_sem_mev = safe_sem_from_sum_sumsq(
        tallies.edep_vs_depth_mev,
        tallies.edep_vs_depth_mev_sumsq,
        N,
    )
    dose_per_primary_sem_gy = edep_sem_mev * MEV_TO_J / bin_mass_kg

    wall_mass_kg = (
        density_g_cm3
        * float(cfg.wall_size_y_cm)
        * float(cfg.wall_size_z_cm)
        * float(cfg.wall_thickness_cm)
        / 1000.0
    )
    total_edep_mev = float(tallies.total_edep_mev_event_sum)
    total_edep_per_primary_mev = total_edep_mev / N
    total_edep_per_primary_sem_mev = _expert_scalar_sem(
        tallies.total_edep_mev_event_sum,
        tallies.total_edep_mev_event_sumsq,
        N,
    )
    whole_wall_dose_per_primary_gy = total_edep_per_primary_mev * MEV_TO_J / wall_mass_kg
    whole_wall_dose_per_primary_sem_gy = total_edep_per_primary_sem_mev * MEV_TO_J / wall_mass_kg
    simulated_run_whole_wall_dose_gy = whole_wall_dose_per_primary_gy * N
    simulated_run_whole_wall_dose_sigma_gy = whole_wall_dose_per_primary_sem_gy * N

    results['dose_vs_depth_gy'] = dose_per_primary_gy.copy()
    results['dose_vs_depth_gy_sigma'] = dose_per_primary_sem_gy.copy()
    results['whole_slice_average_dose_per_primary_gy'] = dose_per_primary_gy.copy()
    results['whole_slice_average_dose_per_primary_gy_sigma'] = dose_per_primary_sem_gy.copy()

    # Backward-compatible summary names retain their former run-total meaning,
    # but now carry a correctly scaled uncertainty.
    summary['total_absorbed_dose_gy'] = float(simulated_run_whole_wall_dose_gy)
    summary['total_absorbed_dose_gy_sigma'] = float(simulated_run_whole_wall_dose_sigma_gy)
    summary['whole_wall_average_dose_per_primary_gy'] = float(whole_wall_dose_per_primary_gy)
    summary['whole_wall_average_dose_per_primary_gy_sigma'] = float(whole_wall_dose_per_primary_sem_gy)
    summary['simulated_run_whole_wall_average_dose_gy'] = float(simulated_run_whole_wall_dose_gy)
    summary['simulated_run_whole_wall_average_dose_gy_sigma'] = float(simulated_run_whole_wall_dose_sigma_gy)
    summary['dose_scoring_definition'] = (
        'energy deposited in the Shield volume divided by the mass of the complete '
        'wall cross-sectional depth slice; plotted values are per primary history'
    )

    # Retain data aliases for compatibility, but the redundant wR=1 plot is
    # intentionally not displayed by the corrected plot wrapper below.
    results['equivalent_dose_proxy_vs_depth_sv_wR1'] = dose_per_primary_gy.copy()
    results['equivalent_dose_proxy_vs_depth_sv_wR1_sigma'] = dose_per_primary_sem_gy.copy()
    summary['equivalent_dose_proxy_retired'] = False
    summary['equivalent_dose_proxy_definition'] = (
        'deterministic wR=1 copy of absorbed dose; retained only as a transparent derived display'
    )

    deposited = np.asarray(tallies.edep_vs_depth_mev, dtype=float)
    cumulative = np.cumsum(np.maximum(deposited, 0.0))
    if cumulative.size and cumulative[-1] > 0.0:
        cumulative = cumulative / cumulative[-1]
    results['cumulative_deposited_energy_fraction'] = cumulative

    primary_mean = {}
    primary_sem = {}
    all_mean = {}
    all_sem = {}
    for process_name in PRIMARY_GAMMA_PROCESSES:
        ptotal = float(tallies.process_counts_primary_total.get(process_name, 0))
        atotal = float(tallies.process_counts_all_tracks.get(process_name, 0))
        primary_mean[process_name] = ptotal / N
        all_mean[process_name] = atotal / N
        primary_sem[process_name] = _expert_scalar_sem(
            ptotal,
            tallies.process_counts_primary_sumsq.get(process_name, 0.0),
            N,
        )
        all_sem[process_name] = _expert_scalar_sem(
            atotal,
            tallies.process_counts_all_tracks_sumsq.get(process_name, 0.0),
            N,
        )
    results['process_counts_primary_mean_per_primary'] = primary_mean
    results['process_counts_primary_sem_per_primary'] = primary_sem
    results['process_counts_all_tracks_mean_per_primary'] = all_mean
    results['process_counts_all_tracks_sem_per_primary'] = all_sem

    partition_counts = {
        'forward_transmit': int(tallies.transmitted_count),
        'backscatter': int(tallies.backscatter_count),
        'removed': int(tallies.removed_count),
        'other_escape': int(tallies.other_escape_count),
    }
    results['partition_fraction_wilson95'] = {
        key: wilson_interval(count, N)
        for key, count in partition_counts.items()
    }
    transmitted_n = int(tallies.transmitted_count)
    results['buildup_primary_fraction_wilson95'] = {
        key: wilson_interval(int(count), transmitted_n)
        for key, count in tallies.buildup_primary_counts.items()
    }

    summary.update({
        'geant4_version': _expert_geant4_version(),
        'physics_list': 'G4EmStandardPhysics_option4',
        'production_cut_mm': float(cfg.production_cut_mm),
        'random_seed': int(getattr(cfg, 'random_seed', 0)),
        'wall_size_y_cm': float(cfg.wall_size_y_cm),
        'wall_size_z_cm': float(cfg.wall_size_z_cm),
        'depth_nbins': int(cfg.depth_nbins),
        'depth_bin_width_cm': float(np.mean(widths_cm)) if widths_cm.size else math.nan,
        'energy_hist_bins': int(cfg.energy_hist_bins),
        'angle_hist_bins': int(cfg.angle_hist_bins),
        'correlation_energy_bins': int(cfg.corr_energy_bins),
        'correlation_angle_bins': int(cfg.corr_angle_bins),
        'beam_half_angle_deg': float(cfg.beam_half_angle_deg),
        'source_to_wall_gap_cm': float(cfg.source_to_wall_gap_cm),
        'scoring_volume': 'Shield only',
        'track_display_scope': 'Shield segment only',
        'process_count_uncertainty_basis': 'per-primary-history mean and SEM',
        'global_energy_balance': (
            'not scored: the current run does not tally the energy of every '
            'particle escaping every world boundary, so no closure claim is made'
        ),
    })

    summary['simulation_thickness_feature_build'] = SIMULATION_THICKNESS_FEATURE_BUILD
    return results


_expert_previous_plot_results = plot_results

def plot_results(results, cfg, sweep_results=None):
    outer_show = plt.show

    def corrected_show(*args, **kwargs):
        figure = plt.gcf()
        axes = [axis for axis in figure.axes if getattr(axis, 'get_title', lambda: '')()]
        axis = axes[0] if axes else (figure.axes[0] if figure.axes else None)
        title = axis.get_title().strip() if axis is not None else ''
        title_lower = title.lower()

        if 'cumulative deposited dose' in title_lower:
            plt.close(figure)
            return None

        if axis is not None and title == 'Dose vs Depth':
            axis.set_title('Whole-Slice Average Absorbed Dose per Primary vs Depth')
            axis.set_ylabel('Absorbed dose per primary (Gy/primary)')

        if axis is not None and 'Equivalent Dose Proxy' in title:
            axis.set_title('Absorbed Dose with Unit Weighting wR=1 (Derived)')
            axis.set_ylabel('Weighted absorbed dose per primary (wR=1; Gy/primary)')

        if axis is not None and title.startswith('3D Photon Tracks'):
            axis.set_title('3D Photon Tracks Within the Shield')

        if axis is not None and title == 'Fluence vs depth (N=%d)' % int(cfg.num_photons):
            axis.set_title('Shield Track-Length Density vs Depth')
            axis.set_ylabel('Photon track length per primary per cm of depth')

        if axis is not None and title == 'Actual Deposited Energy vs Depth':
            axis.set_title('Shield-Only Deposited Energy vs Depth')

        if axis is not None and title == 'Energy deposition proxy vs depth (N=%d)' % int(cfg.num_photons):
            axis.set_title('Shield Primary-Compton Energy-Loss Proxy vs Depth')

        # Sparse non-zero secondary spectra remain visible; the GUI/report
        # assessment explicitly limits their interpretation.

        def replace_bar_values(values, errors=None, lower=None, upper=None, ylabel=None):
            if axis is None:
                return
            values_arr = np.asarray(values, dtype=float)
            bars = list(axis.patches)[: values_arr.size]
            for bar, value in zip(bars, values_arr):
                bar.set_height(float(value))
            for line in list(axis.lines):
                try:
                    line.remove()
                except Exception:
                    pass
            for collection in list(axis.collections):
                try:
                    collection.remove()
                except Exception:
                    pass
            x = np.arange(values_arr.size, dtype=float)
            if lower is not None and upper is not None:
                lo = np.maximum(values_arr - np.asarray(lower, dtype=float), 0.0)
                hi = np.maximum(np.asarray(upper, dtype=float) - values_arr, 0.0)
                yerr = np.vstack([lo, hi])
            else:
                yerr = np.asarray(errors, dtype=float) if errors is not None else None
            if yerr is not None:
                axis.errorbar(
                    x,
                    values_arr,
                    yerr=yerr,
                    fmt='none',
                    capsize=3,
                    ecolor=PLOT_RED,
                    elinewidth=PLOT_ERROR_LINE_WIDTH,
                    capthick=PLOT_ERROR_LINE_WIDTH,
                )
            if ylabel:
                axis.set_ylabel(ylabel)
            axis.relim()
            axis.autoscale_view()
            axis.set_ylim(bottom=0.0)

        if title == 'Primary-Gamma Process Counts':
            names = list(results['process_counts_primary_total'].keys())
            replace_bar_values(
                [results['process_counts_primary_mean_per_primary'][key] for key in names],
                errors=[results['process_counts_primary_sem_per_primary'][key] for key in names],
                ylabel='Mean Shield interactions per primary ±1σ SEM',
            )
            axis.set_title('Primary-Gamma Shield Process Rate per Primary')

        if title == 'All-Track Photon Process Counts':
            names = list(results['process_counts_all_tracks_total'].keys())
            replace_bar_values(
                [results['process_counts_all_tracks_mean_per_primary'][key] for key in names],
                errors=[results['process_counts_all_tracks_sem_per_primary'][key] for key in names],
                ylabel='Mean Shield interactions per primary ±1σ SEM',
            )
            axis.set_title('All-Track Shield Photon Process Rate per Primary')

        if title == 'Primary-Gamma Eventwise Process Incidence':
            names = list(results['process_counts_primary_eventwise'].keys())
            counts = [int(results['process_counts_primary_eventwise'][key]) for key in names]
            values = [count / max(1, int(cfg.num_photons)) for count in counts]
            intervals = [wilson_interval(count, int(cfg.num_photons)) for count in counts]
            replace_bar_values(
                values,
                lower=[item[0] for item in intervals],
                upper=[item[1] for item in intervals],
                ylabel='Fraction of primary histories (95% Wilson CI)',
            )

        if title == 'Forward / Backscatter / Removal Partition':
            labels = ('forward_transmit', 'backscatter', 'removed', 'other_escape')
            counts = (
                int(results['summary']['transmitted_count']),
                int(results['summary']['backscatter_count']),
                int(results['summary']['removed_count']),
                int(results['summary']['other_escape_count']),
            )
            values = [count / max(1, int(cfg.num_photons)) for count in counts]
            intervals = [results['partition_fraction_wilson95'][key] for key in labels]
            replace_bar_values(
                values,
                lower=[item[0] for item in intervals],
                upper=[item[1] for item in intervals],
                ylabel='Fraction of primaries (95% Wilson CI)',
            )

        if title == 'Primary-Gamma Buildup Decomposition':
            labels = list(results['buildup_primary_counts'].keys())
            values = [results['buildup_primary_fraction_of_transmitted'][key] for key in labels]
            intervals = [results['buildup_primary_fraction_wilson95'][key] for key in labels]
            replace_bar_values(
                values,
                lower=[item[0] for item in intervals],
                upper=[item[1] for item in intervals],
                ylabel='Fraction of transmitted primaries (95% Wilson CI)',
            )

        style_black_red_figure(figure)
        try:
            figure.tight_layout()
        except Exception:
            pass
        return outer_show(*args, **kwargs)

    plt.show = corrected_show
    try:
        _expert_previous_plot_results(results, cfg, sweep_results=sweep_results)
    finally:
        plt.show = outer_show

    # In single-thickness GUI mode, retain an attenuation diagnostic without
    # launching hidden sweep simulations.  This is one measured point at the
    # manually selected thickness, not a thickness-dependence curve.
    if sweep_results is None:
        summary = results.get('summary', {})
        thickness = float(summary.get('wall_thickness_cm', cfg.wall_thickness_cm))
        mu_eff = float(summary.get('effective_mu_per_cm', math.nan))
        mu_sigma = float(summary.get('effective_mu_per_cm_sigma', math.nan))
        if math.isfinite(thickness) and math.isfinite(mu_eff):
            fig, ax = plt.subplots()
            yerr = mu_sigma if math.isfinite(mu_sigma) and mu_sigma >= 0.0 else None
            ax.errorbar(
                [thickness],
                [mu_eff],
                yerr=None if yerr is None else [yerr],
                fmt='o',
                capsize=4,
                color=PLOT_RED,
                ecolor=PLOT_RED,
                markerfacecolor=PLOT_RED,
                markeredgecolor=PLOT_RED,
                linewidth=PLOT_LINE_WIDTH,
                elinewidth=PLOT_ERROR_LINE_WIDTH,
            )
            span = max(1.0, 0.08 * max(thickness, 1.0))
            ax.set_xlim(max(0.0, thickness - span), thickness + span)
            ax.set_ylim(bottom=0.0)
            ax.set_xlabel('Selected wall thickness (cm)')
            ax.set_ylabel('Effective attenuation coefficient (1/cm)')
            ax.set_title('Effective Attenuation Coefficient at Selected Thickness')
            style_black_red_figure(fig)
            fig.tight_layout()
            fig.savefig(
                os.path.join(cfg.plot_dir, 'effective_mu_at_selected_thickness.png'),
                dpi=cfg.plot_dpi,
                bbox_inches='tight',
                facecolor=PLOT_BLACK,
                edgecolor=PLOT_BLACK,
            )
            plt.show()

    depth_cm = np.asarray(results['depth_cm'], dtype=float)
    cumulative = np.asarray(results.get('cumulative_deposited_energy_fraction', []), dtype=float)
    if depth_cm.size and cumulative.size == depth_cm.size:
        fig, ax = plt.subplots()
        ax.plot(depth_cm, cumulative, color=PLOT_RED, linewidth=PLOT_LINE_WIDTH)
        ax.set_xlabel('Depth x (cm)')
        ax.set_ylabel('Cumulative fraction of Shield-deposited energy')
        ax.set_ylim(0.0, 1.02)
        ax.set_title('Cumulative Deposited-Energy Fraction vs Depth')
        style_black_red_figure(fig)
        fig.tight_layout()
        fig.savefig(
            os.path.join(cfg.plot_dir, 'cumulative_deposited_energy_fraction_vs_depth.png'),
            dpi=cfg.plot_dpi,
            bbox_inches='tight',
            facecolor=PLOT_BLACK,
            edgecolor=PLOT_BLACK,
        )
        plt.show()


# Replace the legacy summary printer so the console and GUI output use the
# corrected definitions and do not present the retired wR=1 proxy as an
# independent result.
def print_summary(results: Dict[str, Any]) -> None:
    s = results['summary']

    def _console_optional_number(value: Any, fmt: str = '.3f', unavailable: str = 'NOT ESTABLISHED') -> str:
        try:
            x = float(value)
        except (TypeError, ValueError):
            return unavailable
        if not math.isfinite(x):
            return unavailable
        return format(x, fmt)

    print('\n==================== RESULTS ====================')
    # v7.3i.2c.2: the required concrete wall thickness is the primary result.
    # Only the supervisory audit may define it; the plot-bearing batch and
    # legacy Wilson/transmission fields remain diagnostic.
    _audit = None
    for _key in ('v73i2_photon_validation_audit', 'v73h_photon_validation_audit', 'v73f_photon_validation_audit'):
        _candidate = s.get(_key)
        if isinstance(_candidate, dict):
            _audit = _candidate
            break
    if isinstance(_audit, dict):
        _min_t = _console_optional_number(_audit.get('minimum_verified_compliant_thickness_cm'))
        _verified_t = _console_optional_number(_audit.get('verified_passing_thickness_cm'))
        _final_t = _console_optional_number(_audit.get('final_tested_thickness_cm'))
        _state = str(_audit.get('final_state', 'STATISTICALLY_UNRESOLVED')).replace('_', ' ')
        print('')
        print('******** PRIMARY SHIELDING THICKNESS RESULT ********')
        if _min_t != 'NOT ESTABLISHED':
            print(f'REQUIRED CONCRETE WALL THICKNESS = {_min_t} cm  [VERIFIED MINIMUM COMPLIANT]')
        elif _verified_t != 'NOT ESTABLISHED':
            print('REQUIRED CONCRETE WALL THICKNESS = NOT ESTABLISHED')
            print(f'Verified passing thickness      = {_verified_t} cm; thinner candidate(s) remain unresolved')
        else:
            print('REQUIRED CONCRETE WALL THICKNESS = NOT ESTABLISHED')
        print(f'PRIMARY RESULT STATUS           = {_state}')
        if _final_t != 'NOT ESTABLISHED':
            print(f'Final supervisory tested thickness = {_final_t} cm')
        print('*****************************************************')
        print('')
    for key in (
        'beam_half_angle_deg',
        'source_to_wall_gap_cm',
        'beam_mev_input',
        'source_energy_mev_internal',
        'source_mode',
        'spectrum_case',
        'spectrum_case_origin',
        'source_energy_hist_max_mev',
        'wall_thickness_cm',
        'num_photons',
        'material_name_used',
        'material_density_g_cm3',
        'geant4_version',
        'physics_list',
        'transport_backend',
        'available_logical_cpus',
        'transport_threads',
        'thread_policy',
        'native_worker_build',
        'native_adapter_build',
        'production_cut_mm',
        'random_seed',
        'scoring_volume',
    ):
        if key in s and s[key] is not None:
            print(f'{key:40s} = {s[key]}')

    print('')
    for key in (
        'transmitted_count',
        'uncollided_transmitted_count',
        'backscatter_count',
        'removed_count',
        'other_escape_count',
    ):
        print(f'{key:40s} = {int(s.get(key, 0))}')

    print('')
    for key in (
        'transmitted_fraction',
        'uncollided_transmitted_fraction',
        'backscatter_fraction',
        'removed_fraction',
        'other_escape_fraction',
    ):
        sigma_key = key + '_sigma'
        print(f'{key:40s} = {float(s.get(key, 0.0)):.6e} ± {float(s.get(sigma_key, 0.0)):.2e}')

    _bf = s.get('buildup_factor_estimate')
    _bfs = s.get('buildup_factor_sigma')
    try:
        _bf_ok = math.isfinite(float(_bf)) and math.isfinite(float(_bfs))
    except (TypeError, ValueError):
        _bf_ok = False
    if _bf_ok:
        print(f"{'buildup_factor_estimate':40s} = {float(_bf):.6e} ± {float(_bfs):.2e}")
    else:
        print(f"{'buildup_factor_estimate':40s} = NOT AVAILABLE — no finite uncollided-transmission reference")
    _mu=s.get('effective_mu_per_cm'); _mus=s.get('effective_mu_per_cm_sigma')
    try: _mu_ok=math.isfinite(float(_mu)) and math.isfinite(float(_mus))
    except (TypeError,ValueError): _mu_ok=False
    if _mu_ok:
        print(f"{'effective_mu_per_cm':40s} = {float(_mu):.6e} ± {float(_mus):.2e}")
    else:
        print(f"{'effective_mu_per_cm':40s} = NOT ESTABLISHED — zero/sparse transmission censors the log-transmission diagnostic")

    print('')
    print(
        f"{'whole_wall_average_dose_per_primary_gy':40s} = "
        f"{float(s.get('whole_wall_average_dose_per_primary_gy', 0.0)):.6e} ± "
        f"{float(s.get('whole_wall_average_dose_per_primary_gy_sigma', 0.0)):.2e}"
    )
    print(
        f"{'simulated_run_whole_wall_average_dose_gy':40s} = "
        f"{float(s.get('simulated_run_whole_wall_average_dose_gy', 0.0)):.6e} ± "
        f"{float(s.get('simulated_run_whole_wall_average_dose_gy_sigma', 0.0)):.2e}"
    )
    print(f"{'dose_scoring_definition':40s} = {s.get('dose_scoring_definition', 'not available')}")
    print(f"{'global_energy_balance':40s} = {s.get('global_energy_balance', 'not scored')}")

    requirement = s.get('simulation_concrete_requirement', {})
    if isinstance(requirement, dict) and requirement:
        print('')
        print('-------- GEANT4 CONCRETE THICKNESS SWEEP REQUIREMENT --------')
        print(f"{'target_transmission':40s} = {float(requirement.get('target_transmission', math.nan)):.6e}")
        print(f"{'sweep_range_cm':40s} = {float(requirement.get('sweep_min_thickness_cm', math.nan)):.3f} to {float(requirement.get('sweep_max_thickness_cm', math.nan)):.3f}")
        print(f"{'sweep_points':40s} = {int(requirement.get('sweep_point_count', 0))}")
        print(f"{'sweep_photons_per_point':40s} = {int(requirement.get('sweep_num_photons_per_point', 0))}")
        print(f"{'zero_count_wilson95_upper':40s} = {float(requirement.get('zero_count_wilson95_upper_at_sweep_n', math.nan)):.6e}")
        print(f"{'required_concrete_95_cm':40s} = {_console_optional_number(requirement.get('required_thickness_95_cm'))}")
        print(f"{'requirement_relation':40s} = {requirement.get('required_thickness_relation', 'not available')}")
        print(f"{'first_compliant_grid_cm':40s} = {_console_optional_number(requirement.get('first_compliant_simulated_thickness_cm'))}")
        print(f"{'legacy_primary_history_fields_role':40s} = DIAGNOSTIC ONLY — photon compliance uses the cumulative all-photon dry-air weighted supervisory response")
        _k12d_summary = results.get('summary', {}) if isinstance(results, dict) else {}
        print(f"{'final_supervisory_tested_thickness_cm':40s} = {float(requirement.get('final_supervisory_tested_thickness_cm', requirement.get('final_validation_thickness_cm', math.nan))):.3f}")
        print(f"{'cumulative_plotting_primary_histories':40s} = {int(_k12d_summary.get('v73k12d_plot_aggregate_histories', requirement.get('final_validation_primary_photons', 0)) or 0)}")
        print(f"{'cumulative_plotting_supervisory_batches':40s} = {int(_k12d_summary.get('v73k12d_plot_aggregate_batch_count', 0) or 0)}")
        print(f"{'trajectory_visualization_source_histories':40s} = {int(_k12d_summary.get('v73k12d_trajectory_sample_histories', 0) or 0)}")
        print(f"{'extra_plot_only_geant4_histories':40s} = {int(_k12d_summary.get('v73k12d_plot_extra_native_histories', 0) or 0)}")
        print(f"{'quantitative_plot_dataset_role':40s} = {_k12d_summary.get('v73k12e_plot_dataset_role', 'CUMULATIVE SUPERVISORY REUSE — same histories as decision estimator')}")
        print(f"{'trajectory_visualization_role':40s} = BOUNDED FIRST FULL-FEATURE SUPERVISORY BATCH")
        print(f"{'simulation_thickness_status':40s} = {requirement.get('status', 'not available')}")
        print('TVL may seed the editable search range only; TVL is not used for Geant4 PASS/FAIL, interpolation, extrapolation, or final thickness selection.')
        print('-------------------------------------------------------------')

    print('')
    print('process_counts_primary_total           =', results.get('process_counts_primary_total', {}))
    print('process_counts_primary_mean_per_primary=', results.get('process_counts_primary_mean_per_primary', {}))
    print('process_counts_primary_sem_per_primary =', results.get('process_counts_primary_sem_per_primary', {}))
    print('process_counts_primary_eventwise       =', results.get('process_counts_primary_eventwise', {}))
    print('process_counts_all_tracks_total        =', results.get('process_counts_all_tracks_total', {}))
    print('process_counts_all_tracks_mean_per_primary=', results.get('process_counts_all_tracks_mean_per_primary', {}))
    print('process_counts_all_tracks_sem_per_primary =', results.get('process_counts_all_tracks_sem_per_primary', {}))
    print('buildup_primary_counts                 =', results.get('buildup_primary_counts', {}))
    print('secondary_exit_counts_forward          =', results.get('secondary_exit_counts_forward', {}))
    print('')
    print(f"{'tracks_saved':40s} = {int(s.get('tracks_saved', 0))}")
    print(f"{'tracks_rejected_short':40s} = {int(s.get('tracks_rejected_short', 0))}")
    print(f"{'runtime_s':40s} = {float(s.get('runtime_s', 0.0)):.3f}")
    print('=================================================\n')



# ---------------------------------------------------------------------------
# Native C++ photon backend
# ---------------------------------------------------------------------------
# The original Python/geant4_pybind implementation remains in this file for
# validation and an explicit developer fallback.  Normal GUI runs use the
# native worker so per-step and per-event scoring no longer crosses the
# Python/C++ binding boundary.
_PYTHON_PHOTON_RUN_SIMULATION = run_simulation
NATIVE_PHOTON_ADAPTER_BUILD = '2026-08-05-native-photon-adapter-v2-safe-half-cpu-mt'



def _available_logical_cpu_count() -> int:
    """Return logical CPUs available to this process, respecting Linux affinity."""
    try:
        affinity = os.sched_getaffinity(0)
        if affinity:
            return max(1, len(affinity))
    except (AttributeError, NotImplementedError, OSError):
        pass
    return max(1, int(os.cpu_count() or 1))


def _resolve_transport_threads(cfg) -> Tuple[int, int]:
    """Use half the available logical CPUs unless an explicit safe override is set."""
    available = _available_logical_cpu_count()
    configured = getattr(cfg, 'transport_threads', None)

    for name in ('PHOTON_WALL_TRANSPORT_THREADS', 'GEANT4_TRANSPORT_THREADS'):
        value = os.environ.get(name, '').strip()
        if value:
            try:
                configured = int(value)
            except ValueError as exc:
                raise ValueError(f'{name} must be a positive integer') from exc
            break

    if configured is None or int(configured) <= 0:
        threads = max(1, available // 2)
    else:
        threads = int(configured)

    threads = max(1, min(threads, available))
    cfg.transport_threads = threads
    return available, threads

def _native_photon_runner_path() -> str:
    from pathlib import Path

    configured = os.environ.get('PHOTON_WALL_NATIVE_RUNNER', '').strip()
    if configured:
        candidate = Path(configured).expanduser().resolve()
    else:
        candidate = (
            Path(__file__).resolve().parent
            / 'native_photon_worker'
            / 'run_native_photon_worker.sh'
        )
    if not candidate.is_file():
        raise FileNotFoundError(
            'Native photon-worker runner was not found: '
            f'{candidate}. Run install_photon_wall_app.sh.'
        )
    return str(candidate)


def _native_photon_input(cfg) -> Dict[str, Any]:
    seed = getattr(cfg, 'random_seed', None)
    if seed is None:
        seed = int(time.time_ns() % 2147483647) or 1
        cfg.random_seed = seed
    seed = int(seed)
    available_logical_cpus, transport_threads = _resolve_transport_threads(cfg)

    source_mode = get_source_mode(cfg)
    spectrum_energies: List[float] = []
    spectrum_probabilities: List[float] = []
    importance_enabled = bool(getattr(cfg, '_v73k_inline_importance_enabled', False))
    importance_original_probabilities: List[float] = []
    importance_group_index_by_bin: List[int] = []
    importance_statistical_batches = int(getattr(cfg, '_v73k_inline_statistical_batches', 1) or 1)
    if source_mode == 'digitized_spectrum':
        # v7.3k.7 performs the grouped source-importance proposal inside ONE
        # native Geant4 process per pilot/production replicate.  The full
        # physical spectrum p(E) and proposal q(E) are both passed to C++; C++
        # samples q(E) event-by-event and accumulates exact p/q-weighted source
        # and exit sufficient statistics.  Ordinary coarse/fine/analog runs
        # continue to sample the unmodified physical spectrum directly.
        if importance_enabled:
            custom_e = getattr(cfg, '_v73k_inline_spectrum_energies_mev', None)
            custom_p = getattr(cfg, '_v73k_inline_original_probabilities', None)
            custom_q = getattr(cfg, '_v73k_inline_proposal_probabilities', None)
            custom_g = getattr(cfg, '_v73k_inline_group_index_by_bin', None)
            e = np.asarray(custom_e, dtype=float).ravel()
            p_phys = np.asarray(custom_p, dtype=float).ravel()
            q_prop = np.asarray(custom_q, dtype=float).ravel()
            gidx = np.asarray(custom_g, dtype=int).ravel()
            if e.size == 0 or not (e.size == p_phys.size == q_prop.size == gidx.size):
                raise RuntimeError('v7.3k.7 inline importance arrays are empty or length-inconsistent')
            if np.any(~np.isfinite(e)) or np.any(~np.isfinite(p_phys)) or np.any(~np.isfinite(q_prop)):
                raise RuntimeError('v7.3k.7 inline importance arrays contain nonfinite values')
            if np.any(e <= 0.0) or np.any(p_phys < 0.0) or np.any(q_prop < 0.0):
                raise RuntimeError('v7.3k.7 inline importance arrays contain invalid values')
            if not float(np.sum(p_phys)) > 0.0 or not float(np.sum(q_prop)) > 0.0:
                raise RuntimeError('v7.3k.7 inline p(E)/q(E) has zero total probability mass')
            p_support = p_phys > 0.0
            q_support = q_prop > 0.0
            if not np.array_equal(p_support, q_support):
                raise RuntimeError('v7.3k.7 inline proposal support must exactly match physical source support')
            if np.any(gidx[p_support] < 0):
                raise RuntimeError('v7.3k.7 every positive-support source bin requires an importance group')
            if importance_statistical_batches <= 0 or importance_statistical_batches > int(cfg.num_photons):
                raise RuntimeError('v7.3k.7 invalid inline statistical batch count')
            spectrum_energies = e.tolist()
            spectrum_probabilities = q_prop.tolist()
            importance_original_probabilities = p_phys.tolist()
            importance_group_index_by_bin = gidx.tolist()
        else:
            custom_e = getattr(cfg, '_v73k_spectrum_energies_mev', None)
            custom_p = getattr(cfg, '_v73k_spectrum_probabilities', None)
            if custom_e is not None or custom_p is not None:
                spectrum_energies = np.asarray(custom_e, dtype=float).ravel().tolist()
                spectrum_probabilities = np.asarray(custom_p, dtype=float).ravel().tolist()
                if not spectrum_energies or len(spectrum_energies) != len(spectrum_probabilities):
                    raise RuntimeError('v7.3k custom native source spectrum arrays are empty or inconsistent')
                if not np.all(np.isfinite(spectrum_energies)) or not np.all(np.isfinite(spectrum_probabilities)):
                    raise RuntimeError('v7.3k custom native source spectrum contains nonfinite values')
                if np.any(np.asarray(spectrum_energies) <= 0.0) or np.any(np.asarray(spectrum_probabilities) < 0.0):
                    raise RuntimeError('v7.3k custom native source spectrum contains invalid values')
                if not float(np.sum(spectrum_probabilities)) > 0.0:
                    raise RuntimeError('v7.3k custom native source spectrum has zero probability mass')
            else:
                spectrum = get_digitized_spectrum(cfg)
                spectrum_energies = np.asarray(
                    spectrum['energy_mev_center'], dtype=float
                ).tolist()
                spectrum_probabilities = np.asarray(
                    spectrum['probability_mass_bin'], dtype=float
                ).tolist()

    source_hist_max = float(source_energy_hist_max_mev(cfg))
    return {
        'source_mode': source_mode,
        'spectrum_case': (
            get_spectrum_case(cfg)
            if source_mode == 'digitized_spectrum'
            else None
        ),
        'spectrum_energies_mev': spectrum_energies,
        'spectrum_probabilities': spectrum_probabilities,
        'importance_sampling_enabled': importance_enabled,
        'importance_original_probabilities': importance_original_probabilities,
        'importance_group_index_by_bin': importance_group_index_by_bin,
        'importance_statistical_batches': importance_statistical_batches if importance_enabled else 1,
        'mono_energy_mev': float(derived_source_energy_mev(cfg)),
        'beam_half_angle_deg': float(cfg.beam_half_angle_deg),
        'wall_thickness_cm': float(cfg.wall_thickness_cm),
        'wall_size_y_cm': float(cfg.wall_size_y_cm),
        'wall_size_z_cm': float(cfg.wall_size_z_cm),
        'source_to_wall_gap_cm': float(cfg.source_to_wall_gap_cm),
        'production_cut_mm': float(cfg.production_cut_mm),
        'concrete_density_g_cm3': (
            None
            if cfg.concrete_density_g_cm3 is None
            else float(cfg.concrete_density_g_cm3)
        ),
        'num_photons': int(cfg.num_photons),
        'random_seed': seed,
        'available_logical_cpus': int(available_logical_cpus),
        'transport_threads': int(transport_threads),
        'thread_policy': 'half of available logical CPUs',
        'max_raw_exit_samples': int(cfg.max_raw_exit_samples),
        'n_visual_tracks': int(cfg.n_visual_tracks),
        'min_pts_after_decimate': int(cfg.min_pts_after_decimate),
        'obj_decimate': int(cfg.obj_decimate),
        'depth_nbins': int(cfg.depth_nbins),
        'max_scatter_bin': int(cfg.max_scatter_bin),
        'energy_hist_bins': int(cfg.energy_hist_bins),
        'angle_hist_bins': int(cfg.angle_hist_bins),
        'corr_energy_bins': int(cfg.corr_energy_bins),
        'corr_angle_bins': int(cfg.corr_angle_bins),
        'secondary_energy_hist_bins': int(cfg.secondary_energy_hist_bins),
        'geant4_control_verbose': int(cfg.geant4_control_verbose),
        'geant4_run_verbose': int(cfg.geant4_run_verbose),
        'geant4_event_verbose': int(cfg.geant4_event_verbose),
        'geant4_tracking_verbose': int(cfg.geant4_tracking_verbose),
        'depth_edges_cm': np.linspace(
            0.0,
            float(cfg.wall_thickness_cm),
            int(cfg.depth_nbins) + 1,
            dtype=float,
        ).tolist(),
        'energy_hist_edges_mev': np.linspace(
            0.0,
            source_hist_max,
            int(cfg.energy_hist_bins) + 1,
            dtype=float,
        ).tolist(),
        'exit_cosine_hist_edges': np.linspace(
            -1.0,
            1.0,
            int(cfg.angle_hist_bins) + 1,
            dtype=float,
        ).tolist(),
        'corr_cos_edges': np.linspace(
            -1.0,
            1.0,
            int(cfg.corr_angle_bins) + 1,
            dtype=float,
        ).tolist(),
        'corr_energy_edges_mev': np.linspace(
            0.0,
            source_hist_max,
            int(cfg.corr_energy_bins) + 1,
            dtype=float,
        ).tolist(),
        'secondary_energy_hist_edges_mev': np.linspace(
            0.0,
            source_hist_max,
            int(cfg.secondary_energy_hist_bins) + 1,
            dtype=float,
        ).tolist(),
    }



def _v7_all_photon_air_kerma_metrics(raw: Dict[str, Any], cfg: Any) -> Dict[str, Any]:
    """Build the all-photon downstream dry-air energy-absorption response and eventwise SEM."""
    n = int(getattr(cfg, 'num_photons', 0) or 0)
    sx = float(raw.get('source_air_kerma_weight_sum', 0.0))
    sx2 = float(raw.get('source_air_kerma_weight_sumsq', 0.0))
    sy = float(raw.get('exit_air_kerma_weight_sum', 0.0))
    sy2 = float(raw.get('exit_air_kerma_weight_sumsq', 0.0))
    sxy = float(raw.get('source_exit_air_kerma_weight_cross_sum', 0.0))
    crossings = int(raw.get('all_photon_forward_crossing_count', 0) or 0)
    max_weight = float(raw.get('max_air_kerma_weight_per_photon', float('nan')))

    response = float('nan')
    sem = float('nan')
    source_mean = float('nan')
    if n > 0 and sx > 0.0 and math.isfinite(sx):
        source_mean = sx / float(n)
        response = sy / sx
        if n > 1:
            mx = sx / float(n)
            my = sy / float(n)
            var_x = max(0.0, (sx2 - float(n) * mx * mx) / float(n - 1))
            var_y = max(0.0, (sy2 - float(n) * my * my) / float(n - 1))
            cov_xy = (sxy - float(n) * mx * my) / float(n - 1)
            variance_mean_ratio = (
                var_y + response * response * var_x - 2.0 * response * cov_xy
            ) / (float(n) * mx * mx)
            sem = math.sqrt(max(0.0, variance_mean_ratio))

    target = float('nan')
    for name in (
        'required_max_photon_transmission', 'photon_transmission_target',
        'target_transmission', 'required_transmission', 'max_transmission',
    ):
        try:
            value = float(getattr(cfg, name))
        except Exception:
            continue
        if math.isfinite(value) and 0.0 < value <= 1.0:
            target = value
            break

    verification_upper = float('nan')
    statistically_verified = False
    zero_probability_upper = float('nan')
    zero_weighted_upper = float('nan')
    if n > 0 and math.isfinite(response) and response >= 0.0:
        if crossings == 0:
            z = 1.959963984540054
            z2 = z * z
            zero_probability_upper = z2 / (float(n) + z2)
            ratio = float('nan')
            if math.isfinite(source_mean) and source_mean > 0.0 and math.isfinite(max_weight) and max_weight > 0.0:
                ratio = max(1.0, max_weight / source_mean)
            if math.isfinite(ratio):
                zero_weighted_upper = zero_probability_upper * ratio
            # A zero weighted score never certifies the shielding target because
            # no conditional exit-photon weight distribution was sampled.
            verification_upper = float('nan')
        elif crossings >= 20 and math.isfinite(sem) and sem >= 0.0:
            verification_upper = response + 1.6448536269514722 * sem
        # Sparse non-zero scores are intentionally not declared resolved by a
        # Gaussian approximation. They remain real Monte Carlo observations,
        # but the shielding target is not statistically verified.
        if math.isfinite(target) and math.isfinite(verification_upper):
            statistically_verified = bool(response <= target and verification_upper <= target)

    return {
        'all_photon_air_kerma_response': response,
        'all_photon_air_kerma_response_sem': sem,
        'all_photon_air_kerma_statistical_verification_upper': verification_upper,
        'all_photon_air_kerma_statistically_verified': statistically_verified,
        'all_photon_forward_crossing_count': crossings,
        'source_air_kerma_weight_mean': source_mean,
        'max_air_kerma_weight_per_photon': max_weight,
        'zero_crossing_wilson95_probability_upper': zero_probability_upper,
        'zero_crossing_weighted_response_upper': zero_weighted_upper,
        'all_photon_air_kerma_response_definition': (
            'sum E*(mu_en/rho)_dry-air for every unique forward photon track crossing '
            'the downstream shield boundary, including secondaries, divided by the '
            'incident-primary source-weight sum'
        ),
        'air_kerma_response_reference': str(raw.get('air_kerma_response_reference', 'NIST dry air mu_en/rho')),
        'photonuclear_physics_status': (
            'ENABLED in native worker: EM Option 4 + G4EmExtraPhysics gamma-nuclear + '
            'FTFP_BERT_HP hadron/neutron transport'
        ),
    }


def _native_photon_fill_tallies(cfg, raw: Dict[str, Any]) -> SimulationTallies:
    tallies = SimulationTallies(cfg)

    float_arrays = (
        'fluence_tracklen_cm',
        'fluence_tracklen_cm_sumsq',
        'edep_vs_depth_mev',
        'edep_vs_depth_mev_sumsq',
        'primary_gamma_deltaE_vs_depth_mev',
        'primary_gamma_deltaE_vs_depth_mev_sumsq',
        'rayleigh_vs_depth_counts',
        'rayleigh_vs_depth_counts_sumsq',
    )
    int_arrays = (
        'first_interaction_hist',
        'last_interaction_hist',
        'scatter_hist',
        'rayleigh_hist',
        'energy_hist_counts',
        'exit_cosine_hist_counts',
    )
    for name in float_arrays:
        setattr(tallies, name, np.asarray(raw[name], dtype=np.float64))
    for name in int_arrays:
        setattr(tallies, name, np.asarray(raw[name], dtype=np.int64))

    tallies.exit_energy_angle_hist = np.asarray(
        raw['exit_energy_angle_hist_flat'], dtype=np.int64
    ).reshape(int(cfg.corr_angle_bins), int(cfg.corr_energy_bins))

    tallies.secondary_exit_hist_counts = {
        name: np.asarray(values, dtype=np.int64)
        for name, values in raw['secondary_exit_hist_counts'].items()
    }
    tallies.secondary_exit_counts_forward = {
        name: int(value)
        for name, value in raw['secondary_exit_counts_forward'].items()
    }
    tallies.process_counts_all_tracks = {
        name: int(value)
        for name, value in raw['process_counts_all_tracks'].items()
    }
    tallies.process_counts_primary_total = {
        name: int(value)
        for name, value in raw['process_counts_primary_total'].items()
    }
    tallies.process_counts_primary_eventwise = {
        name: int(value)
        for name, value in raw['process_counts_primary_eventwise'].items()
    }
    tallies.process_counts_all_tracks_sumsq = {
        name: float(value)
        for name, value in raw['process_counts_all_tracks_sumsq'].items()
    }
    tallies.process_counts_primary_sumsq = {
        name: float(value)
        for name, value in raw['process_counts_primary_sumsq'].items()
    }
    tallies.buildup_primary_counts = {
        name: int(value)
        for name, value in raw['buildup_primary_counts'].items()
    }

    for name in (
        'transmitted_count',
        'uncollided_transmitted_count',
        'backscatter_count',
        'removed_count',
        'other_escape_count',
        'tracks_saved',
        'tracks_rejected_short',
    ):
        setattr(tallies, name, int(raw[name]))

    tallies.total_edep_mev_event_sum = float(raw['total_edep_mev_event_sum'])
    tallies.total_edep_mev_event_sumsq = float(raw['total_edep_mev_event_sumsq'])
    tallies.all_photon_forward_crossing_count = int(raw.get('all_photon_forward_crossing_count', 0))
    tallies.source_air_kerma_weight_sum = float(raw.get('source_air_kerma_weight_sum', 0.0))
    tallies.source_air_kerma_weight_sumsq = float(raw.get('source_air_kerma_weight_sumsq', 0.0))
    tallies.exit_air_kerma_weight_sum = float(raw.get('exit_air_kerma_weight_sum', 0.0))
    tallies.exit_air_kerma_weight_sumsq = float(raw.get('exit_air_kerma_weight_sumsq', 0.0))
    tallies.source_exit_air_kerma_weight_cross_sum = float(raw.get('source_exit_air_kerma_weight_cross_sum', 0.0))
    tallies.material_name_used = str(raw['material_name_used'])
    tallies.material_density_g_cm3 = float(raw['material_density_g_cm3'])
    tallies.transmitted_energy_sample_mev = [
        float(value) for value in raw['transmitted_energy_sample_mev']
    ]
    tallies.exit_cosine_sample = [
        float(value) for value in raw['exit_cosine_sample']
    ]
    tallies.exit_sample_seen = int(tallies.transmitted_count)
    tallies.tracks_cm = [
        (
            [float(value) for value in track[0]],
            [float(value) for value in track[1]],
            [float(value) for value in track[2]],
        )
        for track in raw['tracks_cm']
    ]

    partition_total = (
        tallies.transmitted_count
        + tallies.backscatter_count
        + tallies.removed_count
        + tallies.other_escape_count
    )
    if partition_total != int(cfg.num_photons):
        raise RuntimeError(
            'Native photon-worker particle partition does not equal the '
            f'requested histories: {partition_total} != {cfg.num_photons}'
        )
    return tallies


def _native_photon_print_obj_info(cfg, obj_info) -> None:
    if obj_info is None:
        return
    print(f'[OBJ EXPORT] {cfg.obj_filepath}')
    for key in (
        'mtl_filepath',
        'tracks_written',
        'segment_prisms_written',
        'vertices',
        'faces',
        'wall_included',
        'track_radius_cm',
    ):
        if key in obj_info:
            print(f'  {key:24s} = {obj_info[key]}')
    print(f'  {"scale":24s} = {cfg.obj_scale}')
    print(f'  {"decimate":24s} = {cfg.obj_decimate}')
    print('')


def run_simulation(cfg) -> Dict[str, Any]:
    if os.environ.get('PHOTON_WALL_USE_PYTHON_BACKEND', '').strip() == '1':
        print(
            'Photon backend: explicit legacy Python/geant4_pybind fallback',
            flush=True,
        )
        return _PYTHON_PHOTON_RUN_SIMULATION(cfg)

    import json
    import subprocess
    import tempfile
    from pathlib import Path

    runner = _native_photon_runner_path()
    worker_input = _native_photon_input(cfg)
    print(f'Photon core build: {NATIVE_PHOTON_ADAPTER_BUILD}', flush=True)
    print(f'Native photon worker runner: {runner}', flush=True)
    print(
        'Photon transport backend: native C++ Geant4; full Geant4 report retained',
        flush=True,
    )
    print(
        'Photon CPU policy: '
        f"{worker_input['transport_threads']} Geant4 worker threads from "
        f"{worker_input['available_logical_cpus']} available logical CPUs "
        '(safe half-CPU default)',
        flush=True,
    )

    with tempfile.TemporaryDirectory(prefix='photon-wall-native-') as temporary:
        temporary_path = Path(temporary)
        input_path = temporary_path / 'photon_input.json'
        output_path = temporary_path / 'photon_raw_results.json'
        input_path.write_text(
            json.dumps(worker_input, indent=2),
            encoding='utf-8',
        )

        process = subprocess.Popen(
            [runner, str(input_path), str(output_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=os.environ.copy(),
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end='', flush=True)
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(
                f'Native photon worker exited with status {return_code}'
            )
        if not output_path.is_file():
            raise RuntimeError(
                'Native photon worker completed without creating its result file'
            )
        raw = json.loads(output_path.read_text(encoding='utf-8'))

    cfg.random_seed = int(raw['random_seed'])
    tallies = _native_photon_fill_tallies(cfg, raw)

    obj_info = None
    if (
        cfg.export_obj
        and (
            cfg.include_wall_in_obj
            or len(tallies.tracks_cm) > 0
        )
    ):
        obj_info = export_scene_to_obj(
            tracks_cm=tallies.tracks_cm,
            filepath=cfg.obj_filepath,
            scale=cfg.obj_scale,
            decimate=cfg.obj_decimate,
            wall_thickness_cm=cfg.wall_thickness_cm,
            wall_size_y_cm=cfg.wall_size_y_cm,
            wall_size_z_cm=cfg.wall_size_z_cm,
            include_wall=cfg.include_wall_in_obj,
            track_radius_cm=cfg.obj_track_radius_cm,
        )

    results = build_results(
        cfg,
        tallies,
        float(raw['runtime_s']),
        obj_info,
    )
    summary = results.setdefault('summary', {})
    summary['geant4_version'] = str(raw['geant4_version'])
    summary['physics_list'] = str(raw['physics_list'])
    summary['transport_backend'] = str(raw['transport_backend'])
    summary['native_worker_build'] = str(raw['native_worker_build'])
    summary['native_adapter_build'] = NATIVE_PHOTON_ADAPTER_BUILD
    summary['available_logical_cpus'] = int(raw.get('available_logical_cpus', worker_input['available_logical_cpus']))
    summary['transport_threads'] = int(raw.get('transport_threads', worker_input['transport_threads']))
    summary['thread_policy'] = str(raw.get('thread_policy', 'half of available logical CPUs'))
    summary['full_geant4_report_retained'] = True
    summary['photon_boundary_tally_schema'] = str(raw.get('photon_boundary_tally_schema', 'not available'))
    summary['boundary_classification_error_count'] = int(raw.get('boundary_classification_error_count', -1))
    summary['photon_boundary_tally_consistent'] = bool(raw.get('photon_boundary_tally_consistent', False))
    summary['transport_subrun_id'] = str(getattr(cfg, 'transport_subrun_id', 'not assigned'))
    summary['transport_subrun_purpose'] = str(getattr(cfg, 'transport_subrun_purpose', 'photon transport'))
    summary['transport_subrun_seed'] = int(raw.get('random_seed', getattr(cfg, 'random_seed', 0) or 0))
    summary['transport_subrun_histories'] = int(getattr(cfg, 'num_photons', 0) or 0)
    summary['transport_subrun_thickness_cm'] = float(getattr(cfg, 'wall_thickness_cm', math.nan))
    _v7_metrics = _v7_all_photon_air_kerma_metrics(raw, cfg)
    summary.update(_v7_metrics)
    results.update({key: value for key, value in _v7_metrics.items() if key.startswith('all_photon_')})

    if cfg.verbose:
        print_summary(results)
        _native_photon_print_obj_info(cfg, obj_info)
    return results



# ==================== v7.2.4 correctness + auditability override ====================
# 2026-08-11-v7.2.4b-primary-history-semantics-hotfix
import hashlib as _v724_hashlib
from copy import deepcopy as _v724_deepcopy

ADAPTIVE_REFINEMENT_STEPS_CM = (5.0,)
_v724_original_all_photon_metrics = _v7_all_photon_air_kerma_metrics
_v724_original_native_photon_fill_tallies = _native_photon_fill_tallies


def _v7_all_photon_air_kerma_metrics(raw: Dict[str, Any], cfg: Any) -> Dict[str, Any]:
    out = _v724_original_all_photon_metrics(raw, cfg)
    transmitted = int(raw.get('transmitted_count', 0) or 0)
    primary = int(raw.get('primary_photon_forward_crossing_count', -1))
    secondary = int(raw.get('secondary_photon_forward_crossing_count', -1))
    all_gamma = int(raw.get('all_photon_forward_crossing_count', 0) or 0)
    histories = int(raw.get('histories_with_any_forward_gamma_crossing', -1))
    boundary_errors = int(raw.get('boundary_classification_error_count', -1))
    consistent = bool(raw.get('photon_boundary_tally_consistent', False))

    if primary < 0 or secondary < 0 or histories < 0 or boundary_errors < 0:
        raise RuntimeError('v7.2.4 native photon reconciliation fields are missing; rebuild the native photon worker.')
    if primary != transmitted:
        raise RuntimeError(
            f'Photon tally invariant failed: primary forward crossings {primary} != transmitted histories {transmitted}'
        )
    if all_gamma != primary + secondary:
        raise RuntimeError(
            f'Photon tally invariant failed: all crossings {all_gamma} != primary + secondary {primary + secondary}'
        )
    if all_gamma < transmitted:
        raise RuntimeError(
            f'Photon tally invariant failed: all crossings {all_gamma} < transmitted histories {transmitted}'
        )
    if histories < transmitted or histories > all_gamma:
        raise RuntimeError(
            f'Photon tally invariant failed: histories-with-gamma-crossing {histories} is outside [{transmitted}, {all_gamma}]'
        )
    if boundary_errors != 0 or not consistent:
        raise RuntimeError(
            f'Photon tally invariant failed: boundary errors={boundary_errors}, native consistency flag={consistent}'
        )

    source_energy = float(raw.get('source_primary_energy_mev_sum', 0.0) or 0.0)
    global_edep = float(raw.get('global_edep_mev_sum', 0.0) or 0.0)
    world_escape = float(raw.get('world_escape_kinetic_energy_mev_sum', 0.0) or 0.0)
    residual = float(raw.get('global_energy_residual_mev_sum', 0.0) or 0.0)
    residual_fraction = float(raw.get('global_energy_residual_fraction', 0.0) or 0.0)

    out.update({
        'primary_photon_forward_crossing_count': primary,
        'secondary_photon_forward_crossing_count': secondary,
        'histories_with_any_forward_gamma_crossing': histories,
        'boundary_classification_error_count': boundary_errors,
        'photon_boundary_tally_consistent': consistent,
        'global_energy_source_total_mev': source_energy,
        'global_energy_deposit_mev': global_edep,
        'global_energy_world_escape_kinetic_mev': world_escape,
        'global_energy_residual_mev': residual,
        'global_energy_residual_fraction': residual_fraction,
        'global_energy_residual_abs_max_mev': float(raw.get('global_energy_residual_abs_max_mev', 0.0) or 0.0),
        'global_energy_accounting_definition': str(raw.get('global_energy_accounting_definition', 'not available')),
        'photon_tally_reconciliation_status': 'PASS',
    })
    return out


def _native_photon_fill_tallies(cfg, raw: Dict[str, Any]) -> SimulationTallies:
    # v7.3k.5: fail early and explicitly if the native binary was not rebuilt
    # with the boundary-tally schema required by the v7.2.4+ reconciliation
    # contract.  The previous generic invariant error obscured a packaging /
    # binary-source mismatch and looked like a transport-physics failure.
    _required_boundary_fields = (
        'primary_photon_forward_crossing_count',
        'secondary_photon_forward_crossing_count',
        'histories_with_any_forward_gamma_crossing',
        'boundary_classification_error_count',
        'photon_boundary_tally_consistent',
        'photon_boundary_tally_schema',
    )
    _missing_boundary_fields = [key for key in _required_boundary_fields if key not in raw]
    if _missing_boundary_fields:
        raise RuntimeError(
            'Native photon result is missing the v7.3k.6 primary-fate/re-entry boundary schema fields: '
            + ', '.join(_missing_boundary_fields)
            + '. Rebuild native_photon_worker from the installed source before running transport.'
        )
    if str(raw.get('photon_boundary_tally_schema', '')) != 'v7.3k.6-primary-ever-downstream-boundary-v3':
        raise RuntimeError(
            'Native photon boundary-tally schema mismatch: '
            f"{raw.get('photon_boundary_tally_schema')!r}; expected "
            "'v7.3k.6-primary-ever-downstream-boundary-v3'. Rebuild native_photon_worker."
        )
    tallies = _v724_original_native_photon_fill_tallies(cfg, raw)
    tallies.primary_photon_forward_crossing_count = int(raw.get('primary_photon_forward_crossing_count', -1))
    tallies.secondary_photon_forward_crossing_count = int(raw.get('secondary_photon_forward_crossing_count', -1))
    tallies.histories_with_any_forward_gamma_crossing = int(raw.get('histories_with_any_forward_gamma_crossing', -1))
    tallies.boundary_classification_error_count = int(raw.get('boundary_classification_error_count', -1))
    tallies.photon_boundary_tally_consistent = bool(raw.get('photon_boundary_tally_consistent', False))
    tallies.global_energy_source_total_mev = float(raw.get('source_primary_energy_mev_sum', 0.0) or 0.0)
    tallies.global_energy_deposit_mev = float(raw.get('global_edep_mev_sum', 0.0) or 0.0)
    tallies.global_energy_world_escape_kinetic_mev = float(raw.get('world_escape_kinetic_energy_mev_sum', 0.0) or 0.0)
    tallies.global_energy_residual_mev = float(raw.get('global_energy_residual_mev_sum', 0.0) or 0.0)
    tallies.global_energy_residual_fraction = float(raw.get('global_energy_residual_fraction', 0.0) or 0.0)
    tallies.global_energy_residual_abs_max_mev = float(raw.get('global_energy_residual_abs_max_mev', 0.0) or 0.0)

    if tallies.primary_photon_forward_crossing_count != int(tallies.transmitted_count):
        raise RuntimeError(
            'v7.3k.6 primary-forward/transmitted invariant failed after native result import: '
            f'primary={tallies.primary_photon_forward_crossing_count}, transmitted={int(tallies.transmitted_count)}, '
            f'all={int(tallies.all_photon_forward_crossing_count)}, secondary={tallies.secondary_photon_forward_crossing_count}. '
            'This indicates a native boundary-classification defect, not a Geant4 transport failure.'
        )
    if int(tallies.all_photon_forward_crossing_count) != (
        tallies.primary_photon_forward_crossing_count + tallies.secondary_photon_forward_crossing_count
    ):
        raise RuntimeError(
            'v7.3k.6 all-photon primary+secondary invariant failed after native result import: '
            f'all={int(tallies.all_photon_forward_crossing_count)}, primary={tallies.primary_photon_forward_crossing_count}, '
            f'secondary={tallies.secondary_photon_forward_crossing_count}, transmitted={int(tallies.transmitted_count)}. '
            'The native worker must count primary and secondary crossings directly at the same downstream boundary event.'
        )
    if tallies.boundary_classification_error_count != 0 or not tallies.photon_boundary_tally_consistent:
        raise RuntimeError(
            'v7.3k.6 native photon boundary consistency flag failed after result import: '
            f'boundary_errors={tallies.boundary_classification_error_count}, '
            f'consistent={tallies.photon_boundary_tally_consistent}.'
        )
    return tallies


def _v724_target(cfg: Any) -> float:
    for name in (
        'required_max_photon_transmission', 'photon_transmission_target', 'target_transmission',
        'required_transmission', 'max_transmission', 'engineering_transmission_target',
    ):
        try:
            value = float(getattr(cfg, name))
        except Exception:
            continue
        if np.isfinite(value) and 0.0 < value <= 1.0:
            return value
    try:
        value = float(os.environ.get('PHOTON_CONCRETE_REQUIRED_TRANSMISSION', 'nan'))
    except Exception:
        value = np.nan
    return value if np.isfinite(value) and 0.0 < value <= 1.0 else np.nan


def _v724_target_for_thickness(cfg: Any, thickness_cm: float, fallback: float = np.nan) -> float:
    """k.10d thickness-specific B_gamma,max for d_wall(t)=face+t."""
    try:
        face=float(getattr(cfg,'_v73k10d_wall_face_distance_m',getattr(cfg,'wall_face_distance_m',None)))
        k=float(getattr(cfg,'_v73k10d_photon_target_constant_m2'))
        t=float(thickness_cm)
        if np.isfinite(face) and face>0.0 and np.isfinite(k) and k>0.0 and np.isfinite(t) and t>=0.0:
            return float(min(1.0,k*(face+t/100.0)**2))
    except Exception:
        pass
    try:
        f=float(fallback)
    except Exception:
        f=np.nan
    return f if np.isfinite(f) and 0.0<f<=1.0 else _v724_target(cfg)


def _v724_seed(base_seed: int, stage: str, thickness_cm: float) -> int:
    payload = f'{int(base_seed)}|{stage}|{float(thickness_cm):.12g}'.encode('utf-8')
    value = int.from_bytes(_v724_hashlib.sha256(payload).digest()[:8], 'big')
    return 1 + (value % 2147483645)


def _v724_first_observed_candidate_bracket(records: Dict[float, Dict[str, Any]], target: float, cfg: Any = None):
    if not np.isfinite(target) and cfg is None:
        return None
    ordered = sorted(records)
    candidate_position = None
    for i, thickness in enumerate(ordered):
        response = float(records[thickness].get('all_photon_air_kerma_response', np.nan))
        target_i = target
        if cfg is not None:
            # Self-contained on purpose: inherited k.8 regression AST-extracts
            # this helper without adjacent k.10d helper functions.
            try:
                face=float(getattr(cfg,'_v73k10d_wall_face_distance_m',getattr(cfg,'wall_face_distance_m',None)))
                k=float(getattr(cfg,'_v73k10d_photon_target_constant_m2'))
                tt=float(thickness)
                if np.isfinite(face) and face>0.0 and np.isfinite(k) and k>0.0 and np.isfinite(tt) and tt>=0.0:
                    target_i=float(min(1.0,k*(face+tt/100.0)**2))
            except Exception:
                target_i=target
        if np.isfinite(response) and np.isfinite(target_i) and response <= target_i:
            candidate_position = i
            break
    if candidate_position is None:
        return None
    hi = ordered[candidate_position]
    lo = ordered[candidate_position - 1] if candidate_position > 0 else hi
    return (float(lo), float(hi))


def _v724_refinement_points(lo: float, hi: float, max_step_cm: float) -> List[float]:
    lo = float(lo)
    hi = float(hi)
    step = float(max_step_cm)
    if not (np.isfinite(lo) and np.isfinite(hi) and np.isfinite(step) and step > 0.0 and hi > lo):
        return []
    if hi - lo <= step + 1.0e-9:
        return []
    points: List[float] = []
    value = lo + step
    while value < hi - 1.0e-9:
        points.append(float(value))
        value += step
    return points


def run_thickness_sweep(cfg) -> Dict[str, np.ndarray]:
    """v7.3k.3: 25-cm coarse scan plus C-25..C+15 high-stat 5-cm refinement."""
    # Local k.10d target helper is embedded here so inherited AST-extraction
    # regressions execute the same candidate-dependent distance logic.
    def _target_at(thickness_cm, fallback):
        try:
            face=float(getattr(cfg,'_v73k10d_wall_face_distance_m',getattr(cfg,'wall_face_distance_m',None)))
            k=float(getattr(cfg,'_v73k10d_photon_target_constant_m2'))
            tt=float(thickness_cm)
            if np.isfinite(face) and face>0.0 and np.isfinite(k) and k>0.0 and np.isfinite(tt) and tt>=0.0:
                return float(min(1.0,k*(face+tt/100.0)**2))
        except Exception:
            pass
        try:
            return float(fallback)
        except Exception:
            return np.nan
    print('\n================ CONCRETE THICKNESS SWEEP v7.3k.8 ================')
    print('Every thickness is a direct full-physical-spectrum Geant4 run; coarse points in the fine window are rerun/promoted when FINAL_PHOTONS is higher.')
    print('Adaptive resolution sequence: 25 cm coarse grid -> high-stat 5 cm window from C-25 cm through C+15 cm.')
    print('No interpolation, TVL, suffix-maximum forcing, or extrapolation selects thickness.')
    print('=================================================================\n')

    requested = np.asarray(cfg.sweep_thicknesses_cm, dtype=float).ravel()
    if requested.size == 0:
        raise ValueError('sweep_thicknesses_cm must contain at least one thickness')
    if np.any(~np.isfinite(requested)):
        raise ValueError('sweep_thicknesses_cm contains a non-finite value')
    requested = np.unique(np.round(requested, 12))
    requested.sort()

    # >>> v7.3b STANDARD 50-400 CM COARSE GRID >>>
    # The automatic concrete search contract is 50..400 cm in 25 cm coarse
    # increments, followed directly by 5 cm refinement in the first observed
    # coarse candidate bracket. Custom/nonstandard requested grids are not
    # silently rewritten.
    if requested.size >= 2:
        _v73b_diffs = np.diff(requested)
        _v73b_standard_25 = (
            abs(float(requested[0]) - 50.0) <= 1.0e-9
            and np.all(np.isfinite(_v73b_diffs))
            and np.all(np.abs(_v73b_diffs - 25.0) <= 1.0e-9)
        )
        if _v73b_standard_25 and not bool(getattr(cfg, 'sweep_exact_requested_grid', False)) and float(requested[-1]) < 400.0 - 1.0e-9:
            requested = np.arange(50.0, 400.0 + 1.0e-9, 25.0, dtype=float)
            try:
                cfg.sweep_thicknesses_cm = requested.copy()
            except Exception:
                pass
    # <<< v7.3b STANDARD 50-400 CM COARSE GRID <<<
    n_point = int(cfg.sweep_num_photons_per_point)
    if n_point <= 0:
        raise ValueError('sweep_num_photons_per_point must be positive')
    exact_requested_grid = bool(getattr(cfg, 'sweep_exact_requested_grid', False))
    stage_override = str(getattr(cfg, 'sweep_stage_override', '') or '').strip()
    purpose_override = str(getattr(cfg, 'sweep_purpose_override', '') or '').strip()

    target = _v724_target(cfg)
    base_seed = int(getattr(cfg, 'random_seed', 1) or 1)
    parent_run_id = str(
        getattr(cfg, 'run_id', '') or getattr(cfg, '_ncrp_run_id', '') or
        os.environ.get('NCRP_RUN_ID', '') or f'photon-sweep-seed-{base_seed}'
    )
    records: Dict[float, Dict[str, Any]] = {}
    # v7.3k.8 keeps an immutable execution ledger.  `records` remains the
    # authoritative per-thickness realization (higher-stat fine runs may
    # supersede a coarse value for design), while every actually executed
    # Geant4 subrun is retained below for provenance/accounting.
    run_ledger: List[Dict[str, Any]] = []
    coarse_records: Dict[float, Dict[str, Any]] = {}
    stage_records: List[Dict[str, Any]] = []

    def run_point(thickness: float, stage: str, purpose: str, histories: int = None, allow_replace: bool = False):
        key = float(round(float(thickness), 12))
        n_run = int(n_point if histories is None else histories)
        if n_run <= 0:
            raise ValueError('Photon sweep point histories must be positive')
        if key in records and not allow_replace:
            return records[key]
        if key in records and allow_replace and n_run <= int(records[key].get('primaries', 0) or 0):
            return records[key]
        seed = _v724_seed(base_seed, stage, key)
        point_cfg = _v724_deepcopy(cfg)
        point_cfg.wall_thickness_cm = key
        # v7.3k.7 correctness: use the requested stage history count, not the
        # coarse default.  This is the direct-call counterpart of the GUI's
        # v7.3k.3 high-stat fine-stage promotion.
        point_cfg.num_photons = n_run
        point_cfg.random_seed = seed
        # Coarse/fine thickness scans are ALWAYS direct full-physical-spectrum
        # analog runs.  The VR supervisor is supplemental and cannot replace
        # any concrete scan point.
        point_cfg._v73k_internal_photon = True
        point_cfg._v73k_inline_importance_enabled = False
        # Explicitly clear the legacy per-group custom-source hook as well.
        # Even if a caller accidentally carries VR scratch attributes, a
        # concrete sweep point must reload the configured full physical source.
        point_cfg._v73k_spectrum_energies_mev = None
        point_cfg._v73k_spectrum_probabilities = None
        point_cfg.n_visual_tracks = 0
        point_cfg.export_obj = False
        point_cfg.verbose = False
        point_cfg.run_thickness_sweep = False
        subrun_id = f'{parent_run_id}:photon-{stage}:{key:g}cm'
        try:
            point_cfg.run_id = subrun_id
            point_cfg._ncrp_run_id = subrun_id
        except Exception:
            pass

        print(f'\n===== PHOTON {stage.upper()} SUBRUN: {key:g} cm, N={n_run:,}, seed={seed} =====', flush=True)
        point_result = run_simulation(point_cfg)
        summary = point_result.get('summary', {})
        if not bool(summary.get('photon_boundary_tally_consistent', False)):
            raise RuntimeError(f'{subrun_id}: native photon boundary tally did not pass v7.2.4 reconciliation.')
        native_seed = summary.get('random_seed', seed)
        try:
            native_seed_int = int(native_seed)
        except Exception as exc:
            raise RuntimeError(
                f'{subrun_id}: native result random_seed is not parseable: {native_seed!r}.'
            ) from exc
        if native_seed_int != seed:
            raise RuntimeError(f'{subrun_id}: requested seed {seed} but native result reports {native_seed}.')

        record = {
            'thickness_cm': key,
            'primaries': n_run,
            'transmitted_fraction': float(summary.get('transmitted_fraction', np.nan)),
            'transmitted_fraction_sigma': float(summary.get('transmitted_fraction_sigma', np.nan)),
            'transmitted_count': int(summary.get('transmitted_count', 0) or 0),
            'effective_mu_per_cm': float(summary.get('effective_mu_per_cm', np.nan)),
            'effective_mu_per_cm_sigma': float(summary.get('effective_mu_per_cm_sigma', np.nan)),
            'all_photon_air_kerma_response': float(summary.get('all_photon_air_kerma_response', np.nan)),
            'all_photon_air_kerma_response_sem': float(summary.get('all_photon_air_kerma_response_sem', np.nan)),
            'all_photon_forward_crossing_count': int(summary.get('all_photon_forward_crossing_count', 0) or 0),
            'primary_photon_forward_crossing_count': int(summary.get('primary_photon_forward_crossing_count', -1)),
            'secondary_photon_forward_crossing_count': int(summary.get('secondary_photon_forward_crossing_count', -1)),
            'histories_with_any_forward_gamma_crossing': int(summary.get('histories_with_any_forward_gamma_crossing', -1)),
            'boundary_classification_error_count': int(summary.get('boundary_classification_error_count', -1)),
            'photon_boundary_tally_consistent': bool(summary.get('photon_boundary_tally_consistent', False)),
            'photon_boundary_tally_schema': str(summary.get('photon_boundary_tally_schema', 'not available')),
            'source_air_kerma_weight_mean': float(summary.get('source_air_kerma_weight_mean', np.nan)),
            'max_air_kerma_weight_per_photon': float(summary.get('max_air_kerma_weight_per_photon', np.nan)),
            'statistical_verification_upper': float(summary.get('all_photon_air_kerma_statistical_verification_upper', np.nan)),
            'statistically_verified': bool(summary.get('all_photon_air_kerma_statistically_verified', False)),
            'random_seed': seed,
            'subrun_id': subrun_id,
            'subrun_stage': stage,
            'subrun_purpose': purpose,
            'source_sampling_role': 'FULL_PHYSICAL_SPECTRUM_ANALOG',
        }
        if record['primary_photon_forward_crossing_count'] != record['transmitted_count']:
            raise RuntimeError(f'{subrun_id}: primary forward crossings != transmitted primary histories.')
        if record['all_photon_forward_crossing_count'] != (
            record['primary_photon_forward_crossing_count'] + record['secondary_photon_forward_crossing_count']
        ):
            raise RuntimeError(f'{subrun_id}: all-photon crossings != primary + secondary crossings.')
        if record['boundary_classification_error_count'] != 0:
            raise RuntimeError(f'{subrun_id}: boundary-classification error count is nonzero.')
        # Preserve every physical subrun exactly once in the immutable ledger.
        ledger_record = dict(record)
        ledger_record['execution_index'] = len(run_ledger) + 1
        run_ledger.append(ledger_record)
        if 'coarse' in str(stage).lower():
            coarse_records[key] = dict(record)
        records[key] = record
        print(
            f"SWEEP RESULT: thickness={key:g} cm | response={record['all_photon_air_kerma_response']:.8e} | "
            f"SEM={record['all_photon_air_kerma_response_sem']:.3e} | "
            f"gamma crossings={record['all_photon_forward_crossing_count']:,} "
            f"(primary={record['primary_photon_forward_crossing_count']:,}, secondary={record['secondary_photon_forward_crossing_count']:,}) | "
            f"primary-history T={record['transmitted_fraction']:.8e}",
            flush=True,
        )
        return record

    coarse_like = requested.size >= 2 and float(np.median(np.diff(requested))) >= 20.0 - 1.0e-9
    initial_stage = stage_override or ('coarse-25cm' if coarse_like else 'requested-grid')
    initial_purpose = purpose_override or ('coarse candidate bracketing' if coarse_like else 'explicit requested sweep grid')
    for thickness in requested:
        run_point(float(thickness), initial_stage, initial_purpose)

    if coarse_like and np.isfinite(target) and not exact_requested_grid:
        # v7.3k.3 direct-call fallback: same production policy as the GUI
        # orchestrator.  The first coarse observed candidate C anchors a local
        # 5-cm window from C-25 through C+15.  Fine-stage N is FINAL_PHOTONS
        # when that exceeds the coarse N.  Existing coarse points in the window
        # are rerun independently and their higher-N realization becomes the
        # authoritative sweep value.
        # Candidate routing MUST use only actual low-stat coarse realizations.
        # Fine-stage values are never allowed to create or move the coarse anchor.
        bracket_before = _v724_first_observed_candidate_bracket(coarse_records, target, cfg)
        if bracket_before is None:
            stage_records.append({
                'stage': 'refine-5cm-highstat',
                'requested_step_cm': 5.0,
                'coarse_candidate_cm': None,
                'fine_window_cm': None,
                'fine_points_cm': [],
                'rerun_coarse_points_cm': [],
                'status': 'not run: no observed coarse candidate',
            })
        else:
            candidate = float(bracket_before[1])
            coarse_candidate_record = coarse_records.get(float(round(candidate, 12)), {})
            lower_bound = max(float(np.min(requested)), candidate - 25.0)
            upper_bound = min(float(np.max(requested)), candidate + 15.0)
            fine_n = max(n_point, int(getattr(cfg, 'fine_screen_num_photons_per_point', getattr(cfg, 'fine_num_photons_per_point', getattr(cfg, 'final_num_photons', getattr(cfg, 'num_photons', n_point)))) or n_point))
            step_cm = 5.0
            tol = 1.0e-9
            k_lo = int(math.ceil((lower_bound - candidate - tol) / step_cm))
            k_hi = int(math.floor((upper_bound - candidate + tol) / step_cm))
            fine_points = [float(round(candidate + k * step_cm, 12)) for k in range(k_lo, k_hi + 1)]
            fine_points = [x for x in fine_points if x >= lower_bound - tol and x <= upper_bound + tol]
            rerun_points = [x for x in fine_points if x in records and fine_n > int(records[x].get('primaries', 0) or 0)]
            for thickness in fine_points:
                run_point(
                    float(thickness),
                    'refine-5cm-highstat',
                    'high-stat 5 cm candidate-centered refinement: 25 cm below through 15 cm above coarse candidate',
                    histories=fine_n,
                    allow_replace=True,
                )
            stage_records.append({
                'stage': 'refine-5cm-highstat',
                'requested_step_cm': step_cm,
                'coarse_candidate_cm': candidate,
                'fine_window_cm': [lower_bound, upper_bound],
                'fine_points_cm': fine_points,
                'rerun_coarse_points_cm': rerun_points,
                'coarse_histories_per_point': n_point,
                'fine_histories_per_point': fine_n,
                'coarse_candidate_source_subrun_id': str(coarse_candidate_record.get('subrun_id', 'not available')),
                'coarse_candidate_source_response': float(coarse_candidate_record.get('all_photon_air_kerma_response', np.nan)),
                'status': 'complete',
            })

    ordered = sorted(records)
    def arr(name, dtype=float):
        return np.asarray([records[t][name] for t in ordered], dtype=dtype)

    return {
        'thicknesses_cm': np.asarray(ordered, dtype=float),
        'primaries': arr('primaries', np.int64),
        'transmitted_count': arr('transmitted_count', np.int64),
        'transmitted_fraction': arr('transmitted_fraction', float),
        'transmitted_fraction_sigma': arr('transmitted_fraction_sigma', float),
        'effective_mu_per_cm': arr('effective_mu_per_cm', float),
        'effective_mu_per_cm_sigma': arr('effective_mu_per_cm_sigma', float),
        'all_photon_air_kerma_response': arr('all_photon_air_kerma_response', float),
        'all_photon_air_kerma_response_sem': arr('all_photon_air_kerma_response_sem', float),
        'all_photon_forward_crossing_count': arr('all_photon_forward_crossing_count', np.int64),
        'primary_photon_forward_crossing_count': arr('primary_photon_forward_crossing_count', np.int64),
        'secondary_photon_forward_crossing_count': arr('secondary_photon_forward_crossing_count', np.int64),
        'histories_with_any_forward_gamma_crossing': arr('histories_with_any_forward_gamma_crossing', np.int64),
        'boundary_classification_error_count': arr('boundary_classification_error_count', np.int64),
        'photon_boundary_tally_consistent': arr('photon_boundary_tally_consistent', bool),
        'photon_boundary_tally_schema': np.asarray([records[t]['photon_boundary_tally_schema'] for t in ordered], dtype=object),
        'source_air_kerma_weight_mean': arr('source_air_kerma_weight_mean', float),
        'max_air_kerma_weight_per_photon': arr('max_air_kerma_weight_per_photon', float),
        'statistical_verification_upper': arr('statistical_verification_upper', float),
        'statistically_verified': arr('statistically_verified', bool),
        'random_seed': arr('random_seed', np.int64),
        'subrun_id': np.asarray([records[t]['subrun_id'] for t in ordered], dtype=object),
        'subrun_stage': np.asarray([records[t]['subrun_stage'] for t in ordered], dtype=object),
        'subrun_purpose': np.asarray([records[t]['subrun_purpose'] for t in ordered], dtype=object),
        'source_sampling_role': np.asarray([records[t]['source_sampling_role'] for t in ordered], dtype=object),
        # Compatibility aliases used by shielding/report merge code.
        'transport_subrun_ids': np.asarray([records[t]['subrun_id'] for t in ordered], dtype=object),
        'transport_subrun_purposes': np.asarray([records[t]['subrun_purpose'] for t in ordered], dtype=object),
        'transport_subrun_seeds': arr('random_seed', np.int64),
        'adaptive_stage_records': stage_records,
        'subrun_ledger': [dict(row) for row in run_ledger],
        'coarse_subrun_ledger': [dict(row) for row in run_ledger if 'coarse' in str(row.get('subrun_stage','')).lower()],
        'fine_subrun_ledger': [dict(row) for row in run_ledger if ('fine' in str(row.get('subrun_stage','')).lower() or 'refine' in str(row.get('subrun_stage','')).lower())],
        'coarse_candidate_cm': (float(stage_records[-1].get('coarse_candidate_cm')) if stage_records and stage_records[-1].get('coarse_candidate_cm') is not None else None),
        'coarse_candidate_source_subrun_id': (str(stage_records[-1].get('coarse_candidate_source_subrun_id')) if stage_records and stage_records[-1].get('coarse_candidate_source_subrun_id') else None),
        'adaptive_strategy': (
            'exact requested grid only; adaptive orchestration is handled by the GUI workflow'
            if exact_requested_grid else
            '25 cm coarse grid -> high-stat 5 cm candidate-centered C-25..C+15 refinement'
        ),
        'sweep_exact_requested_grid': exact_requested_grid,
        'adaptive_target': target,
        'adaptive_target_is_thickness_dependent': bool(getattr(cfg,'_v73k10d_dynamic_wall_target',False)),
        'required_max_transmission_by_thickness': np.asarray([_target_at(t,target) for t in ordered],dtype=float),
        'protected_distance_m_by_thickness': np.asarray([(float(getattr(cfg,'_v73k10d_wall_face_distance_m'))+float(t)/100.0) if hasattr(cfg,'_v73k10d_wall_face_distance_m') else np.nan for t in ordered],dtype=float),
        'wall_face_distance_m': float(getattr(cfg,'_v73k10d_wall_face_distance_m')) if hasattr(cfg,'_v73k10d_wall_face_distance_m') else np.nan,
        'distance_geometry_policy': 'd_wall(t)=wall_face+t; thickness-specific B_gamma,max(t)',
        'physical_selection_tally': 'all_photon_downstream_dry_air_kerma_response',
        'primary_history_fraction_role': 'diagnostic only',
        'uses_interpolation': False,
        'uses_tvl': False,
        'uses_extrapolation': False,
        'v724_feature_marker': '2026-08-11-v7.2.4b-primary-history-semantics-hotfix',
    }
# ================== end v7.2.4 correctness + auditability override ==================



# >>> v7.3d SUPERVISORY SEARCH + INTERVAL HARDENING PATCH >>>
# High-statistics supervisory validation only. The validated native canonical
# downstream photon scorer is not changed.
from v73d_supervisory import (
    finite as _v73d_finite,
    should_supervise_photon as _v73d_should_supervise_photon,
    supervise_photon_validation as _v73d_supervise_photon_validation,
)

_v73d_previous_air_kerma_metrics = _v7_all_photon_air_kerma_metrics


def _v7_all_photon_air_kerma_metrics(raw, cfg):
    # Expose the already-scored additive sufficient statistics so independent
    # validation batches can be pooled exactly without changing native scoring.
    out = _v73d_previous_air_kerma_metrics(raw, cfg)
    out.update({
        'v73d_source_air_kerma_weight_sum': float(raw.get('source_air_kerma_weight_sum', 0.0)),
        'v73d_source_air_kerma_weight_sumsq': float(raw.get('source_air_kerma_weight_sumsq', 0.0)),
        'v73d_exit_air_kerma_weight_sum': float(raw.get('exit_air_kerma_weight_sum', 0.0)),
        'v73d_exit_air_kerma_weight_sumsq': float(raw.get('exit_air_kerma_weight_sumsq', 0.0)),
        'v73d_source_exit_air_kerma_weight_cross_sum': float(raw.get('source_exit_air_kerma_weight_cross_sum', 0.0)),
        # v7.3k.2 native rare-event influence audit.  These are the largest
        # per-primary downstream dry-air event weights retained by the native
        # worker; they are diagnostic only and never alter transport/scoring.
        'v73k2_top_exit_air_kerma_event_weights': [
            float(v) for v in raw.get('top_exit_air_kerma_event_weights', [])
            if np.isfinite(float(v)) and float(v) > 0.0
        ],
        'v73k2_top_exit_air_kerma_event_weight_cap': int(raw.get('top_exit_air_kerma_event_weight_cap', 0) or 0),
        # v7.3k.7 inline native source-importance payload.  Response moments
        # above are already p/q-weighted by the native worker; these additional
        # fields retain group/support/ESS/batch auditability without per-group
        # Geant4 process launches.
        'v73k7_importance_sampling_enabled': bool(raw.get('importance_sampling_enabled', False)),
        'v73k7_importance_sampling_schema': str(raw.get('importance_sampling_schema', 'ANALOG_PHYSICAL_SOURCE')),
        'v73k7_importance_group_count': int(raw.get('importance_group_count', 0) or 0),
        'v73k7_importance_statistical_batch_count': int(raw.get('importance_statistical_batch_count', 0) or 0),
        'v73k7_importance_likelihood_weight_sum': float(raw.get('importance_likelihood_weight_sum', 0.0) or 0.0),
        'v73k7_importance_likelihood_weight_sumsq': float(raw.get('importance_likelihood_weight_sumsq', 0.0) or 0.0),
        'v73k7_importance_min_likelihood_weight': float(raw.get('importance_min_likelihood_weight', 1.0) or 1.0),
        'v73k7_importance_max_likelihood_weight': float(raw.get('importance_max_likelihood_weight', 1.0) or 1.0),
        'v73k7_importance_group_stats': list(raw.get('importance_group_stats', []) or []),
        'v73k7_importance_statistical_batches': list(raw.get('importance_statistical_batches', []) or []),
    })
    return out


_v73d_previous_run_simulation = run_simulation


def run_simulation(cfg):
    if not _v73d_should_supervise_photon(cfg):
        return _v73d_previous_run_simulation(cfg)

    target = None
    if '_v724_target' in globals():
        try:
            target = _v73d_finite(_v724_target(cfg))
        except Exception:
            target = None
    if target is None:
        try:
            target = _v73d_finite(os.environ.get('PHOTON_CONCRETE_REQUIRED_TRANSMISSION'))
        except Exception:
            target = None
    if target is None:
        # Fail safely: without a frozen target there is no supervisory decision.
        return _v73d_previous_run_simulation(cfg)

    results, audit = _v73d_supervise_photon_validation(
        _v73d_previous_run_simulation,
        cfg,
        target,
    )
    print('\n================ v7.3d PHOTON SUPERVISORY RESULT ================', flush=True)
    print('Policy:', audit.get('policy'), flush=True)
    print('Verified passing thickness:', audit.get('verified_passing_thickness_cm'), flush=True)
    print('Minimum verified compliant thickness:', audit.get('minimum_verified_compliant_thickness_cm'), flush=True)
    print('Unresolved thinner candidates:', audit.get('unresolved_thinner_candidate_thicknesses_cm'), flush=True)
    print('Range exhausted:', audit.get('range_exhausted'), flush=True)
    print('==================================================================\n', flush=True)
    return results

V73D_PHOTON_VALIDATION_CONTRACT = (
    'cumulative 1x -> 2x -> 4x independent high-statistics validation at each actual candidate (added batches 1x,+1x,+2x); '
    'if unresolved/fail then next thicker 5 cm point; hard maximum 400 cm'
)
# <<< v7.3d SUPERVISORY SEARCH + INTERVAL HARDENING PATCH <<<


# >>> v7.3e ENERGY-CONSISTENT SUPERVISORY CORRECTNESS PATCH >>>
from v73e_correctness import should_supervise_photon as _v73e_should_supervise_photon, supervise_photon_validation as _v73e_supervise_photon_validation, finite as _v73e_finite

# Bypass the v7.3d supervisor while retaining the validated pre-supervisory native adapter.
_v73e_native_photon_run = _v73d_previous_run_simulation

def run_simulation(cfg):
    if not _v73e_should_supervise_photon(cfg):
        return _v73e_native_photon_run(cfg)
    target=None
    if '_v724_target' in globals():
        try: target=_v73e_finite(_v724_target(cfg))
        except Exception: target=None
    if target is None:
        try: target=_v73e_finite(os.environ.get('PHOTON_CONCRETE_REQUIRED_TRANSMISSION'))
        except Exception: target=None
    if target is None: return _v73e_native_photon_run(cfg)
    results,audit=_v73e_supervise_photon_validation(_v73e_native_photon_run,cfg,target)
    print('\n================ v7.3e PHOTON SUPERVISORY RESULT ================',flush=True)
    print('Verified passing thickness:',audit.get('verified_passing_thickness_cm'),flush=True)
    print('Minimum verified compliant thickness:',audit.get('minimum_verified_compliant_thickness_cm'),flush=True)
    print('Unresolved thinner candidates:',audit.get('unresolved_thinner_candidate_thicknesses_cm'),flush=True)
    print('Plot-bearing native batch:',audit.get('plot_bearing_batch_thickness_cm'),'cm /',audit.get('plot_bearing_batch_histories'),'histories',flush=True)
    print('Cumulative decision histories:',audit.get('final_cumulative',{}).get('histories'),flush=True)
    print('==================================================================\n',flush=True)
    return results

V73E_PHOTON_CONTRACT='25 cm coarse -> direct 5 cm refinement; adaptive cumulative high-stat validation; plot-bearing first batch retained; no 15/10 cm stages'
# <<< v7.3e ENERGY-CONSISTENT SUPERVISORY CORRECTNESS PATCH <<<


# >>> v7.3f TEST-BUDGET + REPORT-INTEGRATION CORRECTNESS PATCH >>>
import os as _v73f_os
from v73f_correctness import (
    should_supervise_photon as _v73f_should_supervise_photon,
    supervise_photon_validation as _v73f_supervise_photon_validation,
    finite as _v73f_finite,
)

# Deliberately bypass the v7.3e supervisor while retaining the validated
# pre-supervisory native Python adapter and untouched native C++ scorer.
_v73f_native_photon_run = _v73e_native_photon_run

def run_simulation(cfg):
    if not _v73f_should_supervise_photon(cfg):
        return _v73f_native_photon_run(cfg)
    target = None
    if '_v724_target' in globals():
        try:
            target = _v73f_finite(_v724_target(cfg))
        except Exception:
            target = None
    if target is None:
        try:
            target = _v73f_finite(_v73f_os.environ.get('PHOTON_CONCRETE_REQUIRED_TRANSMISSION'))
        except Exception:
            target = None
    if target is None:
        return _v73f_native_photon_run(cfg)
    results, audit = _v73f_supervise_photon_validation(_v73f_native_photon_run, cfg, target)
    # v7.3i.2: make the supervisory audit an explicit part of the returned result.
    # Report generation must never have to infer the final decision from the
    # plot-bearing batch or reload mutable GUI configuration.
    if not isinstance(results, dict):
        raise TypeError('photon supervisor returned a non-dictionary result')
    summary = results.setdefault('summary', {})
    if not isinstance(summary, dict):
        raise TypeError('photon result summary is not a dictionary')
    summary['v73f_photon_validation_audit'] = audit
    summary['v73h_photon_validation_audit'] = audit
    summary['v73i2_photon_validation_audit'] = audit
    fc = audit.get('final_cumulative', {}) if isinstance(audit, dict) and isinstance(audit.get('final_cumulative'), dict) else {}
    for src, dst in (
        ('response', 'v73e_cumulative_photon_response'),
        ('sem', 'v73e_cumulative_photon_response_sem'),
        ('lower95', 'v73e_cumulative_photon_lower95'),
        ('upper95', 'v73e_cumulative_photon_upper95'),
        ('histories', 'v73e_cumulative_photon_histories'),
        ('crossings', 'v73e_cumulative_photon_crossings'),
    ):
        if src in fc:
            summary[dst] = fc.get(src)
    print('\n================ v7.3f PHOTON SUPERVISORY RESULT ================', flush=True)
    _primary_min = _v73f_finite(audit.get('minimum_verified_compliant_thickness_cm'))
    _primary_verified = _v73f_finite(audit.get('verified_passing_thickness_cm'))
    print('******** PRIMARY RESULT — REQUIRED CONCRETE WALL THICKNESS ********', flush=True)
    if _primary_min is not None:
        print(f'REQUIRED CONCRETE WALL THICKNESS: {_primary_min:.3f} cm — VERIFIED MINIMUM COMPLIANT', flush=True)
    else:
        print('REQUIRED CONCRETE WALL THICKNESS: NOT ESTABLISHED', flush=True)
        if _primary_verified is not None:
            print(f'Verified passing thickness: {_primary_verified:.3f} cm; thinner candidate(s) remain unresolved', flush=True)
    print('PRIMARY RESULT STATUS:', str(audit.get('final_state') or 'STATISTICALLY_UNRESOLVED').replace('_', ' '), flush=True)
    print('********************************************************************', flush=True)
    print('Initial supervisory candidate after coarse/fine screening:', audit.get('initial_low_stat_sweep_candidate_thickness_cm'), 'cm', flush=True)
    print('Initial candidate post-validation state:', audit.get('initial_low_stat_candidate_post_validation_state'), flush=True)
    print('Verified passing thickness:', audit.get('verified_passing_thickness_cm'), flush=True)
    print('Minimum verified compliant thickness:', audit.get('minimum_verified_compliant_thickness_cm'), flush=True)
    print('Unresolved thinner candidates:', audit.get('unresolved_thinner_candidate_thicknesses_cm'), flush=True)
    print('Validation minimum / user hard maximum:', audit.get('minimum_validation_histories'), '/', audit.get('user_hard_max_validation_histories'), flush=True)
    print('Total supervisory validation histories executed:', audit.get('total_supervisory_validation_histories_executed'), flush=True)
    print('Final tested thickness:', audit.get('final_tested_thickness_cm'), 'cm', flush=True)
    print('Final state:', audit.get('final_state'), flush=True)
    print('==================================================================\n', flush=True)
    return results

V73F_PHOTON_CONTRACT = (
    '25 cm coarse -> direct 5 cm only; validation starts at sweep/min count and '
    'may double only to user-selected num_photons absolute cumulative maximum; '
    'unresolved at user max is VALIDATION_BUDGET_EXHAUSTED; no 15/10 cm stages'
)
# <<< v7.3f TEST-BUDGET + REPORT-INTEGRATION CORRECTNESS PATCH <<<


# >>> v7.3i.2 WEIGHTED PHOTON FINAL-VALIDATION REPORT CONTRACT >>>
_v73i2_legacy_attach_final_concrete_validation = attach_final_concrete_validation


def _v73i2_finite(value):
    try:
        x = float(value)
    except Exception:
        return None
    return x if math.isfinite(x) else None


def _v73i2_min_zero_count_wilson_histories(target):
    p = _v73i2_finite(target)
    if p is None or not (0.0 < p < 1.0):
        return None
    z2 = 1.959963984540054 ** 2
    return int(math.ceil(z2 * (1.0 - p) / p))


def _v73i2_photon_audit(results):
    summary = results.get('summary', {}) if isinstance(results, dict) else {}
    if not isinstance(summary, dict):
        return {}
    for key in (
        'v73i2_photon_validation_audit',
        'v73h_photon_validation_audit',
        'v73f_photon_validation_audit',
        'v73e_photon_validation_audit',
        'v73d_photon_validation_audit',
    ):
        value = summary.get(key)
        if isinstance(value, dict):
            return value
    return {}


def attach_final_concrete_validation(requirement: Dict[str, Any], results: Dict[str, Any]) -> Dict[str, Any]:
    """Attach the authoritative weighted-response supervisory state.

    Primary-history transmission/Wilson values remain available only as a
    diagnostic.  The physical concrete decision is the all-photon downstream
    dry-air E*(mu_en/rho) response from the supervisory audit.
    """
    output = dict(requirement) if isinstance(requirement, dict) else {}
    summary = results.get('summary', {}) if isinstance(results, dict) else {}
    if not isinstance(summary, dict):
        summary = {}
    audit = _v73i2_photon_audit(results)
    fc = audit.get('final_cumulative', {}) if isinstance(audit.get('final_cumulative'), dict) else {}

    def choose(*values):
        for value in values:
            x = _v73i2_finite(value)
            if x is not None:
                return x
        return None

    target = choose(audit.get('target'), output.get('target_transmission'), output.get('required_max_photon_transmission'))
    final_t = choose(audit.get('final_tested_thickness_cm'), summary.get('wall_thickness_cm'))
    response = choose(fc.get('response'), summary.get('v73e_cumulative_photon_response'), summary.get('all_photon_air_kerma_response'))
    sem = choose(fc.get('sem'), summary.get('v73e_cumulative_photon_response_sem'), summary.get('all_photon_air_kerma_response_sem'))
    lower = choose(fc.get('lower95'), summary.get('v73e_cumulative_photon_lower95'))
    upper = choose(fc.get('upper95'), summary.get('v73e_cumulative_photon_upper95'), summary.get('all_photon_air_kerma_statistical_verification_upper'))
    histories = int(choose(fc.get('histories'), summary.get('v73e_cumulative_photon_histories'), summary.get('num_photons')) or 0)
    crossings = int(choose(fc.get('crossings'), summary.get('v73e_cumulative_photon_crossings'), summary.get('all_photon_forward_crossing_count')) or 0)
    final_state = str(audit.get('final_state') or fc.get('state') or 'STATISTICALLY_UNRESOLVED').upper()
    statistical_state = str(fc.get('state') or audit.get('final_state') or 'STATISTICALLY_UNRESOLVED').upper()
    verified = choose(audit.get('verified_passing_thickness_cm'))
    minimum_verified = choose(audit.get('minimum_verified_compliant_thickness_cm'))
    global_exhausted = bool(audit.get('global_validation_budget_exhausted'))
    validation_exhausted = bool(audit.get('validation_budget_exhausted'))

    # The count-probability Wilson result is a useful scale diagnostic only.  It
    # is not the pass criterion for the weighted photon response.
    total_diag = int(summary.get('num_photons', 0) or 0)
    transmitted_diag = int(summary.get('transmitted_count', 0) or 0)
    if total_diag > 0 and 0 <= transmitted_diag <= total_diag:
        diag_observed = transmitted_diag / total_diag
        diag_low, diag_high = wilson_interval(transmitted_diag, total_diag)
    else:
        diag_observed = math.nan
        diag_low = math.nan
        diag_high = math.nan

    nmin = _v73i2_min_zero_count_wilson_histories(target)
    if verified is not None:
        status = 'VERIFIED PASS — weighted all-photon supervisory upper response is at/below B_gamma,max at an actually simulated thickness'
    elif statistical_state in ('RESOLVED_FAIL', 'FAIL', 'FAIL_RESOLVED'):
        status = 'RESOLVED FAIL — weighted all-photon supervisory lower response remains above B_gamma,max'
    elif global_exhausted:
        status = 'STATISTICALLY UNRESOLVED — global low-count supervisory validation budget exhausted; no compliant thickness established'
    elif validation_exhausted:
        status = 'STATISTICALLY UNRESOLVED — user-selected per-candidate validation history maximum reached; no compliant thickness established'
    elif final_state in ('RANGE_EXHAUSTED', 'SEARCH_LIMIT_EXHAUSTED', 'NOT_REACHED'):
        status = 'NOT REACHED within actually simulated concrete range — no extrapolation'
    else:
        status = final_state.replace('_', ' ')

    output.update({
        'feature_build': '2026-08-14-v7.3i.2-weighted-photon-supervisory-report',
        'method': 'all-photon downstream dry-air E*(mu_en/rho) response; actual Geant4 points; 25 cm coarse -> direct 5 cm refinement; cumulative user-bounded supervisory validation; no TVL/interpolation/extrapolation',
        'target_transmission': target if target is not None else math.nan,
        'required_max_photon_transmission': target if target is not None else math.nan,
        'minimum_histories_for_zero_count_to_bound_target_95': nmin,
        'required_thickness_95_cm': verified if verified is not None else math.nan,
        'first_compliant_simulated_thickness_cm': minimum_verified if minimum_verified is not None else math.nan,
        'final_validation_thickness_cm': final_t if final_t is not None else math.nan,
        'final_supervisory_tested_thickness_cm': final_t,
        'final_supervisory_cumulative_primary_photons': histories,
        'final_supervisory_all_photon_forward_crossings': crossings,
        'final_supervisory_weighted_response': response,
        'final_supervisory_weighted_response_sem': sem,
        'final_supervisory_weighted_lower95': lower,
        'final_supervisory_weighted_upper95': upper,
        'final_supervisory_statistical_state': statistical_state,
        'verified_passing_thickness_cm': verified,
        'minimum_verified_compliant_thickness_cm': minimum_verified,
        'global_validation_budget_exhausted': global_exhausted,
        'validation_budget_exhausted': validation_exhausted,
        'final_validation_meets_target_95': bool(verified is not None),
        'legacy_primary_history_diagnostic_primary_photons': total_diag,
        'legacy_primary_history_diagnostic_transmitted_photons': transmitted_diag,
        'legacy_primary_history_diagnostic_fraction': diag_observed,
        'legacy_primary_history_diagnostic_wilson95_low': float(diag_low),
        'legacy_primary_history_diagnostic_wilson95_high': float(diag_high),
        'status': status,
        'note': 'Concrete compliance uses the weighted all-photon downstream dry-air response and its supervisory uncertainty. Primary-history transmission/Wilson is diagnostic only. Actual simulated points only; no TVL, interpolation, or extrapolation.',
        'uses_wilson_for_selection': False,
        'uses_interpolation': False,
        'uses_tvl': False,
        'uses_extrapolation': False,
    })
    # Preserve legacy keys for old consumers, but make their diagnostic nature
    # explicit and never let them drive status/pass selection.
    output['final_validation_primary_photons'] = total_diag
    output['final_validation_transmitted_photons'] = transmitted_diag
    output['final_validation_observed_transmission'] = diag_observed
    output['final_validation_wilson95_low'] = float(diag_low)
    output['final_validation_wilson95_high'] = float(diag_high)
    return output

V73I2_PHOTON_REPORT_CONTRACT = 'weighted supervisory result authoritative; primary-history Wilson diagnostic only; immutable audit exported'
# <<< v7.3i.2 WEIGHTED PHOTON FINAL-VALIDATION REPORT CONTRACT <<<

# >>> v7.3k.1 ADAPTIVE DEEP-PENETRATION SOURCE-IMPORTANCE + STUDENT-T PATCH >>>
# The native C++ physics/scorer remains unchanged.  Production final-validation
# calls use source-spectrum importance sampling at the orchestration layer and
# combine the native eventwise sufficient statistics with exact likelihood
# weights.  Debug workflows retain the validated v7.3f path.
from v73k_variance_reduction import (
    V73K_FEATURE as V73K_PHOTON_FEATURE,
    should_use_v73k as _v73k_should_use,
    supervise_photon_variance_reduced as _v73k_supervise,
)

_v73k_previous_run_simulation = run_simulation
_v73k_native_photon_run = _v73f_native_photon_run
_v73k_previous_attach_final_concrete_validation = attach_final_concrete_validation


def _v73k_attach_audit(results, audit):
    if not isinstance(results, dict):
        raise TypeError('v7.3k photon supervisor returned a non-dictionary result')
    summary = results.setdefault('summary', {})
    if not isinstance(summary, dict):
        raise TypeError('v7.3k photon result summary is not a dictionary')
    for key in (
        'v73k_photon_validation_audit',
        'v73i2_photon_validation_audit',
        'v73h_photon_validation_audit',
        'v73f_photon_validation_audit',
        'v73e_photon_validation_audit',
        'v73d_photon_validation_audit',
    ):
        summary[key] = audit
    fc = audit.get('final_cumulative', {}) if isinstance(audit, dict) and isinstance(audit.get('final_cumulative'), dict) else {}
    for src, dst in (
        ('response', 'v73e_cumulative_photon_response'),
        ('sem', 'v73e_cumulative_photon_response_sem'),
        ('lower95', 'v73e_cumulative_photon_lower95'),
        ('upper95', 'v73e_cumulative_photon_upper95'),
        ('histories', 'v73e_cumulative_photon_histories'),
        ('crossings', 'v73e_cumulative_photon_crossings'),
    ):
        if src in fc:
            summary[dst] = fc.get(src)
    summary['v73k_variance_reduction_feature'] = V73K_PHOTON_FEATURE
    summary['v73k1_variance_reduction_feature'] = V73K_PHOTON_FEATURE
    summary['v73k1_pass_support_certified'] = fc.get('pass_support_certified')
    summary['v73k1_failure_support_certified'] = fc.get('failure_support_certified')
    summary['v73k2_failure_support_prerequisites_met'] = fc.get('failure_support_prerequisites_met')
    summary['v73k2_estimator_validation_internal'] = fc.get('estimator_validation_internal')
    summary['v73k2_deep_analog_validation'] = fc.get('deep_analog_validation')
    summary['v73k2_replicate_consistency'] = fc.get('replicate_consistency')
    summary['v73k2_group_support'] = fc.get('group_support')
    summary['v73k2_event_influence'] = fc.get('event_influence')
    summary['v73k2_tvl_diagnostic'] = fc.get('tvl_diagnostic')
    summary['v73k2_external_benchmark'] = fc.get('external_benchmark')
    summary['v73k2_monotonicity_diagnostic'] = fc.get('monotonicity_diagnostic')
    summary['v73k1_batch_count'] = fc.get('batch_count')
    summary['v73k1_nonzero_batch_count'] = fc.get('nonzero_batch_count')
    summary['v73k1_adaptive_proposal'] = fc.get('adaptive_proposal')
    summary['v73k1_energy_diagnostics'] = fc.get('energy_diagnostics')
    summary['v73k10d_wall_face_distance_m'] = getattr(audit,'get',lambda *_:None)('wall_face_distance_m') if isinstance(audit,dict) else None
    summary['v73k10d_target_is_thickness_dependent'] = bool(audit.get('target_is_thickness_dependent')) if isinstance(audit,dict) else False
    summary['v73k10d_distance_geometry_policy'] = 'd_wall(t)=wall_face+t; candidate-specific B_gamma,max(t)'
    summary['v73k_compliance_estimator_role'] = (
        'AUTHORITATIVE for production final validation: source-spectrum importance-sampled original-primary event moments; '
        'plot-bearing analog batch remains diagnostic/visualization only'
    )
    return results


def run_simulation(cfg):
    if bool(getattr(cfg, '_v73k_internal_photon', False)):
        return _v73k_native_photon_run(cfg)
    if not _v73k_should_use(cfg):
        # Debug 10k/100k retains the inherited lightweight supervisor.  When a
        # completed k.10d design exists, freeze its scalar compatibility target
        # at the dynamic observed candidate. If that legacy debug supervisor
        # marches thicker, the frozen candidate target is conservative because
        # B_gamma,max(t) only increases with thickness. Production v7.3k uses
        # the exact target separately at every thickness.
        design=getattr(cfg,'v6_photon_design',None)
        if bool(getattr(cfg,'_v73k10d_dynamic_wall_target',False)) and isinstance(design,dict):
            cand=design.get('first_observed_candidate_thickness_cm')
            try:
                cand=float(cand); dyn=_v724_target_for_thickness(cfg,cand,_v724_target(cfg))
            except Exception:
                dyn=np.nan
            if np.isfinite(dyn) and 0.0<dyn<=1.0:
                try:
                    cfg=_v724_deepcopy(cfg)
                    for _name in ('required_max_photon_transmission','photon_transmission_target','target_transmission','required_transmission','max_transmission'):
                        setattr(cfg,_name,float(dyn))
                    setattr(cfg,'_v73k10d_debug_target_thickness_cm',float(cand))
                except Exception:
                    pass
        return _v73k_previous_run_simulation(cfg)
    results, audit = _v73k_supervise(
        _v73k_native_photon_run,
        cfg,
        get_digitized_spectrum,
    )
    return _v73k_attach_audit(results, audit)


def attach_final_concrete_validation(requirement, results):
    output = _v73k_previous_attach_final_concrete_validation(requirement, results)
    audit = _v73i2_photon_audit(results)
    if not str(audit.get('version', '')).strip().lower().startswith('v7.3k'):
        return output
    state = str(audit.get('final_state') or 'STATISTICALLY_UNRESOLVED').upper()
    output['feature_build'] = V73K_PHOTON_FEATURE
    output['method'] = (
        'all-photon downstream dry-air E*(mu_en/rho) response; 25 cm coarse -> direct 5 cm refinement; '
        'production final validation uses one pilot native process then four independent frozen-proposal native production replicates; each native run samples the complete grouped q(E) internally with exact eventwise p/q and likelihood-weighted independent-primary '
        'native event moments, separate Student-t PASS/FAIL support, and a cumulative batch-means cross-check; targeted candidate validation; no TVL/interpolation/extrapolation'
    )
    output['v73k_variance_reduction'] = True
    output['v73k1_adaptive_variance_reduction'] = True
    output['v73k_source_importance_alpha'] = audit.get('importance_alpha')
    output['v73k_source_importance_beta'] = audit.get('importance_beta')
    output['v73k_total_supervisory_native_jobs'] = audit.get('total_supervisory_native_jobs')
    output['v73k_source_likelihood_ess'] = (
        audit.get('final_cumulative', {}).get('source_likelihood_ess')
        if isinstance(audit.get('final_cumulative'), dict) else None
    )
    if state == 'TARGET_BRACKET_UNRESOLVED_AT_PER_CANDIDATE_MAX':
        output['status'] = (
            'STATISTICALLY UNRESOLVED — targeted candidate reached the user-selected per-candidate history maximum; '
            'search stopped at that unresolved thickness rather than marching through thicker walls'
        )
    elif state == 'CONFIGURED_CONCRETE_RANGE_EXHAUSTED_RESOLVED_FAIL':
        output['status'] = 'CONFIGURED CONCRETE RANGE EXHAUSTED — final sampled thickness is a statistically resolved failure; no extrapolation'
    elif state == 'VR_ANALOG_CROSSCHECK_FAILED_CLOSED':
        output['status'] = 'FAIL-CLOSED — source-importance estimator did not reproduce the analog reference within the declared cross-check tolerance'
    elif state == 'VERIFIED_PASS_MINIMUM_COMPLIANT':
        output['status'] = 'VERIFIED PASS — source-importance primary-moment upper is within B_gamma,max and all required thinner sampled candidates are resolved'
    output['note'] = (
        'Concrete compliance uses the v7.3k.7 pilot-trained/frozen inline full-proposal source-importance likelihood-weighted original-primary all-photon downstream dry-air response after same-thickness deep analog estimator validation. Every coarse/fine scan point is a direct full-physical-spectrum analog run. '
        'Raw/split/track crossing counts are diagnostic only; the returned analog detailed batch supplies plots and does not override the VR decision. '
        'Actual simulated thicknesses only; no TVL, interpolation, extrapolation, room geometry, maze credit, or detector dimensions.'
    )
    return output

V73K_PHOTON_CONTRACT = (
    'production deep-penetration validation: full-support adaptive source-spectrum importance sampling inside one native run per replicate; exact eventwise likelihood-weighted native '
    'eventwise sufficient-statistic fold; original-primary statistical unit; PASS >=100 nonzero + <=10% relative SEM + >=16 cumulative independent batches; '
    'FAIL uses a separate Student-t lower-bound gate with >=8 batches; same-thickness deep analog/VR validation plus group/ESS/replicate/event-influence gates; targeted stop-on-unresolved search'
)
# <<< v7.3k DEEP-PENETRATION SOURCE-IMPORTANCE + TARGETED VALIDATION PATCH <<<


# >>> v7.3k.11 SINGLE-OWNER WALL SEARCH CORE ROUTING >>>
V73K11_CORE_WALL_SEARCH_FEATURE = "2026-08-19-v7.3k.11-single-owner-wall-search-per-candidate-validation"
_v73k11_previous_run_thickness_sweep = run_thickness_sweep


def _v73k11_stage_for_requested_grid(values):
    """Classify the single-owner exact grid without rejecting a clipped range end.

    The GUI intentionally includes the user's configured MAX thickness even when
    that endpoint is not exactly one full coarse/fine step above the preceding
    point.  Example: the e14 TVL-seeded 275--340 cm wall range requests
    [275, 300, 325, 340].  This is still one COARSE_25CM owner with a 15 cm
    terminal clip; it is not a competing EXACT_REQUESTED_GRID owner.

    Only the terminal interval may be clipped.  Arbitrary irregular interior
    spacing remains EXACT_REQUESTED_GRID and therefore continues to fail closed
    under the inherited v7.3k.11a report contract.
    """
    arr = np.asarray(values, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    arr = np.unique(np.round(arr, 10))
    arr.sort()
    if arr.size <= 1:
        return "FINE_5CM"

    diffs = np.diff(arr)
    tol = 1.0e-9

    def _regular_or_terminally_clipped(step):
        if diffs.size == 0:
            return False
        if np.any(diffs <= tol):
            return False
        if diffs.size == 1:
            return bool(diffs[0] <= step + tol)
        return bool(
            np.all(np.abs(diffs[:-1] - step) <= tol)
            and diffs[-1] <= step + tol
        )

    # Fine must be checked first for a short two-point interval (for example
    # [335, 340]); otherwise any <=25 cm pair would be mislabeled coarse.
    if _regular_or_terminally_clipped(5.0):
        return "FINE_5CM"
    if _regular_or_terminally_clipped(25.0):
        return "COARSE_25CM"
    return "EXACT_REQUESTED_GRID"


def run_thickness_sweep(cfg):
    """Execute exactly the grid requested by the GUI; never spawn another grid.

    v7.3k.11 makes the GUI/supervisory layer the sole wall-search orchestrator.
    The core remains a pure exact-grid Geant4 executor.  This wrapper deliberately
    disables the inherited C-25..C+15 internal refinement path.
    """
    requested = [float(x) for x in np.asarray(getattr(cfg, "sweep_thicknesses_cm", []), dtype=float).ravel()]
    if not requested:
        raise ValueError("v7.3k.11 exact wall sweep requires at least one requested thickness.")
    if any(not np.isfinite(x) for x in requested):
        raise ValueError("v7.3k.11 requested wall sweep contains a non-finite thickness.")

    stage = _v73k11_stage_for_requested_grid(requested)
    purpose = {
        "COARSE_25CM": "25 cm coarse screening — single-owner wall search",
        "FINE_5CM": "fresh 5 cm fine screening including bracket endpoints",
        "EXACT_REQUESTED_GRID": "explicit exact wall screening grid",
    }[stage]

    saved = {}
    for name in ("sweep_exact_requested_grid", "sweep_stage_override", "sweep_purpose_override"):
        saved[name] = (hasattr(cfg, name), getattr(cfg, name, None))
    try:
        cfg.sweep_exact_requested_grid = True
        cfg.sweep_stage_override = stage
        cfg.sweep_purpose_override = purpose
        print(
            f"v7.3k.11 WALL SEARCH EXECUTOR: stage={stage} | exact requested grid | "
            f"points={len(requested)} | STAGE_PHOTONS={int(getattr(cfg, 'sweep_num_photons_per_point', 0) or 0):,}",
            flush=True,
        )
        result = _v73k11_previous_run_thickness_sweep(cfg)
    finally:
        for name, (existed, value) in saved.items():
            if existed:
                setattr(cfg, name, value)
            else:
                try:
                    delattr(cfg, name)
                except Exception:
                    pass

    actual = np.asarray(result.get("thicknesses_cm", []), dtype=float).ravel() if isinstance(result, dict) else np.asarray([], dtype=float)
    req = np.unique(np.round(np.asarray(requested, dtype=float), 9))
    got = np.unique(np.round(actual, 9))
    if req.size != got.size or not np.allclose(req, got, rtol=0.0, atol=1.0e-8):
        raise RuntimeError(
            "v7.3k.11 fail-closed: core executed a thickness grid different from the exact GUI request; "
            f"requested={req.tolist()} actual={got.tolist()}"
        )
    if isinstance(result, dict):
        result["v73k11_wall_search_executor"] = V73K11_CORE_WALL_SEARCH_FEATURE
        result["v73k11_wall_search_stage"] = stage
        result["v73k11_exact_requested_grid"] = True
    return result
# <<< v7.3k.11 SINGLE-OWNER WALL SEARCH CORE ROUTING <<<

# >>> v7.3k.11a CORE LOAD-BOUNDARY ACTIVATION MARKER >>>
V73K11A_CORE_LOAD_FEATURE = "2026-08-19-v7.3k.11a-core-load-boundary-single-owner-hotfix"
# The GUI's load_core() compiles only the source prefix before
# ``t_total_start = time.perf_counter()``.  This marker and the executor
# attribute therefore MUST live before that sentinel together with the active
# v7.3k.11 run_thickness_sweep() wrapper.
run_thickness_sweep._v73k11_exact_grid_executor = True
# <<< v7.3k.11a CORE LOAD-BOUNDARY ACTIVATION MARKER <<<

# >>> v7.3k.12 10 MV DERIVED-SPECTRUM TRANSFORM + HARD PREFLIGHT >>>
V73K12_SPECTRUM_FEATURE = "2026-08-19-v7.3k.12-supervisory-minimum-proof-spectrum-preflight"
V73K12_10MV_ORIGIN = "project_derived_from_6MV_energy_axis_scaled_10_over_6_probability_conserving_rebin"
V73K12_10MV_TRANSFORM = "E10=E6*(10/6); probability mass linearly split onto the existing 0.1-MeV common-center grid; no invented tail"
V73K12_10MV_CLASSIFICATION = "DERIVED_PROJECT_SOURCE_NOT_MACHINE_COMMISSIONED"
_v73k12_previous_get_spectrum_provenance_note = get_spectrum_provenance_note
_v73k12_previous_get_spectrum_revision = get_spectrum_revision


def _v73k12_probability_conserving_rebin(source_energy, source_probability, target_energy, scale=1.0):
    src_e = np.asarray(source_energy, dtype=float).ravel()
    src_p = np.asarray(source_probability, dtype=float).ravel()
    dst_e = np.asarray(target_energy, dtype=float).ravel()
    if src_e.size == 0 or src_e.size != src_p.size or dst_e.size < 2:
        raise ValueError("v7.3k.12 spectrum rebin requires nonempty equal source arrays and >=2 target bins")
    if np.any(~np.isfinite(src_e)) or np.any(~np.isfinite(src_p)) or np.any(~np.isfinite(dst_e)):
        raise ValueError("v7.3k.12 spectrum rebin received non-finite values")
    if np.any(np.diff(dst_e) <= 0.0):
        raise ValueError("v7.3k.12 target energy grid must be strictly increasing")
    if not np.isfinite(float(scale)) or float(scale) <= 0.0:
        raise ValueError("v7.3k.12 energy-axis scale must be finite and positive")
    src_p = np.maximum(src_p, 0.0)
    total = float(np.sum(src_p))
    if total <= 0.0:
        raise ValueError("v7.3k.12 source spectrum has no positive probability mass")
    src_p = src_p / total
    scaled = src_e * float(scale)
    out = np.zeros(dst_e.shape, dtype=float)
    tol = 1.0e-10
    for energy, mass in zip(scaled, src_p):
        if mass <= 0.0:
            continue
        if energy < dst_e[0] - tol or energy > dst_e[-1] + tol:
            raise ValueError(
                f"v7.3k.12 scaled spectrum energy {energy:g} MeV lies outside target grid "
                f"{dst_e[0]:g}..{dst_e[-1]:g} MeV"
            )
        if energy <= dst_e[0] + tol:
            out[0] += mass
            continue
        if energy >= dst_e[-1] - tol:
            out[-1] += mass
            continue
        right = int(np.searchsorted(dst_e, energy, side="left"))
        if right < dst_e.size and abs(dst_e[right] - energy) <= tol:
            out[right] += mass
            continue
        left = right - 1
        span = float(dst_e[right] - dst_e[left])
        w_right = float((energy - dst_e[left]) / span)
        w_left = 1.0 - w_right
        out[left] += mass * w_left
        out[right] += mass * w_right
    out_sum = float(np.sum(out))
    if out_sum <= 0.0 or not np.isfinite(out_sum):
        raise ValueError("v7.3k.12 transformed spectrum has invalid total probability")
    out /= out_sum
    return out


def _v73k12_expected_10mv_probability():
    base = SPECTRUM_LIBRARY.get("6MV_40x40")
    if not isinstance(base, dict):
        raise RuntimeError("v7.3k.12 requires the existing 6MV_40x40 project source")
    return _v73k12_probability_conserving_rebin(
        base.get("energy_mev_center", []),
        base.get("probability_mass_bin", []),
        E_40x40_COMMON,
        scale=10.0 / 6.0,
    )


def _v73k12_install_10mv_derived_source():
    global P_10MV_40x40_common
    expected = _v73k12_expected_10mv_probability()
    P_10MV_40x40_common = expected.tolist()
    SPECTRUM_LIBRARY["10MV_40x40_interpolated"] = _make_spectrum(
        case="10MV_40x40_interpolated",
        beam_mv=10,
        field_cm="40x40",
        peak_mev=0.0,
        energy_mev_center=E_40x40_COMMON,
        probability_mass_bin=P_10MV_40x40_common,
        origin=V73K12_10MV_ORIGIN,
        interpolation_method="energy_axis_scaled_10_over_6_probability_conserving_rebin_from_6MV_40x40",
    )
    SPECTRUM_LIBRARY["10MV_40x40_interpolated"]["source_base_case"] = "6MV_40x40"
    SPECTRUM_LIBRARY["10MV_40x40_interpolated"]["source_validation_class"] = V73K12_10MV_CLASSIFICATION
    SPECTRUM_LIBRARY["10MV_40x40_interpolated"]["transformation_contract"] = V73K12_10MV_TRANSFORM


_v73k12_install_10mv_derived_source()


def get_spectrum_provenance_note(case: str) -> str:
    if str(case) == "10MV_40x40_interpolated":
        return (
            "Project-derived 10 MV source from the existing 6MV_40x40 project spectrum by explicit "
            "10/6 energy-axis scaling and probability-conserving rebinning onto the 0.1-MeV common grid; "
            "not asserted as a commissioned machine spectrum or an independent published 10 MV measurement"
        )
    return _v73k12_previous_get_spectrum_provenance_note(case)


def get_spectrum_revision(case: str) -> str:
    if str(case) == "10MV_40x40_interpolated":
        return V73K12_SPECTRUM_FEATURE
    return _v73k12_previous_get_spectrum_revision(case)


def v73k12_spectrum_preflight(cfg):
    case = get_spectrum_case(cfg)
    if case not in SPECTRUM_LIBRARY:
        raise RuntimeError(f"v7.3k.12 SPECTRUM PREFLIGHT FAILED: unknown case {case!r}")
    entry = SPECTRUM_LIBRARY[case]
    energy = np.asarray(entry.get("energy_mev_center", []), dtype=float).ravel()
    probability = np.asarray(entry.get("probability_mass_bin", []), dtype=float).ravel()
    errors = []
    if energy.size == 0 or energy.size != probability.size:
        errors.append("energy/probability arrays are empty or length-mismatched")
    if energy.size and (np.any(~np.isfinite(energy)) or np.any(np.diff(energy) <= 0.0)):
        errors.append("energy grid is not finite and strictly increasing")
    if probability.size and (np.any(~np.isfinite(probability)) or np.any(probability < 0.0)):
        errors.append("probability array contains non-finite or negative mass")
    raw_sum = float(np.sum(probability)) if probability.size else 0.0
    if not np.isfinite(raw_sum) or raw_sum <= 0.0:
        errors.append("probability mass is non-positive")
    if errors:
        raise RuntimeError("v7.3k.12 SPECTRUM PREFLIGHT FAILED: " + "; ".join(errors))

    spec = get_digitized_spectrum(cfg)
    p = np.asarray(spec["probability_mass_bin"], dtype=float)
    e = np.asarray(spec["energy_mev_center"], dtype=float)
    norm = float(np.sum(p))
    if abs(norm - 1.0) > 1.0e-12:
        raise RuntimeError(f"v7.3k.12 SPECTRUM PREFLIGHT FAILED: normalized probability sum={norm:.17g}")

    transform_consistency = "NOT APPLICABLE"
    base_mean = None
    expected_mean = None
    if case == "10MV_40x40_interpolated":
        expected = _v73k12_expected_10mv_probability()
        actual = np.asarray(SPECTRUM_LIBRARY[case]["probability_mass_bin"], dtype=float)
        if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0.0, atol=2.0e-15):
            raise RuntimeError(
                "v7.3k.12 SPECTRUM PREFLIGHT FAILED: 10 MV array does not equal the declared "
                "6 MV energy-axis-scaled probability-conserving transform"
            )
        base_cfg = deepcopy(cfg)
        base_cfg.spectrum_case = "6MV_40x40"
        base_spec = get_digitized_spectrum(base_cfg)
        base_mean = float(base_spec["mean_energy_mev"])
        expected_mean = base_mean * (10.0 / 6.0)
        actual_mean = float(spec["mean_energy_mev"])
        if abs(actual_mean - expected_mean) > max(1.0e-12, abs(expected_mean) * 1.0e-10):
            raise RuntimeError(
                f"v7.3k.12 SPECTRUM PREFLIGHT FAILED: transformed mean={actual_mean:.12g} MeV "
                f"but 10/6-scaled base mean={expected_mean:.12g} MeV"
            )
        transform_consistency = "PASS"

    return {
        "feature": V73K12_SPECTRUM_FEATURE,
        "case": case,
        "beam_mv": int(entry.get("beam_mv", 0) or 0),
        "origin": str(entry.get("origin", "unknown")),
        "classification": str(entry.get("source_validation_class", "PROJECT_SOURCE")),
        "transformation": str(entry.get("transformation_contract", entry.get("interpolation_method", "none"))),
        "base_case": str(entry.get("source_base_case", "not applicable")),
        "probability_sum": norm,
        "mean_energy_mev": float(spec["mean_energy_mev"]),
        "base_mean_energy_mev": base_mean,
        "expected_scaled_mean_energy_mev": expected_mean,
        "support_min_mev": float(e[0]),
        "support_max_mev": float(spec["max_energy_mev"]),
        "last_positive_center_mev": float(e[-1]),
        "array_sha256": get_spectrum_array_sha256(case),
        "transform_consistency": transform_consistency,
        "status": "PASS",
        "machine_commissioned": False if case == "10MV_40x40_interpolated" else None,
    }


V73K12_SPECTRUM_CONTRACT = (
    "10MV_40x40_interpolated is rebuilt from 6MV_40x40 by explicit E*10/6 scaling and "
    "probability-conserving rebinning; hard runtime preflight; no fabricated high-energy tail"
)
# <<< v7.3k.12 10 MV DERIVED-SPECTRUM TRANSFORM + HARD PREFLIGHT <<<

# >>> v7.3k.12e FINE/FULL-ANALOG RESULT RETENTION + VR PLOT PROVENANCE >>>
V73K12E_CORE_FEATURE = "2026-08-20-v7.3k.12e-production-vr-k12-integration-shared-fine-final"
_v73k12e_previous_run_thickness_sweep = run_thickness_sweep


def _v73k12e4_stage_from_metadata(point_cfg, parent_cfg, result=None):
    tokens = [
        getattr(point_cfg, "sweep_stage_override", ""),
        getattr(parent_cfg, "sweep_stage_override", ""),
        getattr(point_cfg, "transport_subrun_id", ""),
        getattr(point_cfg, "transport_subrun_purpose", ""),
    ]
    if isinstance(result, dict):
        summary = result.get("summary", {}) if isinstance(result.get("summary"), dict) else {}
        tokens.extend((
            summary.get("transport_subrun_id", ""),
            summary.get("transport_subrun_purpose", ""),
            summary.get("sweep_stage", ""),
        ))
    joined = " ".join(str(x or "") for x in tokens).upper().replace("-", "_")
    if "FINE_5CM" in joined:
        return "FINE_5CM"
    if "COARSE_25CM" in joined:
        return "COARSE_25CM"
    return ""


def run_thickness_sweep(cfg):
    """Retain full FINE_5CM analog results for later zero-cost plotting/VR validation reuse."""
    cache = getattr(cfg, "_v73k12e_full_spectrum_result_cache", None)
    if not isinstance(cache, dict):
        cache = {}
        setattr(cfg, "_v73k12e_full_spectrum_result_cache", cache)
    original_run_simulation = globals()["run_simulation"]
    parent_tracks = int(getattr(cfg, "n_visual_tracks", 0) or 0)
    parent_samples = int(getattr(cfg, "max_raw_exit_samples", 0) or 0)

    def capturing_run(point_cfg):
        stage = _v73k12e4_stage_from_metadata(point_cfg, cfg)
        if stage == "FINE_5CM" and hasattr(point_cfg, "n_visual_tracks"):
            # Bounded trajectory reservoir only; no additional histories and no
            # per-fine-point OBJ export. The selected result can be exported later.
            point_cfg.n_visual_tracks = parent_tracks
        if stage == "FINE_5CM" and hasattr(point_cfg, "max_raw_exit_samples") and parent_samples > 0:
            point_cfg.max_raw_exit_samples = parent_samples
        result = original_run_simulation(point_cfg)
        post_stage = _v73k12e4_stage_from_metadata(point_cfg, cfg, result)
        if post_stage:
            stage = post_stage
        if stage == "FINE_5CM" and isinstance(result, dict):
            summary = result.get("summary", {}) if isinstance(result.get("summary"), dict) else {}
            try:
                t = round(float(summary.get("wall_thickness_cm", getattr(point_cfg, "wall_thickness_cm"))), 9)
                n = int(summary.get("num_photons", getattr(point_cfg, "num_photons", 0)) or 0)
            except Exception:
                return result
            cache[t] = {
                "result": result,
                "histories": n,
                "source": "FINE_5CM_FULL_PHYSICAL_SPECTRUM_ANALOG_PRE_SENTINEL_CACHE",
                "thickness_cm": t,
                "subrun_id": str(summary.get("transport_subrun_id", getattr(point_cfg, "transport_subrun_id", ""))),
            }
        return result

    globals()["run_simulation"] = capturing_run
    try:
        return _v73k12e_previous_run_thickness_sweep(cfg)
    finally:
        globals()["run_simulation"] = original_run_simulation


# >>> v7.3k.12e4 CORE CACHE LOAD-BOUNDARY ACTIVATION >>>
V73K12E4_CORE_LOAD_FEATURE = "2026-08-25-v7.3k.12e4-core-cache-load-boundary-vr-plot-reuse-hotfix"
# load_core() compiles only the source before t_total_start.  The active
# k12e FINE-result retaining run_thickness_sweep wrapper must therefore be
# the last pre-sentinel wrapper.  Preserve the k11 exact-grid ownership flag
# because this wrapper delegates to that executor rather than replacing its
# search semantics.
run_thickness_sweep._v73k11_exact_grid_executor = True
run_thickness_sweep._v73k12e_full_result_cache_executor = True
run_thickness_sweep._v73k12e4_cache_stage_detection = "FINE_5CM subrun-aware pre/post transport"
# <<< v7.3k.12e4 CORE CACHE LOAD-BOUNDARY ACTIVATION <<<

# Production VR attaches its audit in this core layer. Synchronize the k12d
# console/report aliases to the truthful full-spectrum analog plot source.
if "_v73k_attach_audit" in globals():
    _v73k12e_previous_attach_audit = _v73k_attach_audit

    def _v73k_attach_audit(results, audit):
        out = _v73k12e_previous_attach_audit(results, audit)
        if isinstance(out, dict) and isinstance(audit, dict) and audit.get("v73k12e_feature") == V73K12E_CORE_FEATURE:
            summary = out.setdefault("summary", {})
            plot_n = int(audit.get("plot_aggregate_histories", 0) or 0)
            summary["v73k12d_plot_aggregate_histories"] = plot_n
            summary["v73k12d_plot_aggregate_batch_count"] = int(audit.get("plot_aggregate_batch_count", 0) or 0)
            summary["v73k12d_trajectory_sample_histories"] = int(audit.get("trajectory_sample_histories", 0) or 0)
            summary["v73k12d_plot_extra_native_histories"] = int(audit.get("plot_extra_native_histories", 0) or 0)
            summary["v73k12e_plot_dataset_role"] = (
                "VR SAME-THICKNESS FULL-PHYSICAL-SPECTRUM ANALOG REUSE — shared FINE/FINAL N; "
                "formal compliance estimator is separately likelihood-weighted"
            )
            summary["v73k12e_plot_source_kind"] = str(audit.get("plot_source_kind", "NOT AVAILABLE"))
            summary["v73k12e_formal_decision_histories"] = int(audit.get("formal_production_decision_histories", 0) or 0)
            summary["v73k12e_pilot_histories_training_only"] = int(audit.get("proposal_training_pilot_histories", 0) or 0)
        return out
# <<< v7.3k.12e FINE/FULL-ANALOG RESULT RETENTION + VR PLOT PROVENANCE <<<


# >>> v7.3k.12e9 CONSOLIDATED DEPOSITION / ABSORBED-DOSE DEPTH PLOT >>>
V73K12E9_CORE_FEATURE = "2026-08-27-v7.3k.12e9-consolidated-deposition-absorbed-dose-depth-plot"
V73K12E9_CONSOLIDATED_DEPTH_TITLE = "Shield Energy Deposition / Whole-Slice Absorbed Dose vs Depth"
_v73k12e9_previous_plot_results = plot_results


def _v73k12e9_depth_plot_title(figure):
    axes = [axis for axis in figure.axes if getattr(axis, "get_title", lambda: "")()]
    axis = axes[0] if axes else (figure.axes[0] if figure.axes else None)
    return (axis.get_title().strip() if axis is not None else ""), axis


def _v73k12e9_remove_redundant_depth_plot_files(cfg):
    try:
        folder = os.path.abspath(str(cfg.plot_dir))
        patterns = (
            "*_actual_deposited_energy_vs_depth.png",
            "*_dose_vs_depth.png",
            "*_equivalent_dose_proxy_vs_depth.png",
        )
        import glob as _v73k12e9_glob
        for pattern in patterns:
            for path in _v73k12e9_glob.glob(os.path.join(folder, pattern)):
                try:
                    os.remove(path)
                except OSError:
                    pass
    except Exception:
        pass


def _v73k12e9_make_consolidated_depth_plot(results, cfg):
    depth = np.asarray(results.get("depth_cm", []), dtype=float).ravel()
    edep = np.asarray(results.get("edep_per_primary_per_cm", []), dtype=float).ravel()
    edep_sem = np.asarray(results.get("edep_per_primary_per_cm_sigma", []), dtype=float).ravel()
    dose = np.asarray(results.get("dose_vs_depth_gy", []), dtype=float).ravel()

    n = min(depth.size, edep.size, edep_sem.size, dose.size)
    if n <= 0:
        return None
    depth = depth[:n]
    edep = edep[:n]
    edep_sem = edep_sem[:n]
    dose = dose[:n]
    finite = np.isfinite(depth) & np.isfinite(edep) & np.isfinite(edep_sem) & np.isfinite(dose)
    if not np.any(finite):
        return None

    fig, ax = plt.subplots()
    ax.plot(
        depth,
        edep,
        color=PLOT_RED,
        linewidth=PLOT_LINE_WIDTH,
        label="Deposited-energy tally",
    )
    ax.fill_between(
        depth,
        np.maximum(edep - edep_sem, 0.0),
        edep + edep_sem,
        color=PLOT_RED,
        alpha=PLOT_UNCERTAINTY_ALPHA,
        linewidth=0.0,
        label="±1σ SEM",
    )
    ax.set_xlabel("Depth x (cm)")
    ax.set_ylabel("Deposited energy (MeV per primary per cm)")
    ax.set_title(V73K12E9_CONSOLIDATED_DEPTH_TITLE)
    ax.legend()

    ratio_mask = finite & (np.abs(edep) > 0.0)
    dose_per_edep_density = math.nan
    if np.any(ratio_mask):
        ratios = dose[ratio_mask] / edep[ratio_mask]
        ratios = ratios[np.isfinite(ratios) & (ratios > 0.0)]
        if ratios.size:
            dose_per_edep_density = float(np.median(ratios))

    if math.isfinite(dose_per_edep_density) and dose_per_edep_density > 0.0:
        scale = dose_per_edep_density
        secondary = ax.secondary_yaxis(
            "right",
            functions=(lambda value, _s=scale: value * _s, lambda value, _s=scale: value / _s),
        )
        secondary.set_ylabel("Whole-slice absorbed dose per primary (Gy/primary)")
        secondary.tick_params(axis="y", colors=PLOT_RED)
        secondary.yaxis.label.set_color(PLOT_RED)
        try:
            secondary.spines["right"].set_color(PLOT_RED)
        except Exception:
            pass

    style_black_red_figure(fig)
    try:
        fig.tight_layout()
    except Exception:
        pass
    output_path = os.path.join(
        str(cfg.plot_dir),
        "v73k12e9_consolidated_deposition_absorbed_dose_vs_depth.png",
    )
    fig.savefig(
        output_path,
        dpi=cfg.plot_dpi,
        bbox_inches="tight",
        facecolor=PLOT_BLACK,
        edgecolor=PLOT_BLACK,
    )
    return fig


def plot_results(results, cfg, sweep_results=None):
    """v7.3k.12e9: replace three deterministic-equivalent depth panels with one plot.

    The underlying deposited-energy, absorbed-dose and wR=1 data arrays are retained
    unchanged for numerical/report compatibility.  Only redundant visual presentation
    is consolidated; no transport, scoring, normalization or statistical quantity is
    altered.
    """
    outer_show = plt.show

    def e9_show(*args, **kwargs):
        figure = plt.gcf()
        title, _axis = _v73k12e9_depth_plot_title(figure)
        title_lower = title.lower()
        if title_lower in {
            "shield-only deposited energy vs depth",
            "actual deposited energy vs depth",
            "whole-slice average absorbed dose per primary vs depth",
            "dose vs depth",
            "absorbed dose with unit weighting wr=1 (derived)",
            "equivalent dose proxy vs depth",
        }:
            plt.close(figure)
            return None
        return outer_show(*args, **kwargs)

    plt.show = e9_show
    try:
        _v73k12e9_previous_plot_results(results, cfg, sweep_results=sweep_results)
        _v73k12e9_remove_redundant_depth_plot_files(cfg)
        consolidated = _v73k12e9_make_consolidated_depth_plot(results, cfg)
        if consolidated is not None:
            plt.show()
    finally:
        plt.show = outer_show

    summary = results.setdefault("summary", {}) if isinstance(results, dict) else {}
    if isinstance(summary, dict):
        summary["v73k12e9_consolidated_depth_plot"] = True
        summary["v73k12e9_consolidated_depth_plot_title"] = V73K12E9_CONSOLIDATED_DEPTH_TITLE
        summary["v73k12e9_depth_plot_catalog_delta"] = -2
        summary["equivalent_dose_proxy_retired"] = True
        summary["equivalent_dose_proxy_definition"] = (
            "wR=1 deterministic alias of whole-slice absorbed dose; numerical array retained for compatibility, "
            "separate plot retired in v7.3k.12e9"
        )


plot_results._v73k12e9_consolidated_depth_plot = True
# <<< v7.3k.12e9 CONSOLIDATED DEPOSITION / ABSORBED-DOSE DEPTH PLOT <<<


# >>> v7.3k.12e10 SEPARATE COARSE / FINE-SCREEN / FINAL-VALIDATION HISTORY CONTROLS >>>
V73K12E10_HISTORY_CONTROL_FEATURE = "2026-08-27-v7.3k.12e10-separated-coarse-fine-final-history-controls"
_v73k12e10_previous_run_simulation = run_simulation


def _v73k12e10_positive_int(value, fallback=0):
    try:
        out = int(value)
    except Exception:
        out = int(fallback or 0)
    return out if out > 0 else int(fallback or 0)


def run_simulation(cfg):
    """Decorate the authoritative supervisory audit with separated history controls.

    Transport/scoring/statistical algorithms are unchanged.  FINE_5CM screening
    uses fine_screen_num_photons_per_point; the formal supervisor continues to
    read cfg.num_photons as its per-candidate cumulative maximum.
    """
    results = _v73k12e10_previous_run_simulation(cfg)
    if bool(getattr(cfg, '_v73k_internal_photon', False)) or not isinstance(results, dict):
        return results
    summary = results.setdefault('summary', {})
    if not isinstance(summary, dict):
        return results
    coarse_n = _v73k12e10_positive_int(getattr(cfg, 'coarse_num_photons_per_point', 0))
    fine_n = _v73k12e10_positive_int(
        getattr(cfg, 'fine_screen_num_photons_per_point', getattr(cfg, 'fine_num_photons_per_point', 0))
    )
    final_n = _v73k12e10_positive_int(
        getattr(cfg, 'final_validation_num_photons', getattr(cfg, 'final_num_photons', getattr(cfg, 'num_photons', 0)))
    )
    summary['v73k12e10_history_control_feature'] = V73K12E10_HISTORY_CONTROL_FEATURE
    summary['v73k12e10_coarse_photons'] = coarse_n
    summary['v73k12e10_fine_screen_photons'] = fine_n
    summary['v73k12e10_final_validation_photons'] = final_n
    summary['v73k12e10_history_control_mode'] = 'SEPARATE_COARSE_FINE_SCREEN_FINAL_VALIDATION'
    for key in ('v73i2_photon_validation_audit','v73h_photon_validation_audit','v73f_photon_validation_audit','v73e_photon_validation_audit','v73d_photon_validation_audit'):
        audit = summary.get(key)
        if isinstance(audit, dict):
            audit['v73k12e10_history_control_feature'] = V73K12E10_HISTORY_CONTROL_FEATURE
            audit['configured_coarse_photons'] = coarse_n
            audit['configured_fine_screen_photons'] = fine_n
            audit['configured_final_validation_photons'] = final_n
            # Old report layers still read this field.  From e10 onward it is a
            # compatibility alias for FINAL_VALIDATION_PHOTONS only.
            audit['configured_fine_final_photons'] = final_n
            audit['history_control_mode'] = 'SEPARATE_COARSE_FINE_SCREEN_FINAL_VALIDATION'
            audit['v73k12e10_formal_plot_n_independence_allowed'] = True
    if 'v73k12e_plot_dataset_role' in summary:
        summary['v73k12e_plot_dataset_role'] = (
            'VR SAME-THICKNESS FULL-PHYSICAL-SPECTRUM ANALOG REUSE — independent analog/FINE_SCREEN and '
            'formal FINAL_VALIDATION history budgets; compliance estimator is separately likelihood-weighted'
        )
    summary['v73k12e10_formal_plot_n_independence_allowed'] = True
    return results


run_simulation._v73k12e10_separate_history_controls = True
V73K12E10_HISTORY_CONTROL_CONTRACT = (
    'COARSE_PHOTONS sizes COARSE_25CM only; FINE_SCREEN_PHOTONS sizes every FINE_5CM screening point; '
    'FINAL_VALIDATION_PHOTONS is the independent per-candidate formal supervisory maximum; no preset ratio; '
    'transport physics/scoring/VR/statistical gates unchanged'
)
# <<< v7.3k.12e10 SEPARATE COARSE / FINE-SCREEN / FINAL-VALIDATION HISTORY CONTROLS <<<


# >>> v7.3k.12e11 EXACT CUMULATIVE BOUNDARY-HISTORY RECONCILIATION >>>
V73K12E11_BOUNDARY_RECON_FEATURE = "2026-08-27-v7.3k.12e11-exact-cumulative-boundary-history-reconciliation"
_v73k12e11_previous_run_simulation = run_simulation


def _v73k12e11_int_field(mapping, key, default=None):
    if not isinstance(mapping, dict):
        if default is None:
            raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: missing mapping for {key}")
        return int(default)
    value = mapping.get(key, default)
    if value is None:
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: missing {key}")
    try:
        out = int(value)
    except Exception as exc:
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: non-integer {key}={value!r}") from exc
    if out < 0:
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: negative {key}={out}")
    return out


def _v73k12e11_float_field(mapping, key):
    try:
        out = float(mapping.get(key))
    except Exception as exc:
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: invalid {key}") from exc
    if not math.isfinite(out):
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: non-finite {key}")
    return out


def _v73k12e11_capture_native_boundary(result):
    if not isinstance(result, dict):
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: native supervisory result is not a mapping")
    summary = result.get('summary', {})
    if not isinstance(summary, dict):
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: native supervisory summary is not a mapping")
    record = {
        'histories': _v73k12e11_int_field(summary, 'num_photons'),
        'thickness_cm': _v73k12e11_float_field(summary, 'wall_thickness_cm'),
        'transmitted': _v73k12e11_int_field(summary, 'transmitted_count'),
        'primary': _v73k12e11_int_field(summary, 'primary_photon_forward_crossing_count'),
        'secondary': _v73k12e11_int_field(summary, 'secondary_photon_forward_crossing_count'),
        'all_gamma': _v73k12e11_int_field(summary, 'all_photon_forward_crossing_count'),
        'histories_any': _v73k12e11_int_field(summary, 'histories_with_any_forward_gamma_crossing'),
        'boundary_errors': _v73k12e11_int_field(summary, 'boundary_classification_error_count'),
        'consistent': bool(summary.get('photon_boundary_tally_consistent', False)),
        'schema': str(summary.get('photon_boundary_tally_schema', 'not available')),
        'subrun_id': str(summary.get('transport_subrun_id', summary.get('subrun_id', 'not available'))),
    }
    if record['histories'] <= 0:
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: captured native batch has zero histories")
    if record['primary'] != record['transmitted']:
        raise RuntimeError(
            "v7.3k.12e11 boundary reconciliation: captured native batch violates primary/transmitted invariant "
            f"({record['primary']} != {record['transmitted']})"
        )
    if record['all_gamma'] != record['primary'] + record['secondary']:
        raise RuntimeError(
            "v7.3k.12e11 boundary reconciliation: captured native batch violates all=primary+secondary invariant "
            f"({record['all_gamma']} != {record['primary']}+{record['secondary']})"
        )
    if not (record['transmitted'] <= record['histories_any'] <= record['all_gamma']):
        raise RuntimeError(
            "v7.3k.12e11 boundary reconciliation: captured native batch unique-history count outside valid range "
            f"[{record['transmitted']}, {record['all_gamma']}]: {record['histories_any']}"
        )
    if record['boundary_errors'] != 0 or not record['consistent']:
        raise RuntimeError(
            "v7.3k.12e11 boundary reconciliation: captured native batch failed native boundary consistency "
            f"errors={record['boundary_errors']}, consistent={record['consistent']}"
        )
    return record


def _v73k12e11_latest_analog_audit(summary):
    if not isinstance(summary, dict):
        return None
    # Production VR publishes v73k_photon_validation_audit.  e11 intentionally
    # leaves that estimator/analog split untouched; this repair is for the
    # inherited lightweight exact-additive analog supervisor only.
    if isinstance(summary.get('v73k_photon_validation_audit'), dict):
        return None
    for key in ('v73i2_photon_validation_audit', 'v73h_photon_validation_audit', 'v73f_photon_validation_audit'):
        value = summary.get(key)
        if isinstance(value, dict):
            return value
    return None


def _v73k12e11_attempt_batch_sizes(audit, final_thickness_cm, final_histories):
    attempts = audit.get('attempts', []) if isinstance(audit, dict) else []
    if not isinstance(attempts, list):
        return []
    out = []
    for attempt in attempts:
        if not isinstance(attempt, dict):
            continue
        try:
            t = float(attempt.get('thickness_cm'))
            n = int(attempt.get('batch_histories', 0) or 0)
        except Exception:
            continue
        if math.isfinite(t) and abs(t - final_thickness_cm) <= 1.0e-7 and n > 0:
            out.append(n)
    if out and sum(out) == final_histories:
        return out
    return []


def _v73k12e11_select_final_batches(captures, audit):
    fc = audit.get('final_cumulative', {}) if isinstance(audit.get('final_cumulative'), dict) else {}
    final_histories = _v73k12e11_int_field(fc, 'histories')
    try:
        final_thickness = float(audit.get('final_tested_thickness_cm', audit.get('decision_thickness_cm')))
    except Exception as exc:
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: final tested thickness unavailable") from exc
    if not math.isfinite(final_thickness):
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: final tested thickness is non-finite")

    same = [x for x in captures if abs(float(x['thickness_cm']) - final_thickness) <= 1.0e-7]
    if not same:
        raise RuntimeError(
            f"v7.3k.12e11 boundary reconciliation: no captured native supervisory batches at final thickness {final_thickness:g} cm"
        )

    expected = _v73k12e11_attempt_batch_sizes(audit, final_thickness, final_histories)
    if expected:
        if len(same) < len(expected):
            raise RuntimeError(
                "v7.3k.12e11 boundary reconciliation: fewer captured final-thickness batches than audit attempts"
            )
        selected = same[-len(expected):]
        actual = [int(x['histories']) for x in selected]
        if actual != expected:
            raise RuntimeError(
                "v7.3k.12e11 boundary reconciliation: captured final-thickness batch sizes do not match audit attempts "
                f"actual={actual}, expected={expected}"
            )
    else:
        # Fail closed rather than silently taking an arbitrary subset.  A suffix
        # is allowed only when it closes exactly to final_cumulative.histories.
        selected = []
        running = 0
        for record in reversed(same):
            selected.append(record)
            running += int(record['histories'])
            if running == final_histories:
                selected.reverse()
                break
            if running > final_histories:
                selected = []
                break
        if not selected:
            raise RuntimeError(
                "v7.3k.12e11 boundary reconciliation: captured final-thickness histories cannot be matched exactly "
                f"to final cumulative N={final_histories}"
            )

    selected_histories = sum(int(x['histories']) for x in selected)
    if selected_histories != final_histories:
        raise RuntimeError(
            f"v7.3k.12e11 boundary reconciliation: selected batch N={selected_histories} != final cumulative N={final_histories}"
        )
    return selected, final_thickness, final_histories


def _v73k12e11_existing_nonnegative(summary, key):
    if key not in summary or summary.get(key) is None:
        return None
    try:
        value = int(summary.get(key))
    except Exception as exc:
        raise RuntimeError(f"v7.3k.12e11 boundary reconciliation: aggregate {key} is not an integer") from exc
    return value if value >= 0 else None


def _v73k12e11_reconcile_analog_cumulative(results, captures):
    if not isinstance(results, dict):
        return results
    summary = results.get('summary', {})
    if not isinstance(summary, dict):
        return results
    audit = _v73k12e11_latest_analog_audit(summary)
    if audit is None:
        summary['v73k12e11_boundary_reconciliation_feature'] = V73K12E11_BOUNDARY_RECON_FEATURE
        summary['v73k12e11_boundary_reconciliation_status'] = (
            'NOT APPLICABLE — production-VR or non-cumulative lightweight supervisor path'
        )
        return results
    fc = audit.get('final_cumulative', {}) if isinstance(audit.get('final_cumulative'), dict) else {}
    if _v73k12e11_int_field(fc, 'histories', 0) <= 0:
        summary['v73k12e11_boundary_reconciliation_feature'] = V73K12E11_BOUNDARY_RECON_FEATURE
        summary['v73k12e11_boundary_reconciliation_status'] = 'NOT APPLICABLE — no cumulative supervisory histories'
        return results

    selected, final_thickness, final_histories = _v73k12e11_select_final_batches(captures, audit)
    sums = {
        'transmitted': sum(x['transmitted'] for x in selected),
        'primary': sum(x['primary'] for x in selected),
        'secondary': sum(x['secondary'] for x in selected),
        'all_gamma': sum(x['all_gamma'] for x in selected),
        'histories_any': sum(x['histories_any'] for x in selected),
        'boundary_errors': sum(x['boundary_errors'] for x in selected),
    }
    if sums['primary'] != sums['transmitted']:
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: cumulative primary/transmitted invariant failed")
    if sums['all_gamma'] != sums['primary'] + sums['secondary']:
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: cumulative all=primary+secondary invariant failed")
    if not (sums['transmitted'] <= sums['histories_any'] <= sums['all_gamma']):
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: cumulative unique-history count outside valid range")
    if sums['boundary_errors'] != 0 or not all(bool(x['consistent']) for x in selected):
        raise RuntimeError("v7.3k.12e11 boundary reconciliation: cumulative native boundary consistency failed")

    aggregate_n = _v73k12e11_existing_nonnegative(summary, 'num_photons')
    if aggregate_n is not None and aggregate_n != final_histories:
        raise RuntimeError(
            f"v7.3k.12e11 boundary reconciliation: aggregate num_photons={aggregate_n} != final cumulative N={final_histories}"
        )
    for key, expected in (
        ('transmitted_count', sums['transmitted']),
        ('primary_photon_forward_crossing_count', sums['primary']),
        ('secondary_photon_forward_crossing_count', sums['secondary']),
        ('all_photon_forward_crossing_count', sums['all_gamma']),
    ):
        existing = _v73k12e11_existing_nonnegative(summary, key)
        if existing is not None and existing != expected:
            raise RuntimeError(
                f"v7.3k.12e11 boundary reconciliation: aggregate {key}={existing} disagrees with exact batch sum {expected}"
            )
    existing_errors = _v73k12e11_existing_nonnegative(summary, 'boundary_classification_error_count')
    if existing_errors not in (None, 0):
        raise RuntimeError(
            f"v7.3k.12e11 boundary reconciliation: aggregate boundary error count is nonzero ({existing_errors})"
        )

    # Only after every independent invariant above passes do we synchronize the
    # cumulative boundary bookkeeping.  This repairs stale aggregation metadata;
    # it does not alter transport, weighted response moments, or statistical gates.
    synchronized = {
        'num_photons': final_histories,
        'transmitted_count': sums['transmitted'],
        'primary_photon_forward_crossing_count': sums['primary'],
        'secondary_photon_forward_crossing_count': sums['secondary'],
        'all_photon_forward_crossing_count': sums['all_gamma'],
        'histories_with_any_forward_gamma_crossing': sums['histories_any'],
        'boundary_classification_error_count': 0,
        'photon_boundary_tally_consistent': True,
    }
    for key, value in synchronized.items():
        summary[key] = value
        results[key] = value
    summary['transmitted_fraction'] = float(sums['transmitted']) / float(final_histories)
    summary['v73k12e11_boundary_reconciliation_feature'] = V73K12E11_BOUNDARY_RECON_FEATURE
    summary['v73k12e11_boundary_reconciliation_status'] = (
        f"PASS — exact additive boundary-history sum across {len(selected)} independent supervisory native batches "
        f"at {final_thickness:g} cm; N={final_histories:,}; histories-with-any-forward-gamma={sums['histories_any']:,}"
    )
    summary['v73k12e11_boundary_reconciliation_batches'] = len(selected)
    summary['v73k12e11_boundary_reconciliation_histories'] = final_histories
    summary['v73k12e11_boundary_reconciliation_histories_any'] = sums['histories_any']
    summary['v73k12e11_boundary_reconciliation_primary'] = sums['primary']
    summary['v73k12e11_boundary_reconciliation_secondary'] = sums['secondary']
    summary['v73k12e11_boundary_reconciliation_all_gamma'] = sums['all_gamma']
    audit['v73k12e11_boundary_reconciliation_feature'] = V73K12E11_BOUNDARY_RECON_FEATURE
    audit['v73k12e11_boundary_reconciliation_status'] = summary['v73k12e11_boundary_reconciliation_status']
    audit['v73k12e11_boundary_reconciliation_batches'] = len(selected)
    audit['v73k12e11_boundary_reconciliation_histories_any'] = sums['histories_any']
    return results


def run_simulation(cfg):
    # Production VR uses _v73k_native_photon_run and is intentionally unaffected.
    # The inherited debug/lightweight supervisor calls the global
    # _v73f_native_photon_run; intercept those exact native batch results only
    # long enough to retain additive boundary-history bookkeeping.
    if bool(getattr(cfg, '_v73k_internal_photon', False)):
        return _v73k12e11_previous_run_simulation(cfg)
    original_native = globals().get('_v73f_native_photon_run')
    captures = []
    if not callable(original_native):
        return _v73k12e11_previous_run_simulation(cfg)

    def capturing_native(point_cfg):
        native_result = original_native(point_cfg)
        captures.append(_v73k12e11_capture_native_boundary(native_result))
        return native_result

    globals()['_v73f_native_photon_run'] = capturing_native
    try:
        results = _v73k12e11_previous_run_simulation(cfg)
    finally:
        globals()['_v73f_native_photon_run'] = original_native
    return _v73k12e11_reconcile_analog_cumulative(results, captures)


run_simulation._v73k12e11_exact_boundary_reconciliation = True

# >>> v7.3k.12e12 FINAL WRAPPER PREFLIGHT-MARKER CHAIN HOTFIX >>>
V73K12E12_PREFLIGHT_CHAIN_FEATURE = "2026-08-27-v7.3k.12e12-final-wrapper-preflight-marker-chain-hotfix"
# e11 wraps the e10 run_simulation callable.  Preserve the inherited e10
# activation marker on the final callable inspected by photon_wall_software.pyw.
# This changes no transport, search, VR, statistics, or result data.
run_simulation._v73k12e10_separate_history_controls = True
run_simulation._v73k12e11_exact_boundary_reconciliation = True
run_simulation._v73k12e12_preflight_marker_chain = True
# <<< v7.3k.12e12 FINAL WRAPPER PREFLIGHT-MARKER CHAIN HOTFIX <<<

V73K12E11_BOUNDARY_RECON_CONTRACT = (
    'lightweight multi-batch supervisory boundary metadata are reconstructed only from exact additive native batch tallies; '
    'primary=transmitted, all=primary+secondary, unique-history range, zero boundary errors, batch-history closure, and aggregate '
    'crossing totals all fail closed before synchronization; production VR/transport/scoring/statistical gates unchanged'
)
# <<< v7.3k.12e11 EXACT CUMULATIVE BOUNDARY-HISTORY RECONCILIATION <<<

# >>> v7.3k.12e13 PRODUCTION-VR INDEPENDENT-N LOAD MARKER >>>
V73K12E13_PRODUCTION_VR_FEATURE = "2026-08-27-v7.3k.12e13-separated-fine-final-production-vr-cache-contract"
_v73k12e13_loaded_supervisor = globals().get("_v73k_supervise")
V73K12E13_PRODUCTION_VR_RUNTIME_READY = bool(
    callable(_v73k12e13_loaded_supervisor)
    and getattr(_v73k12e13_loaded_supervisor, "_v73k12e13_independent_fine_final_vr", False)
    and getattr(_v73k12e13_loaded_supervisor, "_v73k12e13_feature", "") == V73K12E13_PRODUCTION_VR_FEATURE
)
# This is a load/preflight marker only. The production-VR implementation lives
# in v73k_variance_reduction.py and is imported as _v73k_supervise above.
# <<< v7.3k.12e13 PRODUCTION-VR INDEPENDENT-N LOAD MARKER <<<



# >>> v7.3k.12e15 USER HISTORY-BUDGET PREFLIGHT >>>
V73K12E15_HISTORY_BUDGET_FEATURE = "2026-09-02-v7.3k.12e15-history-budget-floor-and-vr-pilot-cap"
V73K12E15_MIN_HISTORY = 1_000_000
_v73k12e15_previous_run_thickness_sweep = run_thickness_sweep


def run_thickness_sweep(cfg):
    """Reject sub-1M user coarse/fine/final controls before transport starts."""
    def _positive(*names):
        for name in names:
            try:
                value = int(getattr(cfg, name, 0) or 0)
            except Exception:
                value = 0
            if value > 0:
                return value
        return 0

    coarse_n = _positive("coarse_num_photons_per_point", "sweep_num_photons_per_point")
    fine_n = _positive("fine_screen_num_photons_per_point", "fine_num_photons_per_point", "sweep_num_photons_per_point")
    final_n = _positive("final_validation_num_photons", "final_num_photons", "num_photons")

    if coarse_n < V73K12E15_MIN_HISTORY:
        raise RuntimeError(
            f"v7.3k.12e15 COARSE_PHOTONS={coarse_n:,} is below the project minimum of {V73K12E15_MIN_HISTORY:,}"
        )
    if fine_n < V73K12E15_MIN_HISTORY:
        raise RuntimeError(
            f"v7.3k.12e15 FINE_SCREEN_PHOTONS={fine_n:,} is below the project minimum of {V73K12E15_MIN_HISTORY:,}"
        )
    if final_n < V73K12E15_MIN_HISTORY:
        raise RuntimeError(
            f"v7.3k.12e15 FINAL_VALIDATION_PHOTONS={final_n:,} is below the project minimum of {V73K12E15_MIN_HISTORY:,}"
        )
    if final_n < fine_n:
        raise RuntimeError(
            f"v7.3k.12e15 FINAL_VALIDATION_PHOTONS={final_n:,} must be >= FINE_SCREEN_PHOTONS={fine_n:,}"
        )

    print("\n================ PHOTON HISTORY PLAN ================", flush=True)
    print(f"COARSE: {coarse_n:,} histories per coarse point", flush=True)
    print(f"FINE:   {fine_n:,} histories per fine point", flush=True)
    print(f"FINAL:  {final_n:,} formal histories at the decision candidate", flush=True)
    print(
        f"Production-VR pilot cap: {min(fine_n, final_n):,} histories; training-only and explicitly additional.",
        flush=True,
    )
    print("=====================================================\n", flush=True)
    return _v73k12e15_previous_run_thickness_sweep(cfg)


# Preserve executor/cache markers used by the downstream preflight chain.
for _name, _value in getattr(_v73k12e15_previous_run_thickness_sweep, "__dict__", {}).items():
    setattr(run_thickness_sweep, _name, _value)
run_thickness_sweep._v73k12e15_history_budget = True
run_thickness_sweep._v73k12e15_feature = V73K12E15_HISTORY_BUDGET_FEATURE
# <<< v7.3k.12e15 USER HISTORY-BUDGET PREFLIGHT <<<


# >>> v7.3k.12e16 SEQUENTIAL FINAL / SAFE-DESIGN CORE INTEGRATION >>>
V73K12E16_CORE_FEATURE = "2026-09-03-v7.3k.12e16-sequential-final-safe-design-search"
V73K12E16_CORE_MIN_HISTORY = 1_000_000

# The GUI worker loads photon_wall_simulation_core.py only through the
# t_total_start boundary.  This block is deliberately installed before that
# sentinel so the active workflow sees the e16 supervisor and report semantics.
if not (
    callable(globals().get("_v73k_supervise"))
    and bool(getattr(_v73k_supervise, "_v73k12e16_sequential_safe_design", False))
    and getattr(_v73k_supervise, "_v73k12e16_feature", "") == V73K12E16_CORE_FEATURE
):
    raise RuntimeError(
        "v7.3k.12e16 core preflight failed: the sequential safe-design production-VR supervisor is not active; "
        "restart the application after installing the complete e16 patch"
    )

V73K12E16_PRODUCTION_VR_RUNTIME_READY = True

# Replace the e15 history-plan wrapper rather than stacking contradictory
# console wording on top of it.  The e15 1M floor remains fully enforced here;
# the underlying physical coarse/fine executor is the exact callable saved by
# e15 before its wrapper was installed.
_v73k12e16_previous_history_wrapper = run_thickness_sweep
_v73k12e16_physical_sweep_executor = globals().get(
    "_v73k12e15_previous_run_thickness_sweep", _v73k12e16_previous_history_wrapper
)


def run_thickness_sweep(cfg):
    def _positive(*names):
        for name in names:
            try:
                value = int(getattr(cfg, name, 0) or 0)
            except Exception:
                value = 0
            if value > 0:
                return value
        return 0

    coarse_n = _positive("coarse_num_photons_per_point", "sweep_num_photons_per_point")
    fine_n = _positive(
        "fine_screen_num_photons_per_point", "fine_num_photons_per_point", "sweep_num_photons_per_point"
    )
    final_n = _positive("final_validation_num_photons", "final_num_photons", "num_photons")
    for label, value in (("COARSE", coarse_n), ("FINE_SCREEN", fine_n), ("FINAL_VALIDATION", final_n)):
        if value < V73K12E16_CORE_MIN_HISTORY:
            raise RuntimeError(
                f"v7.3k.12e16 {label}_PHOTONS={value:,} is below the project minimum "
                f"of {V73K12E16_CORE_MIN_HISTORY:,}"
            )
    if final_n < fine_n:
        raise RuntimeError(
            f"v7.3k.12e16 FINAL_VALIDATION_PHOTONS={final_n:,} must be >= FINE_SCREEN_PHOTONS={fine_n:,}"
        )

    requested_pilot_cap = _positive("v73k12e16_pilot_max_histories")
    pilot_cap = requested_pilot_cap or min(final_n, 5_000_000)
    pilot_cap = max(V73K12E16_CORE_MIN_HISTORY, min(final_n, int(pilot_cap)))
    requested_analog_cap = _positive("v73k12e16_deep_analog_max_histories")
    analog_cap = requested_analog_cap or min(final_n, 10_000_000)
    analog_cap = max(V73K12E16_CORE_MIN_HISTORY, min(final_n, int(analog_cap)))

    print("\n================ PHOTON HISTORY PLAN — v7.3k.12e16 ================", flush=True)
    print(f"COARSE: {coarse_n:,} histories per coarse point", flush=True)
    print(f"FINE:   {fine_n:,} histories per 5-cm screening point", flush=True)
    print(f"FINAL:  sequential formal stages from 1,000,000 up to {final_n:,} histories per candidate", flush=True)
    print(f"Pilot:  starts at 1,000,000; adaptive training-only cap {pilot_cap:,}", flush=True)
    print(f"Deep analog: bounded/sequential diagnostic cap {analog_cap:,}; never defaults to FINAL-sized transport", flush=True)
    print("Resolved FINAL stages stop early; unresolved candidates may advance by an actual 5-cm Geant4 step.", flush=True)
    print("Full Geant4 initialization/physics/process/transport/result output remains enabled.", flush=True)
    print("===================================================================\n", flush=True)
    return _v73k12e16_physical_sweep_executor(cfg)


for _name, _value in getattr(_v73k12e16_previous_history_wrapper, "__dict__", {}).items():
    setattr(run_thickness_sweep, _name, _value)
run_thickness_sweep._v73k12e15_history_budget = True
run_thickness_sweep._v73k12e15_feature = V73K12E15_HISTORY_BUDGET_FEATURE
run_thickness_sweep._v73k12e16_sequential_safe_design = True
run_thickness_sweep._v73k12e16_feature = V73K12E16_CORE_FEATURE


# Persist the e16 decision distinction in the result summary.  The final
# cumulative formal estimator remains authoritative; the full-spectrum analog
# result returned by the supervisor remains plotting/diagnostic provenance.
_v73k12e16_previous_attach_audit = _v73k_attach_audit


def _v73k_attach_audit(results, audit):
    out = _v73k12e16_previous_attach_audit(results, audit)
    if not isinstance(out, dict) or not isinstance(audit, dict):
        return out
    if audit.get("v73k12e16_feature") != V73K12E16_CORE_FEATURE:
        return out
    summary = out.setdefault("summary", {})
    fc = audit.get("final_cumulative", {}) if isinstance(audit.get("final_cumulative"), dict) else {}
    summary["v73k12e16_feature"] = V73K12E16_CORE_FEATURE
    summary["v73k12e16_sequential_final_validation"] = True
    summary["v73k12e16_safe_design_search"] = True
    summary["v73k12e16_accepted_verified_design_thickness_cm"] = audit.get("accepted_verified_design_thickness_cm")
    summary["v73k12e16_minimum_verified_compliant_thickness_cm"] = audit.get("minimum_verified_compliant_thickness_cm")
    summary["v73k12e16_exact_minimum_established"] = bool(audit.get("exact_minimum_established"))
    summary["v73k12e16_unresolved_thinner_candidate_thicknesses_cm"] = list(
        audit.get("unresolved_thinner_candidate_thicknesses_cm", []) or []
    )
    summary["v73k12e16_configured_final_validation_max_histories"] = int(
        audit.get("configured_final_validation_photons", 0) or 0
    )
    summary["v73k12e16_actual_final_formal_histories"] = int(
        audit.get("formal_production_decision_histories", fc.get("histories", 0)) or 0
    )
    summary["v73k12e16_numeric_lower95_diagnostic"] = fc.get("numeric_lower95_diagnostic")
    summary["v73k12e16_numeric_upper95_diagnostic"] = fc.get("numeric_upper95_diagnostic")
    summary["v73k12e16_certified_lower95"] = fc.get("lower95")
    summary["v73k12e16_certified_upper95"] = fc.get("upper95")
    summary["v73k12e16_bound_source"] = "AUTHORITATIVE_FINAL_CUMULATIVE_ONLY"
    summary["v73k12e16_resume_cache_enabled"] = bool(audit.get("v73k12e16_resume_cache_enabled", True))
    summary["v73k12e16_resume_cache_used"] = bool(audit.get("v73k12e16_resume_cache_used", False))
    summary["v73k12e16_resume_checkpoint_histories"] = int(audit.get("v73k12e16_resume_checkpoint_histories", 0) or 0)
    summary["v73k12e16_formal_histories_executed_this_run"] = int(audit.get("v73k12e16_formal_histories_executed_this_run", 0) or 0)
    summary["v73k12e16_resume_checkpoint_integrity"] = str(audit.get("v73k12e16_resume_checkpoint_integrity", "NOT_USED_FRESH_FORMAL_ESTIMATOR"))
    return out


_v73k12e16_previous_attach_final_concrete_validation = attach_final_concrete_validation


def attach_final_concrete_validation(requirement, results):
    output = _v73k12e16_previous_attach_final_concrete_validation(requirement, results)
    audit = _v73i2_photon_audit(results)
    if not isinstance(audit, dict) or audit.get("v73k12e16_feature") != V73K12E16_CORE_FEATURE:
        return output

    accepted = _v73i2_finite(audit.get("accepted_verified_design_thickness_cm"))
    minimum = _v73i2_finite(audit.get("minimum_verified_compliant_thickness_cm"))
    unresolved = [
        float(x) for x in (audit.get("unresolved_thinner_candidate_thicknesses_cm", []) or [])
        if _v73i2_finite(x) is not None
    ]
    final_t = _v73i2_finite(audit.get("final_tested_thickness_cm"))
    fc = audit.get("final_cumulative", {}) if isinstance(audit.get("final_cumulative"), dict) else {}
    configured_max = int(audit.get("configured_final_validation_photons", 0) or 0)
    actual_final_n = int(audit.get("formal_production_decision_histories", fc.get("histories", 0)) or 0)

    output["feature_build"] = V73K12E16_CORE_FEATURE
    output["v73k12e16_sequential_final_validation"] = True
    output["v73k12e16_safe_design_search"] = True
    output["accepted_verified_design_thickness_cm"] = accepted
    output["minimum_verified_compliant_thickness_cm"] = minimum
    output["exact_minimum_established"] = bool(minimum is not None)
    output["unresolved_thinner_candidate_thicknesses_cm"] = sorted(set(unresolved))
    output["final_tested_thickness_cm"] = final_t
    output["configured_final_validation_max_histories"] = configured_max
    output["actual_final_formal_histories"] = actual_final_n
    output["resume_cache_used"] = bool(audit.get("v73k12e16_resume_cache_used", False))
    output["resume_checkpoint_histories"] = int(audit.get("v73k12e16_resume_checkpoint_histories", 0) or 0)
    output["formal_histories_executed_this_run"] = int(audit.get("v73k12e16_formal_histories_executed_this_run", actual_final_n) or 0)
    output["resume_checkpoint_integrity"] = str(audit.get("v73k12e16_resume_checkpoint_integrity", "NOT_USED_FRESH_FORMAL_ESTIMATOR"))
    output["formal_bound_source"] = "AUTHORITATIVE_FINAL_CUMULATIVE_ONLY"
    output["analog_bound_fallback_used"] = False
    output["uses_wilson_for_selection"] = False
    output["uses_interpolation"] = False
    output["uses_tvl"] = False
    output["uses_extrapolation"] = False

    if accepted is not None:
        # An independently VERIFIED PASS is a valid conservative design
        # thickness.  It is the exact minimum only when every relevant thinner
        # candidate has been explicitly resolved.
        output["required_thickness_95_cm"] = float(accepted)
        output["verified_passing_thickness_cm"] = float(accepted)
        output["final_validation_meets_target_95"] = True
        if minimum is not None:
            output["status"] = (
                f"VERIFIED MINIMUM COMPLIANT — {float(accepted):.3f} cm; every relevant thinner candidate is resolved"
            )
        else:
            unresolved_text = ", ".join(f"{x:g}" for x in sorted(set(unresolved))) or "one or more thinner points"
            output["status"] = (
                f"VERIFIED COMPLIANT DESIGN THICKNESS — {float(accepted):.3f} cm; exact minimum not established because "
                f"thinner candidate(s) remain statistically unresolved: {unresolved_text} cm"
            )
    else:
        output["required_thickness_95_cm"] = None
        output["verified_passing_thickness_cm"] = None
        output["final_validation_meets_target_95"] = False
        output["status"] = (
            "REQUIRED CONCRETE THICKNESS NOT ESTABLISHED — no VERIFIED PASS was found within the configured concrete range"
        )

    output["method"] = (
        "all-photon downstream dry-air E*(mu_en/rho) response; 25-cm coarse -> direct 5-cm screening; "
        "production validation uses an adaptive training-only pilot followed by cumulative independent four-replicate "
        "1-2-5 formal stages up to FINAL_VALIDATION_PHOTONS; finite-look one-sided alpha spending controls within-candidate "
        "repeated-look family-wise alpha <=0.05, with unchanged strict Student-t PASS/FAIL support, ESS, "
        "group, replicate and event-influence gates; bounded same-thickness analog consistency diagnostics where only a "
        "support-certified inconsistency vetoes the formal decision; an unresolved candidate remains unresolved but does not "
        "prevent testing the next actual 5-cm concrete thickness; compatible exact sufficient-statistic checkpoints can be resumed"
    )
    output["note"] = (
        "Concrete compliance uses only the authoritative likelihood-weighted formal estimator at the actual tested wall "
        "thickness. Full-spectrum analog data are diagnostic/plot provenance and never replace formal confidence bounds. "
        "A VERIFIED PASS may be published as a conservative design thickness even when a thinner point remains unresolved; "
        "an underpowered deep analog diagnostic does not manufacture or veto PASS, while a supported inconsistency fails closed. "
        "The exact minimum is published only after all relevant thinner points are resolved. Actual Geant4 thicknesses only; "
        "no TVL/interpolation/extrapolation, room geometry, maze credit, or detector dimensions are used."
    )
    return output


V73K12E16_CORE_CONTRACT = (
    "pre-t_total_start activation; >=1M coarse/fine/final floor; sequential FINAL maximum with finite-look alpha-spent early statistical stop; "
    "continue through actual 5-cm thicknesses after an unresolved point; distinguish accepted VERIFIED-PASS design from exact minimum; "
    "authoritative formal bounds only; exact sufficient-statistic resume supported; bounded deep analog cannot veto merely for sparse support; "
    "full Geant4 output retained; no room geometry/detector dimensions"
)
# <<< v7.3k.12e16 SEQUENTIAL FINAL / SAFE-DESIGN CORE INTEGRATION <<<


t_total_start = time.perf_counter()

print(
    'beam_mev_input:',
    CFG.mev,
)

print(
    'source_energy_mev_internal:',
    derived_source_energy_mev(CFG),
)

t_sim_start = time.perf_counter()

results = run_simulation(
    CFG
)

t_sim_end = time.perf_counter()

print_summary(
    results
)

t_plot_start = time.perf_counter()

plot_results(
    results,
    CFG,
)

t_plot_end = time.perf_counter()

sweep = None

t_sweep_start = None
t_sweep_end = None
t_sweep_plot_start = None
t_sweep_plot_end = None

if CFG.run_thickness_sweep:
    t_sweep_start = time.perf_counter()

    sweep = run_thickness_sweep(
        CFG
    )

    t_sweep_end = time.perf_counter()

    t_sweep_plot_start = time.perf_counter()

    plot_results(
        results,
        CFG,
        sweep_results=sweep,
    )

    t_sweep_plot_end = time.perf_counter()

print(
    'OBJ export info:',
    results['obj_export'],
)

print(
    'Number of stored tracks:',
    len(results['tracks_cm']),
)

print(
    'Transmitted energy sample size:',
    len(
        results[
            'transmitted_energy_sample_mev'
        ]
    ),
)

print(
    'Exit cosine sample size:',
    len(
        results[
            'exit_cosine_sample'
        ]
    ),
)

t_total_end = time.perf_counter()

print(
    '\n==================== TIMING ===================='
)

print(
    f'Main simulation time          = '
    f'{t_sim_end - t_sim_start:.3f} s'
)

print(
    f'Main plotting time            = '
    f'{t_plot_end - t_plot_start:.3f} s'
)

if (
    CFG.run_thickness_sweep
    and t_sweep_start is not None
    and t_sweep_end is not None
):
    print(
        f'Thickness sweep time          = '
        f'{t_sweep_end - t_sweep_start:.3f} s'
    )

if (
    CFG.run_thickness_sweep
    and t_sweep_plot_start is not None
    and t_sweep_plot_end is not None
):
    print(
        f'Sweep plotting time           = '
        f'{t_sweep_plot_end - t_sweep_plot_start:.3f} s'
    )

print(
    f'TOTAL workflow time           = '
    f'{t_total_end - t_total_start:.3f} s'
)

print(
    '================================================\n'
)

results.keys()
results['summary']
results['process_counts_primary_total']
results['secondary_exit_hist_counts']['e-']
# >>> v7.3b 25CM-TO-5CM + SAME-ENERGY COMBINATION PATCH >>>
V73B_PHOTON_SWEEP_CONTRACT = "50-400 cm standard coarse grid; 25 cm coarse -> 5 cm direct refinement only"
# <<< v7.3b 25CM-TO-5CM + SAME-ENERGY COMBINATION PATCH <<<


# 2026-08-19-v7.3k.12c-report-numeric-canonicalization-lineage-cleanup: console terminology cleanup; transport/scoring unchanged.

# 2026-08-19-v7.3k.12d-cumulative-supervisory-plot-reuse: console plot provenance sync; transport/scoring unchanged.


