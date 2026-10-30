"""Compare the port with outputs generated independently by EMPOWER-1.0 in R."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from npzd_funcs import calc_aphy, calc_li_A93, calc_li_ep85, calc_li_num


pytestmark = pytest.mark.empower
REFERENCE = Path(__file__).parent / "fixtures/empower"
CHOICES = [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 2, 1), (3, 2, 1)]


def test_station_table_matches_original():
    actual = pd.read_csv(Path(__file__).parents[1] / "stations_forcing.txt")
    reference = pd.read_csv(REFERENCE / "station_table.csv")
    np.testing.assert_array_equal(actual.to_numpy(), reference.to_numpy())


@pytest.mark.parametrize("station", [1, 2, 3, 4])
def test_daily_forcing_matches_original(station_model, station):
    model = station_model
    model.settings["flag_stn"] = station
    model.initialize_model()
    reference = pd.read_csv(REFERENCE / "forcing.csv")
    reference = reference[reference.station == station]
    nutrient = model.get_deep_water_nitrate(
        model.MLD, model.station_settings["aN"], model.station_settings["bN"]
    )
    np.testing.assert_allclose(
        np.column_stack([model.MLD, model.SST, nutrient]),
        reference[["MLD", "SST", "N0"]],
        rtol=1e-10,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        np.column_stack([model.day_length_cycle, model.noon_par_cycle]),
        reference[["day_length", "noon_par"]].values[:365],
        rtol=1e-10,
        atol=1e-10,
    )


@pytest.mark.parametrize("irradiance,pi_curve", [(1, 1), (1, 2), (2, 1), (2, 2)])
def test_numeric_light_matches_original(irradiance, pi_curve):
    reference = pd.read_csv(REFERENCE / "light.csv")
    reference = reference[(reference.irrad == irradiance) & (reference.PI == pi_curve)]
    actual = [
        calc_li_num(r.depth, r.Iin, r.Iout, 0.08, 0.15, r.Vp, 12, irradiance, pi_curve)
        for r in reference.itertuples()
    ]
    np.testing.assert_allclose(actual, reference.num, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("method", ["ep85", "a93"])
def test_analytic_light_matches_original(method):
    reference = pd.read_csv(REFERENCE / "light.csv")
    reference = reference[(reference.irrad == 1) & (reference.PI == 1)]
    actual = []
    for row in reference.itertuples():
        args = (row.depth, row.Iin, row.Iout, 0.08, 0.15, row.Vp, 12)
        actual.append(
            calc_li_ep85(*args) if method == "ep85" else calc_li_A93(*args, 0.3)
        )
    np.testing.assert_allclose(actual, reference[method], rtol=1e-10, atol=1e-10)


def test_chlorophyll_absorption_matches_original():
    reference = pd.read_csv(REFERENCE / "absorption.csv")
    actual = [
        calc_aphy(r.ss, r.top, r.top + r.thickness, 0.3) for r in reference.itertuples()
    ]
    np.testing.assert_allclose(actual, reference.value, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("li,irradiance,pi_curve", CHOICES)
@pytest.mark.parametrize("attenuation", [1, 2])
@pytest.mark.parametrize("grazing", [1, 2])
def test_ecosystem_rates_match_original(
    station_model, li, irradiance, pi_curve, attenuation, grazing
):
    model = station_model
    model.settings.update(
        flag_LI=li,
        flag_irrad=irradiance,
        flag_PIcurve=pi_curve,
        flag_atten=attenuation,
        flag_grazing=grazing,
    )
    reference = pd.read_csv(REFERENCE / "rates.csv")
    reference = reference[
        (reference.LI == li)
        & (reference.irrad == irradiance)
        & (reference.PI == pi_curve)
        & (reference.atten == attenuation)
        & (reference.grazing == grazing)
    ]
    assert len(reference) >= 60
    for row in reference.itertuples():
        model.MLD = np.full(366, float(row.depth))
        model.MLD[row.day] += row.change
        model.SST = np.full(366, float(row.temperature))
        actual = model.ode([row.P, row.N, row.Z, row.D], row.day - 1 + 0.5)
        np.testing.assert_allclose(
            actual,
            [row.dP, row.dN, row.dZ, row.dD],
            rtol=1e-10,
            atol=1e-10,
            equal_nan=False,
            err_msg=str(row),
        )


@pytest.mark.slow
def test_five_years_with_original_rk4_schedule(station_model):
    """Hold forcing fixed across each RK4 update, exactly as the R driver does."""
    model = station_model
    reference = pd.read_csv(REFERENCE / "biotrans_5years.csv")
    base_mld = model.MLD.copy()
    state = np.array([model.settings[k] for k in ("Pinit", "Ninit", "Zinit", "Dsinit")])
    daily = [state.copy()]
    for day in range(len(reference) - 1):
        d = day % 365
        change = base_mld[d + 1] - base_mld[d]
        for step in range(10):
            depth = base_mld[d] + step / 10 * change
            # Feed the original within-step depth and deepening to the unchanged ODE.
            model.MLD[d : d + 2] = [depth, depth + change]
            k1 = np.array(model.ode(state, float(day)))
            k2 = np.array(model.ode(state + 0.05 * k1, float(day)))
            k3 = np.array(model.ode(state + 0.05 * k2, float(day)))
            k4 = np.array(model.ode(state + 0.1 * k3, float(day)))
            state = state + 0.1 * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        daily.append(state.copy())
    assert np.all(np.isfinite(daily)) and np.all(np.asarray(daily) >= 0)
    np.testing.assert_allclose(
        daily,
        reference[["P", "N", "Z", "D"]],
        rtol=1e-9,
        atol=1e-9,
    )


@pytest.mark.slow
def test_production_solver_stays_close_to_original(station_model):
    model = station_model
    reference = pd.read_csv(REFERENCE / "biotrans_5years.csv")
    times, solution = model.run()
    assert np.all(np.isfinite(solution)) and np.all(solution >= 0)
    # Python omits the final endpoint and saves 10 samples/day in this configuration.
    np.testing.assert_allclose(times[::10], reference.time.values[:-1], atol=1e-8)
    np.testing.assert_allclose(
        solution[::10],
        reference[["P", "N", "Z", "D"]].values[:-1],
        rtol=0,
        atol=5e-3,
    )


def test_generator_rejects_changed_reference_archive(tmp_path):
    archive = tmp_path / "changed.zip"
    archive.write_bytes(b"not the published supplement")
    output = tmp_path / "output"
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "reference/generate.py"),
            "--archive",
            str(archive),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0 and "SHA-256" in result.stderr
    assert not output.exists()
