import numpy as np
import pandas as pd
import pytest


@pytest.mark.parametrize("next_depth", [90.0, 100.0, 110.0])
@pytest.mark.parametrize("time", [0.0, 0.5, 1.5])
def test_time_series_entrainment_and_dilution(station_model, next_depth, time):
    model = station_model
    model.use_time_series = True
    model.MLD = np.array([100.0, next_depth])
    model.SST = np.array([10.0, 10.0])
    model.station_settings.update(aN=0.0, bN=14.6)
    state = np.array([0.1, 2.0, 0.2, 0.3])
    actual = np.array(model.ode(state, time))

    depth = 100.0 + (next_depth - 100.0) * min(time, 1.0)
    model.MLD = np.full(2, depth)
    without_deepening = np.array(model.ode(state, time))
    deepening = max(next_depth - 100.0, 0.0) if time < 1.0 else 0.0
    expected = deepening / depth * np.array([-0.1, 12.6, -0.2, -0.3])
    np.testing.assert_allclose(actual - without_deepening, expected, atol=1e-12)


@pytest.mark.parametrize(
    "start,time,day",
    [
        (None, 0.5, 1),
        ("2001-01-30", 0.5, 30),
        ("2000-02-29", 0.5, 60),
        ("2000-12-31", 1.5, 1),
    ],
)
def test_time_series_sunlight_matches_calendar_day(station_model, start, time, day):
    model = station_model
    model.MLD = np.full(366, 100.0)
    model.SST = np.full(366, 10.0)
    state = [0.1, 2.0, 0.2, 0.3]
    expected = model.ode(state, day - 1 + 0.5)

    model.use_time_series = True
    if start is not None:
        model.time_coords = pd.date_range(start, periods=2).values
    np.testing.assert_allclose(model.ode(state, time), expected, atol=1e-12)


def test_time_series_rejects_non_daily_forcing(station_model):
    model = station_model
    model.use_time_series = True
    model.time_coords = pd.date_range("2001-01-01", periods=13, freq="MS").values
    with pytest.raises(ValueError, match="daily"):
        model.interpolate_forcing_data()
