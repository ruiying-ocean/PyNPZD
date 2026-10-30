import numpy as np
import pandas as pd

from examples.compare_observations import (
    comparison_metrics,
    load_observations,
    monthly_comparison,
    to_chlorophyll,
)


def test_chlorophyll_conversion_uses_empower_units():
    np.testing.assert_allclose(
        to_chlorophyll(np.array([0.0, 1.0, 2.0]), 75.0), [0.0, 1.06, 2.12]
    )


def test_observations_preserve_stations_and_missing_months():
    observations = load_observations()
    counts = observations.groupby("station").observed_chlorophyll.count().to_dict()
    assert counts == {"biotrans": 12, "india": 9, "papa": 11, "kerfix": 10}
    assert len(observations) == 48
    assert observations.observed_nitrate.notna().all()
    january = observations[
        (observations.station == "biotrans") & (observations.month == 1)
    ].iloc[0]
    assert january.observed_chlorophyll == 0.21719593
    assert np.isclose(january.observed_nitrate, 5.947188645601273)


def test_monthly_metrics_use_means_and_exclude_missing_observations():
    model = pd.DataFrame(
        {
            "station": ["biotrans"] * 4,
            "month": [1, 1, 2, 2],
            "chlorophyll": [1.0, 3.0, 3.0, 5.0],
            "N": [9.0, 11.0, 11.0, 13.0],
        }
    )
    observed = pd.DataFrame(
        {
            "station": ["biotrans"] * 2,
            "month": [1, 2],
            "observed_chlorophyll": [1.0, np.nan],
            "observed_nitrate": [11.0, 11.0],
        }
    )
    paired = monthly_comparison(model, observed)
    np.testing.assert_array_equal(paired.model_chlorophyll, [2.0, 4.0])
    assert np.isnan(paired.observed_chlorophyll.iloc[1])
    metrics = comparison_metrics(paired).set_index("variable")
    assert metrics.loc["chlorophyll", "months"] == 1
    assert metrics.loc["chlorophyll", "rmse"] == 1.0
    assert metrics.loc["nitrate", "rmse"] == 1.0
    assert metrics.loc["nitrate", "bias"] == 0.0
