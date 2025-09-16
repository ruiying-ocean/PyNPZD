import numpy as np
import pandas as pd
import xarray as xr

from forcing import RemoteForcingLoader
from npzd import NPZDModel


def _make_test_dataset(path):
    time = pd.date_range("2018-01-01", periods=24, freq="MS")
    lat = np.array([50.0])
    lon = np.array([-145.0])

    months = np.tile(np.arange(1, 13, dtype=float), 2)
    sst = months.reshape(-1, 1, 1)
    mld = (months * 10.0).reshape(-1, 1, 1)
    irr = (months * 100.0).reshape(-1, 1, 1)
    nut = (months * 2.0).reshape(-1, 1, 1)

    dataset = xr.Dataset(
        data_vars={
            "sst": (("time", "lat", "lon"), sst),
            "mld": (("time", "lat", "lon"), mld),
            "par": (("time", "lat", "lon"), irr),
            "nut": (("time", "lat", "lon"), nut),
        },
        coords={"time": time, "lat": lat, "lon": lon},
    )
    dataset.to_netcdf(path)


def test_remote_forcing_loader_returns_monthly_cycles(tmp_path):
    dataset_path = tmp_path / "test_forcing.nc"
    _make_test_dataset(dataset_path)

    config = {
        "latitude": 50.0,
        "longitude": -145.0,
        "time_start": "2018-01-01",
        "time_end": "2019-12-31",
        "cache_dir": tmp_path,
        "datasets": [
            {"target": "SST", "url": dataset_path, "variable": "sst", "engine": "netcdf4"},
            {"target": "MLD", "url": dataset_path, "variable": "mld", "engine": "netcdf4"},
            {"target": "IRRADIANCE", "url": dataset_path, "variable": "par", "engine": "netcdf4"},
            {"target": "NUTRIENT", "url": dataset_path, "variable": "nut", "engine": "netcdf4"},
        ],
    }

    loader = RemoteForcingLoader(config)
    data = loader.load()

    expected_sst = np.concatenate([np.arange(1, 13, dtype=float), [1.0]])
    expected_mld = np.concatenate([np.arange(1, 13, dtype=float) * 10.0, [10.0]])

    np.testing.assert_allclose(data["SST"], expected_sst)
    np.testing.assert_allclose(data["MLD"], expected_mld)
    assert data["IRRADIANCE"].shape == (13,)
    np.testing.assert_allclose(data["NUTRIENT"], np.concatenate([np.arange(1, 13, dtype=float) * 2.0, [2.0]]))


def test_model_initializes_with_remote_forcing(tmp_path):
    dataset_path = tmp_path / "test_forcing.nc"
    _make_test_dataset(dataset_path)

    config_text = f"""
[PARAMETER]
Vp0 = 2.5
alpha = 0.15
kN = 0.85
mP = 0.015
mP2 = 0.025
Imax = 1.0
kz = 0.6
phiP = 0.67
phiD = 0.33
betaz = 0.69
kNz = 0.75
mz = 0.02
mz2 = 0.34
VD = 6.43
mD = 0.06
wmix = 0.13
CtoChl = 75.0
albedo = 0.04
solarconst = 1368
parrac = 0.43
kw = 0.04
kc = 0.03

[SETUP]
Pinit = 0.1
Ninit = 10.0
Zinit = 0.1
Dsinit = 0.1
delta_t = 1
nyears = 1

flag_stn = 1
flag_LI = 1
flag_atten = 2
flag_irrad = 2
flag_PIcurve = 1
flag_grazing = 2
flag_force=false
force_SST=15
force_MLD=100
force_N0=5
force_irradiance=200

[FORCING]
source = "remote"

[FORCING.remote]
latitude = 50.0
longitude = -145.0
time_start = "2018-01-01"
time_end = "2019-12-31"
clouds = 6.0
e0 = 12.0
aN = 0.0
bN = 14.6
cache_dir = "{tmp_path}"

[[FORCING.remote.datasets]]
target = "SST"
url = "{dataset_path}"
variable = "sst"
engine = "netcdf4"

[[FORCING.remote.datasets]]
target = "MLD"
url = "{dataset_path}"
variable = "mld"
engine = "netcdf4"

[[FORCING.remote.datasets]]
target = "IRRADIANCE"
url = "{dataset_path}"
variable = "par"
engine = "netcdf4"

[[FORCING.remote.datasets]]
target = "NUTRIENT"
url = "{dataset_path}"
variable = "nut"
engine = "netcdf4"
"""

    config_path = tmp_path / "config.toml"
    config_path.write_text(config_text)

    model = NPZDModel(str(config_path))

    assert model.MLD.shape == (366,)
    assert model.SST.shape == (366,)
    assert model.noon_par_cycle.shape == (366,)
    assert model.remote_nutrient_cycle.shape == (366,)
    mean_nutrient = float(np.mean(model.remote_nutrient_cycle))
    assert np.isclose(mean_nutrient, model.station_settings["bN"], rtol=0.1)
