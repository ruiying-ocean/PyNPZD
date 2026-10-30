"""Utilities for loading external forcing datasets via OPeNDAP/xarray."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

import numpy as np
import xarray as xr


@dataclass
class DatasetConfig:
    """Configuration for a single remote dataset."""

    target: str
    variable: str
    url: str | None = None
    dataset_id: str | None = None
    engine: str | None = None
    time_coord: str = "time"
    lat_coord: str = "lat"
    lon_coord: str = "lon"
    scale: float = 1.0
    offset: float = 0.0
    wrap_longitude: bool = True
    selection: Mapping[str, Any] | None = None
    lat_margin: float = 0.25
    lon_margin: float = 0.25
    depth_min: float | None = None
    depth_max: float | None = None


class RemoteForcingLoader:
    """Load remote forcing data and provide monthly cycles or full time series for the NPZD model."""

    def __init__(self, config: Mapping[str, Any]):
        self.config = dict(config)
        datasets = self.config.get("datasets", [])
        if not datasets:
            raise ValueError("FORCING.remote.datasets must contain at least one entry")
        self.datasets = [self._build_dataset_cfg(item) for item in datasets]
        cache_dir = Path(self.config.get("cache_dir", "data/cache"))
        self.cache_dir = cache_dir.expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def load(self) -> Dict[str, np.ndarray] | Tuple[Dict[str, np.ndarray], np.ndarray]:
        """Fetch all configured datasets and return either monthly cycles or full time series."""
        latitude = self.config["latitude"]
        longitude = self.config["longitude"]
        start = np.datetime64(self.config["time_start"])
        end = np.datetime64(self.config["time_end"])

        result: Dict[str, np.ndarray] = {}
        time_coords = None

        for dataset in self.datasets:
            data, times = self._load_dataset(dataset, latitude, longitude, start, end)
            result[dataset.target.upper()] = data
            if times is not None and time_coords is None:
                time_coords = times

        if self.config.get("use_climatology", True):
            return result
        else:
            return result, time_coords

    def _build_dataset_cfg(self, raw: Mapping[str, Any]) -> DatasetConfig:
        data = dict(raw)
        selection = data.get("selection")
        if selection is not None and not isinstance(selection, Mapping):
            raise ValueError("FORCING.remote.datasets[].selection must be a mapping")
        dataset = DatasetConfig(
            target=data["target"],
            url=str(data["url"]) if "url" in data else None,
            dataset_id=data.get("dataset_id"),
            variable=data["variable"],
            engine=data.get("engine"),
            time_coord=data.get("time_coord", "time"),
            lat_coord=data.get("lat_coord", "lat"),
            lon_coord=data.get("lon_coord", "lon"),
            scale=float(data.get("scale", 1.0)),
            offset=float(data.get("offset", 0.0)),
            wrap_longitude=bool(data.get("wrap_longitude", True)),
            selection=dict(selection) if selection is not None else None,
            lat_margin=float(data.get("lat_margin", 0.25)),
            lon_margin=float(data.get("lon_margin", 0.25)),
            depth_min=data.get("depth_min"),
            depth_max=data.get("depth_max"),
        )

        if dataset.url is None and dataset.dataset_id is None:
            raise ValueError(
                "Each dataset entry must define either 'url' or 'dataset_id'"
            )

        return dataset

    def _load_dataset(
        self,
        dataset: DatasetConfig,
        latitude: float,
        longitude: float,
        start: np.datetime64,
        end: np.datetime64,
    ) -> Tuple[np.ndarray, np.ndarray | None]:
        cache_path = self._cache_path(dataset, start, end, latitude, longitude)
        if cache_path.exists():
            ds = xr.open_dataset(cache_path)
        else:
            if dataset.dataset_id:
                ds = self._open_copernicus_dataset(
                    dataset, latitude, longitude, start, end
                )
            else:
                open_kwargs = {"decode_times": False}
                if dataset.engine:
                    open_kwargs["engine"] = dataset.engine
                ds = xr.open_dataset(dataset.url, **open_kwargs)

                time_var = ds[dataset.time_coord]
                if not np.issubdtype(time_var.dtype, np.datetime64):
                    units = time_var.attrs.get("units")
                    calendar = time_var.attrs.get("calendar", "standard")
                    encoder = xr.coding.times.encode_cf_datetime
                    start_num = encoder(
                        np.array([start], dtype="datetime64[ns]"), units, calendar
                    )[0].item()
                    end_num = encoder(
                        np.array([end], dtype="datetime64[ns]"), units, calendar
                    )[0].item()
                    time_slice = slice(start_num, end_num)
                else:
                    time_slice = slice(start, end)

                ds = ds.sel({dataset.time_coord: time_slice})
                ds = self._select_point(ds, dataset.lat_coord, latitude, False)
                ds = self._select_point(
                    ds, dataset.lon_coord, longitude, dataset.wrap_longitude
                )
                ds = ds.load()
                try:
                    ds = xr.decode_cf(ds)
                except Exception:
                    pass
            if dataset.selection:
                ds = ds.sel(dataset.selection, method="nearest")
            ds.to_netcdf(cache_path)

        data_array = ds[dataset.variable]
        data_array = data_array.squeeze(drop=True)
        data_array = data_array * dataset.scale + dataset.offset

        # Get time coordinates if not using climatology
        time_coords = None
        if not self.config.get("use_climatology", True):
            time_coords = data_array[dataset.time_coord].values

        processed_data = self._to_monthly_cycle(data_array, dataset.time_coord)
        return processed_data, time_coords

    def _cache_path(self, dataset: DatasetConfig, start: np.datetime64, end: np.datetime64, latitude: float = None, longitude: float = None) -> Path:
        cache_key = {
            "url": str(dataset.url) if dataset.url is not None else None,
            "dataset_id": dataset.dataset_id,
            "variable": dataset.variable,
            "time_coord": dataset.time_coord,
            "lat_coord": dataset.lat_coord,
            "lon_coord": dataset.lon_coord,
            "latitude": latitude,
            "longitude": longitude,
            "selection": dataset.selection,
            "start": str(start),
            "end": str(end),
            "lat_margin": dataset.lat_margin,
            "lon_margin": dataset.lon_margin,
            "depth_min": dataset.depth_min,
            "depth_max": dataset.depth_max,
        }
        digest = hashlib.sha1(
            json.dumps(cache_key, sort_keys=True, default=str).encode()
        ).hexdigest()
        filename = f"{dataset.target.lower()}_{digest}.nc"
        return self.cache_dir / filename

    def _select_point(
        self,
        dataset: xr.Dataset,
        coord_name: str,
        value: float,
        wrap_longitude: bool = True,
    ) -> xr.Dataset:
        if coord_name not in dataset.coords:
            raise KeyError(f"Coordinate '{coord_name}' not found in dataset")
        coord = dataset[coord_name]
        target = value
        coord_values = coord.values
        if wrap_longitude and coord.ndim == 1:
            coord_min = np.nanmin(coord_values)
            coord_max = np.nanmax(coord_values)
            if coord_min >= 0 and value < 0:
                target = (value + 360.0) % 360.0
            elif coord_max <= 180 and value > 180:
                target = ((value + 180.0) % 360.0) - 180.0
        return dataset.sel({coord_name: target}, method="nearest")

    def _to_monthly_cycle(self, data_array: xr.DataArray, time_coord: str) -> np.ndarray:
        if time_coord not in data_array.dims:
            raise ValueError(
                f"Expected time coordinate '{time_coord}' in dimensions {data_array.dims}"
            )

        # Check if we want climatology or full time series
        use_climatology = self.config.get("use_climatology", True)

        if use_climatology:
            # Collapse multiple years into a climatological monthly cycle
            monthly = data_array.groupby(f"{time_coord}.month").mean(time_coord)
            monthly = monthly.sortby("month")
            values = monthly.values.astype(float)
            if values.shape[0] < 12:
                raise ValueError("Not enough data to construct a 12 month climatology")
            climatology = values[:12]
            climatology = np.concatenate([climatology, climatology[:1]])
            return climatology
        else:
            # Return full time series
            values = data_array.values.astype(float)
            return values

    def _open_copernicus_dataset(
        self,
        dataset: DatasetConfig,
        latitude: float,
        longitude: float,
        start: np.datetime64,
        end: np.datetime64,
    ) -> xr.Dataset:
        try:
            from copernicusmarine import open_dataset as copernicus_open_dataset
        except ImportError as exc:  # pragma: no cover - optional dependency guard
            raise ImportError(
                "copernicusmarine package is required for dataset entries with 'dataset_id'"
            ) from exc

        lat_min = latitude - dataset.lat_margin
        lat_max = latitude + dataset.lat_margin
        lon_min = longitude - dataset.lon_margin
        lon_max = longitude + dataset.lon_margin

        kwargs: Dict[str, Any] = {
            "dataset_id": dataset.dataset_id,
            "variables": [dataset.variable],
            "minimum_latitude": lat_min,
            "maximum_latitude": lat_max,
            "minimum_longitude": lon_min,
            "maximum_longitude": lon_max,
            "start_datetime": np.datetime_as_string(start, unit="s"),
            "end_datetime": np.datetime_as_string(end, unit="s"),
        }

        if dataset.depth_min is not None or dataset.depth_max is not None:
            kwargs["minimum_depth"] = dataset.depth_min
            kwargs["maximum_depth"] = dataset.depth_max

        ds = copernicus_open_dataset(**kwargs)

        ds = ds.load()
        ds = self._select_point(ds, dataset.lat_coord, latitude, False)
        ds = self._select_point(ds, dataset.lon_coord, longitude, dataset.wrap_longitude)

        if dataset.depth_min is not None or dataset.depth_max is not None:
            coord_name = "depth"
            if dataset.selection and "depth" in dataset.selection:
                coord_name = "depth"
            if coord_name in ds.coords:
                dmin = dataset.depth_min if dataset.depth_min is not None else ds[coord_name].min().item()
                dmax = dataset.depth_max if dataset.depth_max is not None else ds[coord_name].max().item()
                ds = ds.sel({coord_name: slice(dmin, dmax)})
                if coord_name in ds.dims and ds.dims[coord_name] > 1:
                    ds = ds.mean(coord_name)

        return ds


def load_remote_forcing(config: Mapping[str, Any]) -> Dict[str, np.ndarray]:
    """Helper that mirrors :class:`RemoteForcingLoader` for convenience."""

    loader = RemoteForcingLoader(config)
    return loader.load()
