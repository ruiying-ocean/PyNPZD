# Seasonal model–observation comparison

Run from the repository root using the installed project environment:

```bash
uv run --no-sync python examples/compare_observations.py
```

The script runs five years at BIOTRANS, India, Papa, and KERFIX, then plots
chlorophyll and nitrate from the final year against observation-based monthly
climatologies. It saves a 300 dpi PNG, vector PDF, the model's last-year output,
paired monthly values, and RMSE/bias tables under
`data/examples/empower_observations/`. These generated outputs are git-ignored.
Use `--output /path/to/folder` to save them elsewhere. No downloads or
Copernicus credentials are needed to run the example.

## Model setup

`empower_observations.toml` contains the original EMPOWER-1.0 baseline inputs.
The simulation uses the original `stations_forcing.txt`, numerical light
integration, sinusoidal irradiance, the Smith P–I curve, and three-layer light
attenuation. The Python adaptive solver is retained; output spacing is 2.4 h.
The first four years are discarded. No parameters are fitted to the plotted
observations.

The script applies the station-specific values from Table 3 of
[Anderson, Gentleman, and Yool (2015)](https://doi.org/10.5194/gmd-8-2231-2015):

| Station | `Vp0` | `alpha` | `Imax` | `mz` |
| --- | ---: | ---: | ---: | ---: |
| BIOTRANS | 2.5 | 0.15 | 1.0 | 0.02 |
| India | 2.5 | 0.15 | 1.0 | 0.0 |
| Papa | 1.25 | 0.075 | 1.25 | 0.02 |
| KERFIX | 1.25 | 0.075 | 2.0 | 0.02 |

All other parameters and initial concentrations follow the baseline. The
five-year duration follows the reference example; it is not a claim that every
possible parameter set equilibrates in four years.

## Observation source

`data/phydra_stations_verification.csv` is an unchanged copy of the small
[verification table distributed by Phydra](https://github.com/ben1post/phydra/blob/3a69dc90bd5c48526c0b4cd4c2fb6bbddc88cc5e/notebooks/data/stations_verification.csv).
The corresponding [model notebook](https://github.com/ben1post/phydra/blob/3a69dc90bd5c48526c0b4cd4c2fb6bbddc88cc5e/notebooks/2_NPZD_Slab_ModelApplication.ipynb)
and [Post et al. (2024), section 3.2.3](https://doi.org/10.5194/gmd-17-1175-2024)
document how it is used:

- Chlorophyll: MODIS-Aqua monthly climatology, retrieved through August 2022.
- Nitrate: WOA2018 climatology averaged through the mixed layer using IFREMER
  mixed-layer depths.
- Rows 0–11 represent January–December. `Chla` is mg chlorophyll m⁻³; `N` is
  µmol N L⁻¹, numerically equal to the model's mmol N m⁻³.

Downloaded on 2026-10-30 from commit
`3a69dc90bd5c48526c0b4cd4c2fb6bbddc88cc5e`. The CSV SHA-256 is
`2ba64dfbd20a7926fede0faa14fb102021d5db2d6be7ec08deed96e691787812`.
The upstream BSD-3-Clause notice is retained as `data/PHYDRA_LICENSE.txt`.

This is an illustrative seasonal comparison, not an exact recreation of the
2015 figures. That paper used WOA2009 nitrate and selected years of 8-day
SeaWiFS chlorophyll. Here, updated climatologies are compared with the
original EMPOWER forcing, so the periods and mixed-layer averaging are not
identical. Monthly climatologies also smooth blooms whose timing varies
between years. Satellite chlorophyll is a surface estimate, whereas the model
represents a homogeneous mixed layer.

The root `Chl_obs.csv` is not used: its source and averaging period were not
documented, and its values differ from this verifiable table.

## Transformations and figure interpretation

Model phytoplankton is converted to chlorophyll using the same convention as
the original R code and this port:

```text
chlorophyll [mg m^-3] = P [mmol N m^-3] × 6.625 × 12.0 / CtoChl
CtoChl = 75 g C / g chlorophyll
```

The blue curve shows the simulated final-year cycle. Open blue squares show
model monthly means, and orange-red circles show the observational climatology.
Months use a 365-day calendar; model samples are assigned to half-open
calendar-month intervals and averaged at their uniform 2.4-hour spacing.
Observed points are positioned at month centres. The underlying values are
included in `monthly_comparison.csv`.

Missing satellite months remain NaN and are neither filled nor joined by
lines. There are 12, 9, 11, and 10 chlorophyll months at BIOTRANS, India, Papa,
and KERFIX respectively, and 12 nitrate months at each. RMSE and bias use only
paired monthly means; bias is model minus observation. Both have the units of
their variable. No error bars are drawn because the table supplies no
uncertainty estimates. All stations share the same vertical scale within
each variable, and zero is included.

The plot is a diagnostic, not a unit-test acceptance threshold or proof of
model skill. The tests check conversions, month averaging, and missing-data
handling rather than requiring an arbitrary fit to observations.

## Software acknowledgment

Plots use NumPy, pandas, and Matplotlib already installed for this project.
Figure preparation for this example and `plot_station_cycle.py` followed the
scientific-visualization skill from
Kassis, T., Agarwal, V., He, Y., Patel, D., and Brueckner, A. M. (2026),
*Scientific Agent Skills: A Library of Procedural Knowledge for Research
Agents*, [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065).
