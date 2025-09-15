# 0D NPZD model re-implemented in Python
Nutrient-Phytoplankton-Zooplankton-Detritus Dynamic Model, adapted from EMPOWER (Anderson et al. 2015 GMD)

## Quick Start

### Prerequisites
- Python 3.8+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

### Installation with uv (recommended)

```bash
# Clone the repository
git clone <repository-url>
cd NPZD

# Install dependencies
uv sync

# Install with development dependencies for testing
uv sync --extra dev
```

### Installation with pip

```bash
# Clone the repository
git clone <repository-url>
cd NPZD

# Install dependencies
pip install -e .

# For development
pip install -e ".[dev]"
```

### Running the Model

```bash
# With uv
uv run python npzd.py

# With pip/regular Python
python npzd.py
```

### Running Tests

```bash
# With uv
uv run pytest tests/

# With pip/regular Python
python -m pytest tests/
```

### Building the Package

```bash
# Build wheel and source distribution
uv build
```

## Model Description

## Forcings
1. 3 main forcings: SST (T), MLD (Z), irradiance (I0) 
2. seasonal variation

## Derived parameters
Z -> N0 (deep-sea nutrient)
Z -> light (Beer's law)
Integrated photosynthesis rate => needs integration over depth and time

Kpar (light attenuation, Kpar)
option1: simple linear function (one layer)
option2: polynomial function considering depth

Light limitation (L_I) multi-layer
option1: numeric
option2: EP85
option3: Anderson93

P-I curve (implemented both for numeric L_I)
option1: Smith's function
option2: exponential function

## Configuration

The model is configured via `config.toml`. Key parameters include:

- `flag_stn`: Station selection (1-4)
  - Station 1: India (60°N 20°W)
  - Station 2: Biotrans (47°N 20°W)
  - Station 3: Kerfix (50°40'S 68°25'E)
  - Station 4: Papa (50°N 145°W)

- `flag_LI`: Light limitation method
  - 1: Numeric integration
  - 2: Evans & Parslow 1985
  - 3: Anderson 1993

- `flag_atten`: Light attenuation model
  - 1: Simple single-layer
  - 2: Multi-layer spectral

- `flag_force`: Forcing mode
  - `true`: Constant chemostat conditions
  - `false`: Time-varying environmental data

## Model Variables

The NPZD model tracks four state variables:
- **P**: Phytoplankton (mmol N m⁻³)
- **N**: Nutrients (mmol N m⁻³)
- **Z**: Zooplankton (mmol N m⁻³)
- **D**: Detritus (mmol N m⁻³)

## Project Structure

```
NPZD/
├── npzd.py              # Main model class
├── npzd_funcs.py        # Mathematical functions
├── config.toml          # Model configuration
├── stations_forcing.txt # Station forcing data
├── tests/               # Test suite
├── pyproject.toml       # Package configuration
└── CLAUDE.md           # Development guidance
```

## Reference

Anderson, T. R., Gentleman, W. C., and Sinha, B.: Influence of grazing formulations on the emergent properties of a complex ecosystem model in a global ocean general circulation model, Geosci. Model Dev., 8, 1681–1695, https://doi.org/10.5194/gmd-8-1681-2015, 2015.
