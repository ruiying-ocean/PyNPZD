# 0D NPZD model re-implemented in Python
Nutrient-Phytoplankton-Zooplankton-Detritus Dynamic Model, adapted from EMPOWER (Anderson et al. 2015 GMD)

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
