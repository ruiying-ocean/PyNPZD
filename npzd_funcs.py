import numpy as np
from scipy.special import factorial

try:  # Optional acceleration via Numba when available.
    from numba import njit
except ModuleNotFoundError:  # pragma: no cover - fallback for environments without numba
    def njit(*args, **kwargs):  # type: ignore
        def decorator(func):
            return func

        if args and callable(args[0]):
            return args[0]
        return decorator


# Precomputed constants reused inside the light limitation helpers. Using module
# level arrays avoids reallocating the same small vectors for every ODE function
# evaluation, which noticeably reduces the total runtime of the integrator.
_NUM_TIME_STEPS = 10
_TRAPEZOID_WEIGHTS = np.ones(_NUM_TIME_STEPS, dtype=np.float64)
_TRAPEZOID_WEIGHTS[-1] = 0.5
_TRIANGULAR_IRRADIANCE = np.linspace(1.0, float(_NUM_TIME_STEPS), _NUM_TIME_STEPS, dtype=np.float64) / _NUM_TIME_STEPS
_SINUSOIDAL_IRRADIANCE = np.sin(np.arange(1.0, _NUM_TIME_STEPS + 1.0, dtype=np.float64) / _NUM_TIME_STEPS * (np.pi / 2.0))

_SERIES_INDICES = np.arange(1, 17, dtype=np.int64)
_SERIES_COEFFS = ((-1.0) ** (_SERIES_INDICES + 1)) / (_SERIES_INDICES * factorial(_SERIES_INDICES))

_A93_OMEGA = np.array([1.9004, -2.8333E-01, 2.8050E-02, -1.4729E-03, 3.0841E-05], dtype=np.float64)
_A93_POWERS = np.arange(1, 6, dtype=np.int64)

_APHY_COEFF = np.array([
    0.048014, 0.00023779, -0.023074, 0.0031095, -0.0090545,
    0.0027974, 0.00085217, -3.9804E-06, 0.0012398, -0.00061991
], dtype=np.float64)


@njit(cache=True, fastmath=True)
def calc_li_ep85(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow):
    """
    L_I: Evans & Parslow (1985) function for calculating daily photosynthesis
    Light extinction using Beer's Law and assumed triangular irradiance during day_in_year

    ## validated
    """
    tau1 = daylnow / 2.0 / 24.0
    betaEP = Vp * tau1 / (alpha * Iin)
    betaEP2 = Vp * tau1 / (alpha * Iout)
    
    sqrt_betaEP2_tau1 = np.sqrt(betaEP2**2 + tau1**2)
    sqrt_betaEP_tau1 = np.sqrt(betaEP**2 + tau1**2)
    
    Lim_I = 2.0 * Vp / (zdepth * kPARlay) * (
        sqrt_betaEP2_tau1
        - tau1 * np.log((tau1 + sqrt_betaEP2_tau1) / betaEP2)
        - (sqrt_betaEP_tau1 - tau1 * np.log((tau1 + sqrt_betaEP_tau1) / betaEP))
        - betaEP2
        + betaEP
    )
    
    # Divide by Vp to make dimensionless
    Lim_I /= Vp
    
    return Lim_I

@njit(cache=True, fastmath=True)
def calc_li_num(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, choiceirrad, choicePI):
    if choiceirrad == 1:
        irradiance = _TRIANGULAR_IRRADIANCE
    elif choiceirrad == 2:
        irradiance = _SINUSOIDAL_IRRADIANCE
    else:
        raise ValueError("Unsupported irradiance choice")

    prefactor = Vp / (kPARlay * zdepth)

    if choicePI == 1:
        x0 = alpha * irradiance * Iin
        xH = x0 * np.exp(-kPARlay * zdepth)
        sqrt0 = np.sqrt(Vp * Vp + x0 * x0)
        sqrtH = np.sqrt(Vp * Vp + xH * xH)
        sumpsz = prefactor * (np.log(x0 + sqrt0) - np.log(xH + sqrtH))
    else:
        istar1 = alpha * irradiance * Iin / Vp
        istar2 = alpha * irradiance * Iout / Vp
        diff = np.empty((_NUM_TIME_STEPS, _SERIES_INDICES.shape[0]))
        for i in range(_NUM_TIME_STEPS):
            for j in range(_SERIES_INDICES.shape[0]):
                diff[i, j] = istar1[i] ** _SERIES_INDICES[j] - istar2[i] ** _SERIES_INDICES[j]
        coeff_sum = diff @ _SERIES_COEFFS
        sumpsz = coeff_sum * prefactor

    sumps = 0.0
    for i in range(_NUM_TIME_STEPS):
        sumps += sumpsz[i] * _TRAPEZOID_WEIGHTS[i]
    Lim_I = (sumps / _NUM_TIME_STEPS) * daylnow / 24.0
    return Lim_I / Vp



@njit(cache=True, fastmath=True)
def calc_li_A93(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, ahashnow):
    # Calculate alphamax
    alphamax = alpha * 2.602  # alphamax is alpha at wavelength of maximum absorption cross section
    
    # Calculate daily photosynthesis in each layer, mg C m-2 d-1
    V0 = daylnow * Vp / (np.pi * kPARlay)
    V1 = alphamax * ahashnow * Iin / Vp
    V2 = alphamax * ahashnow * Iout / Vp
    
    V_terms = np.empty(_A93_POWERS.shape[0])
    for i in range(_A93_POWERS.shape[0]):
        V_terms[i] = V1 ** _A93_POWERS[i] - V2 ** _A93_POWERS[i]
    Qpsnow = V0 * np.dot(_A93_OMEGA, V_terms)
    
    # Convert to dimensionless units to get 0 <= J <= 1
    Lim_I = Qpsnow / zdepth / 24.0 / Vp
    
    return Lim_I


@njit(cache=True, fastmath=True)
def calc_aphy(ss, ztop, zbottom, aphylast):
    g = _APHY_COEFF

    # Calculate terms for bottom depth
    x = zbottom + 1.0
    xlg = np.log(x)
    termf1a = x * xlg - x
    termf2a = x * xlg**2 - 2 * x * xlg + 2 * x
    termf3a = x * xlg**3 - 3 * x * xlg**2 + 6 * x * xlg - 6 * x
    
    # Calculate terms for top depth
    x = ztop + 1.0
    xlg = np.log(x)
    termf1b = x * xlg - x
    termf2b = x * xlg**2 - 2 * x * xlg + 2 * x
    termf3b = x * xlg**3 - 3 * x * xlg**2 + 6 * x * xlg - 6 * x
    
    # Calculate composite terms
    terma = g[0] + g[1]*ss + g[4]*ss**2 + g[6]*ss**3
    termb = g[2] + g[3]*ss + g[8]*ss**2
    termc = g[5] + g[9]*ss
    
    # Calculate final result
    acalc = ((zbottom + 1.0) * terma + 
             termf1a * termb + 
             termf2a * termc + 
             termf3a * g[7] - 
             ((ztop + 1.0) * terma + 
              termf1b * termb + 
              termf2b * termc + 
              termf3b * g[7]))
    
    return acalc

def day_length(day_in_year, lat_radian):
    """calculate day_in_year length in hours"""
    declin = 23.45 * np.sin(2 * np.pi * (284 + day_in_year) * 0.00274) * np.pi / 180
    daylnow =  24 * np.arccos(-np.tan(lat_radian) * np.tan(declin)) / np.pi
        
    return daylnow

def noon_par(jday, latradians, clouds, e0):
    # Constants
    albedo = 0.04          # albedo
    solarconst = 1368.0    # solar constant, w m-2
    parrac = 0.43          # PAR fraction
    
    # Calculate solar declination angle
    declin = 23.45 * np.sin(2 * np.pi * (284 + jday) * 0.00274) * np.pi / 180
    
    # Calculate cosine of zenith angle
    coszen = np.sin(latradians) * np.sin(declin) + \
             np.cos(latradians) * np.cos(declin)
    
    # Calculate zenith angle in degrees
    zen = np.arccos(coszen) * 180 / np.pi
    
    # Calculate Earth's radius vector
    Rvector = 1 / np.sqrt(1 + 0.033 * np.cos(2 * np.pi * jday * 0.00274))
    
    # Calculate irradiance at ocean surface, clear sky
    Iclear = (solarconst * coszen**2) / \
             (Rvector**2) / \
             (1.2 * coszen + e0 * (1.0 + coszen) * 0.001 + 0.0455)
    
    # Calculate cloud factor (atmospheric transmission)
    cfac = (1 - 0.62 * clouds * 0.125 + 0.0019 * (90 - zen))
    
    # Calculate noon irradiance: total solar
    Inoon = Iclear * cfac * (1 - albedo)
    
    # Calculate final noon PAR
    noonparnow = parrac * Inoon
    
    return noonparnow


def interpolate_monthly_to_daily(monthly_data):
    daily_data = np.zeros(366)
    
    for istep in range(1,366):
        rdiv = 365.0/12.0
        quot = istep // rdiv
        remdr = (istep % rdiv) / rdiv
        
        n1 = int(quot + 1)
        if n1 == 13:
            n1 = 1
        n2 = n1 + 1
        
        iday = istep + 15
        if iday > 365:
            iday = iday - 365
            
        daily_data[iday-1] = monthly_data[n1-1] + (monthly_data[n2-1] - monthly_data[n1-1]) * remdr
    
    ## shift by one place
    daily_data =np.roll(daily_data, 1)

    daily_data[0] = daily_data[-1]
    
    return daily_data
