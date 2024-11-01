import numpy as np
from scipy.special import factorial


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

def calc_li_num(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, choiceirrad, choicePI):
    sumps = 0.0  # for calculating sum of photosynthesis
    range_rad = np.pi / 2.0
    nstepst = 10  # the day from sunrise to midday will be divided into nstepst steps (can be altered if so desired)
    
    for itme in range(1, nstepst + 1):  # time loop (no need to calculate itme = 0 because irradiance zero at sunrise)
        rtrapt = 0.5 if itme == nstepst else 1.0  # trapezoidal rule
        
        if choiceirrad == 1:  # choice of either (1) triangular or (2) sinusoidal irradiance over day
            Icorr = itme / nstepst  # irradiance multiplier for time of day: triangular
        elif choiceirrad == 2:
            Icorr = np.sin(itme / nstepst * range_rad)  # irradiance multiplier for time of day: sinusoidal
        
        if choicePI == 1:  # Smith function (the most straightforward choice in terms of calculating analytic depth integral)
            x0 = alpha * Icorr * Iin
            xH = alpha * Icorr * Iin * np.exp(-kPARlay * zdepth)
            sumpsz = Vp / kPARlay / zdepth * (np.log(x0 + np.sqrt(Vp**2 + x0**2)) - np.log(xH + np.sqrt(Vp**2 + xH**2)))
        else:  # exponential function (Platt et al., 1990)
            choiceanalytic = 1  # 1=analytic solution; 2=numeric solution; set manually in the code
            if choiceanalytic == 1:  # analytic solution
                istar1 = alpha * Icorr * Iin / Vp
                istar2 = alpha * Icorr * Iout / Vp
                # Analytic solution is an infinite series
                ifact = np.arange(1, 17)
                sumpsz = np.sum((-1)**(ifact + 1) / (ifact * factorial(ifact)) * (istar1**ifact - istar2**ifact))
                sumpsz *= Vp / kPARlay / zdepth
            else:  # numeric solution
                nstepsz = 100  # no. of steps: depth
                izz = np.arange(nstepsz + 1)
                rtrapz = np.ones(nstepsz + 1)
                rtrapz[0] = rtrapz[-1] = 0.5  # trapezoidal rule
                
                z = izz / nstepsz * zdepth
                Iz = Iin * Icorr * np.exp(-kPARlay * z)
                psnow = Vp * (1.0 - np.exp(-alpha * Iz / Vp))
                sumpsz = np.sum(psnow * rtrapz) / nstepsz
        
        sumps += sumpsz * rtrapt  # sum photosynthesis over daylight hours
    
    Lim_I = sumps / nstepst * daylnow / 24.0  # take account of zero PS during night
    Lim_I /= Vp  # dimensionless
    return Lim_I



def calc_li_A93(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, ahashnow):
    # Polynomial coefficients for calculating photosynthesis (Platt et al., 1990)
    omeg = np.array([1.9004, -2.8333E-01, 2.8050E-02, -1.4729E-03, 3.0841E-05])
    
    # Calculate alphamax
    alphamax = alpha * 2.602  # alphamax is alpha at wavelength of maximum absorption cross section
    
    # Calculate daily photosynthesis in each layer, mg C m-2 d-1
    V0 = daylnow * Vp / (np.pi * kPARlay)
    V1 = alphamax * ahashnow * Iin / Vp
    V2 = alphamax * ahashnow * Iout / Vp
    
    # Calculate Qpsnow using vectorized operations
    powers = np.arange(1, 6)
    V_terms = V1**powers - V2**powers
    Qpsnow = V0 * np.sum(omeg * V_terms)
    
    # Convert to dimensionless units to get 0 <= J <= 1
    Lim_I = Qpsnow / zdepth / 24.0 / Vp
    
    return Lim_I


def calc_aphy(ss, ztop, zbottom, aphylast):

    # coefficients for calculating a
    g = np.array([
        0.048014, 0.00023779, -0.023074, 0.0031095, -0.0090545,
        0.0027974, 0.00085217, -3.9804E-06, 0.0012398, -0.00061991
    ])
    
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
