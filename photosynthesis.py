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
    g = np.array([0.048014, 0.00023779, -0.023074, 0.0031095, -0.0090545, 
                  0.0027974, 0.00085217, -3.9804E-06, 0.0012398, -0.00061991])

    x_bottom = zbottom + 1.0
    xlg_bottom = np.log(x_bottom)
    termf1a = x_bottom * xlg_bottom - x_bottom
    termf2a = x_bottom * xlg_bottom**2 - 2 * x_bottom * xlg_bottom + 2 * x_bottom
    termf3a = x_bottom * xlg_bottom**3 - 3 * x_bottom * xlg_bottom**2 + 6 * x_bottom * xlg_bottom - 6 * x_bottom

    x_top = ztop + 1.0
    xlg_top = np.log(x_top)
    termf1b = x_top * xlg_top - x_top
    termf2b = x_top * xlg_top**2 - 2 * x_top * xlg_top + 2 * x_top
    termf3b = x_top * xlg_top**3 - 3 * x_top * xlg_top**2 + 6 * x_top * xlg_top - 6 * x_top

    terma = g[0] + g[1] * ss + g[4] * ss**2 + g[6] * ss**3
    termb = g[2] + g[3] * ss + g[8] * ss**2
    termc = g[5] + g[9] * ss

    acalc = (zbottom + 1.0) * terma + termf1a * termb + \
            termf2a * termc + termf3a * g[7] - \
            ((ztop + 1.0) * terma + termf1b * termb + \
            termf2b * termc + termf3b * g[7])

    return acalc
