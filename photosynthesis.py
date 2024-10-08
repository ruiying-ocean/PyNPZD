from scipy.special import factorial
import numpy as np

def fnli_calc_ep85(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow):
    """
    L_I: Evans & Parslow (1985) function for calculating daily photosynthesis
    Light extinction using Beer's Law and assumed triangular irradiance during day
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

def fnaphy(ss, ztop, zbottom, aphylast):
    """
    L_I: Anderson (1993) function for calculating daily photosynthesis                    
    Uses in-built piecewise light extinction and assumes sinusoidal irradiance during day
    Water column divided into three layers: 0-5m, 5-23m and >23m                          
    """
    # coeffs. for calculating a#
    g = np.array([0.048014, 0.00023779, -0.023074, 0.0031095, -0.0090545, 
                  0.0027974, 0.00085217, -3.9804E-06, 0.0012398, -0.00061991])
    
    x_bottom = zbottom + 1.0
    xlg_bottom = np.log(x_bottom)
    
    x_top = ztop + 1.0
    xlg_top = np.log(x_top)
    
    def calculate_terms(x, xlg):
        termf1 = x * xlg - x
        termf2 = x * xlg**2 - 2 * x * xlg + 2 * x
        termf3 = x * xlg**3 - 3 * x * xlg**2 + 6 * x * xlg - 6 * x
        return termf1, termf2, termf3
    
    termf1a, termf2a, termf3a = calculate_terms(x_bottom, xlg_bottom)
    termf1b, termf2b, termf3b = calculate_terms(x_top, xlg_top)
    
    terma = g[0] + g[1] * ss + g[4] * ss**2 + g[6] * ss**3
    termb = g[2] + g[3] * ss + g[8] * ss**2
    termc = g[5] + g[9] * ss
    
    acalc = (x_bottom * terma + termf1a * termb + termf2a * termc + termf3a * g[7]) - \
            (x_top * terma + termf1b * termb + termf2b * termc + termf3b * g[7])
    
    return acalc

def fnli_calc_a93(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, ahashnow):
    """
    Photosynthesis calculated using polynomial approximation (Anderson, 1993)
    This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis
    """
    
    # Polynomial coefficients for calculating photosynthesis (Platt et al., 1990)
    omeg = np.array([1.9004, -2.8333E-01, 2.8050E-02, -1.4729E-03, 3.0841E-05])
    
    # Calculate alphamax
    alphamax = alpha * 2.602  # alphamax is alpha at wavelength of maximum absorption cross section
    
    # Calculate daily photosynthesis in each layer, mg C m-2 d-1
    V0 = daylnow * Vp / (np.pi * kPARlay)
    V1 = alphamax * ahashnow * Iin / Vp
    V2 = alphamax * ahashnow * Iout / Vp
    
    # Calculate powers of V1 and V2
    V1_powers = np.power(V1, np.arange(1, 6))
    V2_powers = np.power(V2, np.arange(1, 6))
    
    # Calculate Qpsnow using vectorized operations
    Qpsnow = V0 * np.sum(omeg * (V1_powers - V2_powers))
    
    # Convert to dimensionless units to get 0 <= J <= 1
    Lim_I = Qpsnow / zdepth / 24.0 / Vp
    
    return Lim_I


def fnli_calc_num(zdepth, Iin, Iout, kPARlay, alpha, Vp, daylnow, choice_irrad, choice_PI):
    """
    L_I: Calculation of daily photosynthesis: numeric integration (time) and analytic integration (depth)
    Diel cycle of irradiance: choice of (1) triangular vs (2) sinusoidal (flag_irrad)
    P-I curve: choice of (1) Smith fn vs (2) exponential fn  (flag_PIcurve)
    numeric integration (time) uses nstepst+1 steps and trapezoidal rule for first and last steps
    (the user may manually change nstepst in the code if so desired)
    """
    sumps = 0.0  # for calculating sum of photosynthesis
    range_rad = np.pi / 2.0
    nstepst = 10  # the day from sunrise to midday will be divided into nstepst steps (can be altered if so desired)

    def calc_icorr(itme, choice):
        if choice == 1:
            return itme / nstepst  # irradiance multiplier for time of day: triangular
        elif choice == 2:
            return np.sin(itme / nstepst * range_rad)  # irradiance multiplier for time of day: sinusoidal

    def smith_function(Icorr):
        x0 = alpha * Icorr * Iin
        xH = alpha * Icorr * Iin * np.exp(-kPARlay * zdepth)
        return Vp / kPARlay / zdepth * (np.log(x0 + np.sqrt(Vp**2 + x0**2)) - np.log(xH + np.sqrt(Vp**2 + xH**2)))

    def exponential_function_analytic(Icorr):
        istar1 = alpha * Icorr * Iin / Vp
        istar2 = alpha * Icorr * Iout / Vp
        sumpsz = np.sum([((-1)**(ifact + 1) / (ifact * factorial(ifact)) * (istar1**ifact - istar2**ifact)) for ifact in range(1, 17)])
        return sumpsz * Vp / kPARlay / zdepth

    def exponential_function_numeric(Icorr):
        nstepsz = 100  # no. of steps: depth
        z = np.linspace(0, zdepth, nstepsz + 1)
        Iz = Iin * Icorr * np.exp(-kPARlay * z)
        psnow = Vp * (1.0 - np.exp(-alpha * Iz / Vp))
        weights = np.ones(nstepsz + 1)
        weights[0] = weights[-1] = 0.5  # trapezoidal rule
        return np.sum(psnow * weights) / nstepsz

    for itme in range(1, nstepst + 1):
        rtrapt = 0.5 if itme == nstepst else 1.0
        Icorr = calc_icorr(itme, choice_irrad)

        if choice_PI == 1:
            sumpsz = smith_function(Icorr)
        else:
            choice_analytic = 1  # 1=analytic solution; 2=numeric solution; set manually in the code
            if choice_analytic == 1:
                sumpsz = exponential_function_analytic(Icorr)
            else:
                sumpsz = exponential_function_numeric(Icorr)

        sumps += sumpsz * rtrapt

    Lim_I = sumps / nstepst * daylnow / 24.0  # take account of zero PS during night
    Lim_I /= Vp  # dimensionless

    return Lim_I
