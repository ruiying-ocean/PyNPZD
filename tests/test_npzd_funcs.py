from pynpzd.npzd_funcs import *

def test_calc_li_num():
    "original: FNLIcalcNu"
    result = round(calc_li_num(100, 1,1,1,1,1,1,1,1), 10)
    assert result == 0.0001945482

    
def test_li_ep85():
    "original: FNLIcalcEP85"
    result = round(calc_li_ep85(100,30,20,1,1,1,1),9)
    assert result == 0.000162144

def test_li_A93():
    "original: FNLIcalcA93"
    result = round(calc_li_A93(100,30,20,1,1,1,1,1),6)
    assert result == 5.598198


def test_aphy():
    "original: FNaphy"
    result = round(calc_aphy(10,10,20,1),6)
    assert result == 3.348472

def test_noon_par():
    "original: FNnoonparcalc"
    result = round(noon_par(300, 0.1,0,1),4)
    assert result == 489.6077


def test_daylength():
    "original: FNdaylcalc"
    result = round(day_length(10,0.6),6)
    assert result == 9.856603
    
def test_interp_monthly():
    "original: FNinterpolate"
    monthly_mld = np.array([242.4, 510.5, 547.3 ,362.7 , 65.2 , 40.2 , 
        22.1 , 36.2 , 52.7 , 68.0  ,88.0 ,144.9, 242.4])
    daily_mld= interpolate_monthly_to_daily(monthly_mld)
    result = round(daily_mld[0], 8)
    assert result==194.31780822