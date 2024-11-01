from npzd_funcs import calc_li_num, calc_li_ep85, calc_li_A93,calc_aphy
from npzd import noon_par

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
    result = round(calc_aphy(10,20,30,10),6)
    assert result == 3.902094

# def test_cal_day():
#     "original: FNdaylcalc"
#     result = calc_daylenth(10,0.5)
#     assert result == 10.29673

def test_noon_par():
    "original: FNnoonparcalc"
    result = round(noon_par(300, 0.1,0,100),4)
    assert result == 420.9175

