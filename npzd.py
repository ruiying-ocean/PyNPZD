import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tomli

from photosynthesis import *

class NPZDModel:
    def __init__(self, config_file):
        self.load_config(config_file)
        self.initialize_model()

    def load_config(self, config_file):
        with open(config_file, "rb") as f:
            config = tomli.load(f)
        self.params = config["Parameter"]
        self.settings = config["Setup"]

    def initialize_model(self):
        self.set_station_settings()
        self.interpolate_forcing_data()

    def set_station_settings(self):
        station = self.settings["flag_stn"]
        # Station settings
        forcing = pd.read_csv("stations_forcing.txt")

        station_settings = {
            1: {"lat": 60.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0074, "bN": 10.85, "MLD": forcing.iloc[:13, 0], "SST": forcing.iloc[:13, 1]},
            2: {"lat": 47.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0174, "bN": 4.0, "MLD": forcing.iloc[:13, 2], "SST": forcing.iloc[:13, 3]},
            3: {"lat": -50.67, "clouds": 6.0, "e0": 12.0, "aN": 0, "bN": 26.1, "MLD": forcing.iloc[:13, 6], "SST": forcing.iloc[:13, 7]},
            4: {"lat": 50.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0, "bN": 14.6, "MLD": forcing.iloc[:13, 4], "SST": forcing.iloc[:13, 5]}
        }
        self.station_settings = station_settings[station]
        self.latitude = np.radians(self.station_settings["lat"])

    def interpolate_forcing_data(self):
        self.MLD = self.interpolate_yearly_data(self.station_settings["MLD"])
        self.SST = self.interpolate_yearly_data(self.station_settings["SST"])

    @staticmethod
    def interpolate_yearly_data(monthly_data):
        daily_data = np.zeros(366)
        for i in range(365):
            month = i // 30
            next_month = (month + 1) % 12
            frac = (i % 30) / 30
            daily_data[i] = monthly_data[month] * (1 - frac) + monthly_data[next_month] * frac
        daily_data[365] = daily_data[0]
        return daily_data

    def day_length(self, day):
        """calculate day length in hours"""
        declin = 23.45 * np.sin(2 * np.pi * (284 + day) * 0.00274) * np.pi / 180
        daylnow =  24 * np.arccos(-np.tan(self.latitude) * np.tan(declin)) / np.pi
        return daylnow


    def noon_par(self, day):
        """Calculation of noon irradiance"""
        albedo = 0.04 # albedo
        solarconst = 1368.0 # solar constant, w m-2
        parrac = 0.43 # PAR fraction

        # solar declination angle
        declin = 23.45 * np.sin(2 * np.pi * (284 + day) * 0.00274) * np.pi / 180
        # cosine of zenith angle
        coszen = np.sin(self.latitude) * np.sin(declin) + np.cos(self.latitude) * np.cos(declin)
        zen = np.arccos(coszen)
        # Earth's radius vector
        Rvector = 1 / np.sqrt(1 + 0.033 * np.cos(2 * np.pi * day * 0.00274))
        # irradiance at ocean surface, clear sky
        Iclear = solarconst * coszen**2 / (Rvector**2) / (1.2 * coszen + self.station_settings["e0"] * (1.0 + coszen) * 0.001 + 0.0455)
        # cloud factor (atmospheric transmission)
        cfac = (1 - 0.62 * self.station_settings["clouds"] / 8 + 0.0019 * (np.pi/2 - zen))
        # noon irradiance considering PAR fraction and cloud factor
        Inoon_par = Iclear * cfac * (1 - albedo) * parrac
        
        return Inoon_par

    def light_limitation(self, Iin, Iout, kPAR, alpha, Vp, daylength):
        Iavg = (Iin + Iout) / 2
        return 2 * (1 - np.exp(-alpha * Iavg / Vp)) * (daylength / 24)

    def npzd_ode(self, X, t):
        P, N, Z, D = X
        day = int(t) % 365
        
        MLDnow = self.MLD[day]

        ## deep-water nitrate concentration
        N0 = self.station_settings["aN"] * MLDnow + self.station_settings["bN"]
        VpT = self.params["Vp0"] * 1.066**self.SST[day]
        noonparnow = self.noon_par(day)
        daylnow = self.day_length(day)

        chl = P * 6.625 * 12.0 / self.params["CtoChl"]
        kPAR = self.params["kw"] + self.params["kc"] * P
        Iin = noonparnow
        Iout = Iin * np.exp(-kPAR * MLDnow)

        L_I = self.light_limitation(Iin, Iout, kPAR, self.params["alpha"], VpT, daylnow)
        L_N = N / (self.params["kN"] + N)

        if self.settings["flag_grazing"] == 1:
            phiPhat, phiDhat = self.params["phiP"], self.params["phiD"]
        else:
            phiPhat, phiDhat = self.params["phiP"] * P, self.params["phiD"] * D
        
        intakespP = self.params["Imax"] * phiPhat * P / (self.params["kz"]**2 + phiPhat * P + phiDhat * D)
        intakespD = self.params["Imax"] * phiDhat * D / (self.params["kz"]**2 + phiPhat * P + phiDhat * D)

        Pgrow = VpT * L_I * L_N * P
        Pgraz = intakespP * Z
        Dgraz = intakespD * Z
        Pmort = self.params["mP"] * P
        Pmort2 = self.params["mP2"] * P * P
        Pmix = self.params["wmix"] * P / MLDnow
        Zmix = self.params["wmix"] * Z / MLDnow
        Nmix = self.params["wmix"] * (N0 - N) / MLDnow
        Zgrow = self.params["betaz"] * self.params["kNz"] * (Pgraz + Dgraz)
        Zexc = self.params["betaz"] * (1.0 - self.params["kNz"]) * (Pgraz + Dgraz)
        Zpel = (1.0 - self.params["betaz"]) * (Pgraz + Dgraz)
        Zmort = self.params["mz"] * Z
        Zmort2 = self.params["mz2"] * Z**2
        Dmix = self.params["wmix"] * D / MLDnow
        Dsink = self.params["VD"] * D / MLDnow
        Dremin = self.params["mD"] * D

        dPdt = Pgrow - Pgraz - Pmort - Pmort2 - Pmix
        dNdt = Nmix + Zexc + Dremin - Pgrow
        dZdt = Zgrow - Zmort - Zmort2 - Zmix
        dDdt = Pmort + Pmort2 + Zmort + Zpel - Dgraz - Dmix - Dsink - Dremin

        return [dPdt, dNdt, dZdt, dDdt]

    def run_simulation(self):
        ndays = int(self.settings["nyears"] * 365)
        t = np.linspace(0, ndays, ndays * int(1/self.settings["delta_t"]))
        X0 = [self.settings["Pinit"], self.settings["Ninit"], 
              self.settings["Zinit"], self.settings["Dsinit"]]
        
        solution = odeint(self.npzd_ode, X0, t)
        return t, solution

    def plot_results(self, t, solution):
        fig, axs = plt.subplots(2, 2)
        variables = ['P', 'N', 'Z', 'D']
        for i, (ax, var) in enumerate(zip(axs.flatten(), variables)):
            ax.plot(t/365, solution[:, i])
            ax.set_xlabel('Years')
            ax.set_ylabel(f'{var} concentration')
            ax.set_title(f'{var} over time')
        plt.tight_layout()
        plt.show()

# Usage
def main():
    model = NPZDModel("npzd_config.toml")
    t, solution = model.run_simulation()
    ## cut the first year
    t = t[365:]
    solution = solution[365:]    
    model.plot_results(t, solution)
    

if __name__ == "__main__":
    main()
