import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tomli

from npzd_funcs import *

class NPZDModel:
    def __init__(self, config_file):
        self.load_config(config_file)
        self.initialize_model()

    def load_config(self, config_file):
        with open(config_file, "rb") as f:
            config = tomli.load(f)
        self.params = config["PARAMETER"]
        self.settings = config["SETUP"]

    def initialize_model(self):
        self.set_station_settings()
        self.interpolate_forcing_data()

    def set_station_settings(self):
        station = self.settings["flag_stn"]
        # Station settings
        forcing = pd.read_csv("stations_forcing.txt")

        ## station-specific settings
        station_settings = {
            1: {"lat": 60.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0074, "bN": 10.85, "MLD": forcing.iloc[:13, 0], "SST": forcing.iloc[:13, 1]},
            2: {"lat": 47.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0174, "bN": 4.0, "MLD": forcing.iloc[:13, 2], "SST": forcing.iloc[:13, 3]},
            3: {"lat": -50.67, "clouds": 6.0, "e0": 12.0, "aN": 0, "bN": 26.1, "MLD": forcing.iloc[:13, 6], "SST": forcing.iloc[:13, 7]},
            4: {"lat": 50.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0, "bN": 14.6, "MLD": forcing.iloc[:13, 4], "SST": forcing.iloc[:13, 5]}
        }
        self.station_settings = station_settings[station]
        self.lat_radian = self.station_settings["lat"]*np.pi/180

    def interpolate_forcing_data(self):
        self.MLD = interpolate_monthly_to_daily(self.station_settings["MLD"])
        self.SST = interpolate_monthly_to_daily(self.station_settings["SST"])
    
    def nutrient_limitation(self, N, kN):
        return N / (kN + N)

    def temp_limitation(self, T):
        return 1.066**T

    def get_deep_water_nitrate(self, depth, aN, bN):
        return aN * depth + bN        

    def calculate_time_step(self, current_step, n_day):
        n_pday= int(1/self.settings["delta_t"])
        day_number = (current_step // n_pday) % n_day + 1
        step_in_day = current_step % n_pday

        return day_number, step_in_day

    def ode(self, X, t):
        P, N, Z, D = X

        # Calculate day and hour based on t
        iday = np.floor(t).astype(int) + 1  # 1-indexed days
        ihour = np.floor((t - iday + 1) * 24).astype(int)  # 0-indexed hours

        day_in_year = int(t) % 365
        hour_in_day = ihour % 24
        
        MLD1 = self.MLD[day_in_year]        
        MLD2 = self.MLD[day_in_year + 1]

        if self.settings["flag_force"]:
            MLDnow = self.settings['force_MLD']
            N0 = self.settings['force_N0']
            SSTnow = self.settings['force_SST']
            noonparnow = self.settings['force_irradiance']
        else:            
            MLDnow = MLD1 + (MLD2 - MLD1) * hour_in_day / 24
            N0 = self.get_deep_water_nitrate(MLDnow, self.station_settings["aN"], self.station_settings["bN"])
            SSTnow = self.SST[day_in_year]        
            noonparnow = noon_par(day_in_year, self.lat_radian, self.station_settings["clouds"], self.station_settings["e0"])            
        
        VpT = self.params["Vp0"] * self.temp_limitation(SSTnow)            
                

        if self.settings["flag_force"]:
            daylnow = 12
        else:
            daylnow = day_length(day_in_year, self.lat_radian)
            


        kPAR, zbase, Ibase, zdep, aphybase = [np.zeros(4) for _ in range(5)]
        chl = P * 6.625 * 12.0 / self.params["CtoChl"] ## used for light attenuation
        ss = np.sqrt(chl)
        
        ## get the light attenuation coefficient
        match self.settings['flag_atten']:
            case 1:
                jnlay = 1
                zbase[0] = 0.0
                zbase[1] = MLDnow
                zdep[1] = MLDnow
                kPAR[1] = self.params['kw'] + self.params['kc']*P
                Ibase[0] = noonparnow
                Ibase[1] = Ibase[0] * np.exp(-kPAR[1] * MLDnow)

                
            case 2:
                jnlay = 2
                Ibase[0] = noonparnow
                kPAR[1] = 0.13096 + 0.030969*ss + 0.042644*ss**2 - 0.013738*ss**3 + 0.0024617*ss**4 - 0.00018059*ss**5
                kPAR[2] =  0.041025 + 0.036211*ss + 0.062297*ss**2 - 0.030098*ss**3 + 0.0062597*ss**4 - 0.00051944*ss**5
                kPAR[3] = 0.021517 + 0.050150*ss + 0.058900*ss**2 - 0.040539*ss**3 + 0.0087586*ss**4 - 0.00049476*ss**5
                zbase[0] = 0.0
                Ibase[0] = noonparnow

                if MLDnow <= 5:
                    jnlay = 1
                    zbase[1] = MLDnow
                    zdep[1] = MLDnow
                    Ibase[1] = Ibase[0] * np.exp(-kPAR[1] * zdep[1])
                elif MLDnow <= 23 and MLDnow > 5:
                    jnlay = 2
                    zbase[1] = 5.0
                    zdep[1] = 5.0
                    Ibase[1] = Ibase[0] * np.exp(-kPAR[1] * 5.0)
                    zbase[2] = MLDnow
                    zdep[2] = MLDnow - 5.0
                    Ibase[2] = Ibase[1] * np.exp(-kPAR[2] * zdep[2])
                else:
                    jnlay = 3
                    zbase[1] = 5.0
                    zdep[1] = 5.0
                    Ibase[1] = Ibase[0] * np.exp(-kPAR[1] * 5.0)
                    zbase[2] = 23.0
                    zdep[2] = 18.0
                    Ibase[2] = Ibase[1] * np.exp(-kPAR[2] * zdep[2])
                    zbase[3] = MLDnow
                    zdep[3] = MLDnow - 23.0
                    Ibase[3] = Ibase[2] * np.exp(-kPAR[3] * zdep[3])
                        

        ## calculate photosynthesis
        L_Isum = 0.0
        for ilay in range(1,jnlay+1):
            match self.settings["flag_LI"]:
                case 1:
                      L_I = calc_li_num(zdep[ilay], Ibase[ilay-1], Ibase[ilay], kPAR[ilay],
                                        self.params["alpha"],VpT,daylnow,
                                        self.settings["flag_irrad"],
                                        self.settings["flag_PIcurve"])
                case 2:
                    L_I = calc_li_ep85(zdep[ilay], Ibase[ilay-1], Ibase[ilay], kPAR[ilay], self.params["alpha"], VpT, daylnow)
                case 3:
                    aphybase[0] = 0.36796 + 0.17537*ss - 0.065276*ss**2 + 0.013528*ss**3 - 0.0011108*ss**4
                    ahash = calc_aphy(ss, zbase[ilay-1],zbase[ilay], aphybase[ilay-1])
                    aphybase[ilay] = aphybase[ilay-1]+ahash
                    aphyav = aphybase[ilay-1]+ahash*0.5
                    L_I = calc_li_A93(zdep[ilay],Ibase[ilay-1],Ibase[ilay],kPAR[ilay],self.params['alpha'],VpT,daylnow,aphyav)
                    
            L_Isum += L_I * zdep[ilay]         
                                
        L_I_avg = L_Isum/MLDnow                   # weighted average for mixed layer

        L_N = self.nutrient_limitation(N, self.params["kN"])

        if self.settings["flag_grazing"] == 1:
            phiPhat, phiDhat = self.params["phiP"], self.params["phiD"]
        else:
            phiPhat, phiDhat = self.params["phiP"] * P, self.params["phiD"] * D
        
        intakespP = self.params["Imax"] * phiPhat * P / (self.params["kz"]**2 + phiPhat * P + phiDhat * D)
        intakespD = self.params["Imax"] * phiDhat * D / (self.params["kz"]**2 + phiPhat * P + phiDhat * D)

        Pgrow = VpT*L_I_avg*L_N* 24/self.params['CtoChl']*P

        Pgraz = intakespP * Z
        Dgraz = intakespD * Z
        Pmort = self.params["mP"] * P
        Pmort2 = self.params["mP2"] * P * P
        Pmix = self.params["wmix"] * P / MLDnow
        Zmix = self.params["wmix"] * Z / MLDnow


        if self.settings["flag_force"]:
            Nmix = self.params["wmix"] * (N0 - N) #rapid mixing for chemostat
        else:
           Nmix = self.params["wmix"] * (N0 - N) / MLDnow

        Zgrow = self.params["betaz"] * self.params["kNz"] * (Pgraz + Dgraz)
        Dremin = self.params["mD"] * D
               
        
        Zexc = self.params["betaz"] * (1.0 - self.params["kNz"]) * (Pgraz + Dgraz)
        Zpel = (1.0 - self.params["betaz"]) * (Pgraz + Dgraz)
        Zmort = self.params["mz"] * Z
        Zmort2 = self.params["mz2"] * Z**2
        Dmix = self.params["wmix"] * D / MLDnow
        Dsink = self.params["VD"] * D / MLDnow
        

        diffMLD = MLD2 - MLD1
        if (diffMLD > 0.0):
            Nadd = diffMLD * N0 / MLDnow
            Ndilute = N * diffMLD / MLDnow
            Pdilute = P * diffMLD / MLDnow
            Zdilute = Z * diffMLD / MLDnow
            Ddilute = D * diffMLD / MLDnow
        else:
            Nadd = 0.0
            Ndilute = 0.0
            Pdilute = 0.0
            Zdilute = 0.0
            Ddilute = 0.0                
               

        dPdt = Pgrow - Pgraz - Pmort - Pmort2 - Pmix - Pdilute
        dNdt = Nmix + Zexc + Dremin - Pgrow + Nadd - Ndilute
        dZdt = Zgrow - Zmort - Zmort2 - Zmix - Zdilute
        dDdt = Pmort + Pmort2 + Zmort + Zpel - Dgraz - Dmix - Dsink - Dremin - Ddilute

        return [dPdt, dNdt, dZdt, dDdt]

    def run(self):
        ndays = int(self.settings["nyears"] * 365)
        nhr_pday = int(24/self.settings["delta_t"])
        t = np.linspace(0, ndays, ndays * nhr_pday)
        X0 = [self.settings["Pinit"], self.settings["Ninit"], 
              self.settings["Zinit"], self.settings["Dsinit"]]


        print("*************************************")
        print("         Model running")
        print("*************************************")
        
        solution = odeint(self.ode, X0, t)

        print("*************************************")
        print("         Model completed")
        print("*************************************")        
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
    model = NPZDModel("config.toml")
    t, solution = model.run()
    ## cut the first year
    t = t[365:]
    solution = solution[365:]    
    model.plot_results(t, solution)
    

if __name__ == "__main__":
    main()
