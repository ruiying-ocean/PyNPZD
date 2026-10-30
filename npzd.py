import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tomli

from forcing import RemoteForcingLoader
from npzd_funcs import day_length, interpolate_monthly_to_daily, noon_par, calc_li_num, calc_li_ep85, calc_li_A93, calc_aphy

class NPZDModel:
    def __init__(self, config_file):
        self.load_config(config_file)
        self.initialize_model()

    def load_config(self, config_file):
        with open(config_file, "rb") as f:
            config = tomli.load(f)
        self.params = config["PARAMETER"]
        self.settings = config["SETUP"]
        self.forcing_config = config.get("FORCING", {})

    def initialize_model(self):
        self.remote_noon_par_monthly = None
        self.remote_nutrient_monthly = None
        self.remote_nutrient_cycle = None
        self.set_station_settings()
        self.interpolate_forcing_data()
        self.precompute_daily_cycles()

    def set_station_settings(self):
        forcing_source = self.forcing_config.get("source", "station").lower()
        if forcing_source == "remote":
            self._set_remote_station_settings()
        else:
            self._set_local_station_settings()

        self.lat_radian = self.station_settings["lat"] * np.pi / 180

    def _set_local_station_settings(self):
        station = self.settings["flag_stn"]
        forcing = pd.read_csv("stations_forcing.txt")

        station_settings = {
            1: {"lat": 60.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0074, "bN": 10.85, "MLD": forcing.iloc[:13, 0], "SST": forcing.iloc[:13, 1]},
            2: {"lat": 47.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0174, "bN": 4.0, "MLD": forcing.iloc[:13, 2], "SST": forcing.iloc[:13, 3]},
            3: {"lat": -50.67, "clouds": 6.0, "e0": 12.0, "aN": 0, "bN": 26.1, "MLD": forcing.iloc[:13, 6], "SST": forcing.iloc[:13, 7]},
            4: {"lat": 50.0, "clouds": 6.0, "e0": 12.0, "aN": 0.0, "bN": 14.6, "MLD": forcing.iloc[:13, 4], "SST": forcing.iloc[:13, 5]},
        }
        self.station_settings = station_settings[station]

        # Local forcing always uses climatology
        self.time_coords = None
        self.use_time_series = False

    def _set_remote_station_settings(self):
        remote_cfg = self.forcing_config.get("remote")
        if remote_cfg is None:
            raise ValueError("FORCING.remote section is required when source='remote'")

        loader = RemoteForcingLoader(remote_cfg)
        result = loader.load()

        # Handle both climatology and time series data
        if isinstance(result, tuple):
            # Time series mode: (data_dict, time_coords)
            remote_data, self.time_coords = result
            self.use_time_series = True
        else:
            # Climatology mode: just data_dict
            remote_data = result
            self.time_coords = None
            self.use_time_series = False

        if "MLD" not in remote_data or remote_data["MLD"] is None:
            raise ValueError("Remote forcing must provide an 'MLD' dataset")
        if "SST" not in remote_data or remote_data["SST"] is None:
            raise ValueError("Remote forcing must provide an 'SST' dataset")

        self.station_settings = {
            "lat": remote_cfg["latitude"],
            "clouds": remote_cfg.get("clouds", 6.0),
            "e0": remote_cfg.get("e0", 12.0),
            "aN": remote_cfg.get("aN", 0.0),
            "bN": remote_cfg.get("bN", 0.0),
            "MLD": remote_data.get("MLD"),
            "SST": remote_data.get("SST"),
        }

        irradiance = remote_data.get("IRRADIANCE")
        if irradiance is not None:
            self.remote_noon_par_monthly = irradiance

        nutrient = remote_data.get("NUTRIENT")
        if nutrient is not None:
            # Use remote nutrient data directly - no normalization needed
            # The real oceanographic data already represents realistic mixed layer nutrients
            self.remote_nutrient_monthly = np.asarray(nutrient, dtype=float)
            if self.remote_nutrient_monthly.size:
                # Override the baseline nutrient concentration with the remote mean
                self.station_settings["bN"] = float(np.nanmean(self.remote_nutrient_monthly))
        else:
            self.remote_nutrient_monthly = None

    def interpolate_forcing_data(self):
        if self.use_time_series:
            if self.time_coords is not None and np.any(
                np.diff(self.time_coords) != np.timedelta64(1, "D")
            ):
                raise ValueError("Time-series forcing must have daily samples")
            # For time series data, use directly without interpolation
            self.MLD = self.station_settings["MLD"]
            self.SST = self.station_settings["SST"]
            if self.remote_nutrient_monthly is not None:
                self.remote_nutrient_cycle = self.remote_nutrient_monthly
            else:
                self.remote_nutrient_cycle = None
        else:
            # For climatology data, interpolate to daily values
            self.MLD = interpolate_monthly_to_daily(self.station_settings["MLD"])
            self.SST = interpolate_monthly_to_daily(self.station_settings["SST"])
            if self.remote_nutrient_monthly is not None:
                self.remote_nutrient_cycle = interpolate_monthly_to_daily(self.remote_nutrient_monthly)
            else:
                self.remote_nutrient_cycle = None
    
    def precompute_daily_cycles(self):
        if self.settings["flag_force"]:
            self.day_length_cycle = None
            self.noon_par_cycle = None
            return

        if self.use_time_series:
            # For time series, we'll compute light values dynamically in the ODE
            self.day_length_cycle = None
            if self.remote_noon_par_monthly is not None:
                self.noon_par_cycle = self.remote_noon_par_monthly  # Use directly
            else:
                self.noon_par_cycle = None
        else:
            # For climatology, precompute annual cycles
            days = np.arange(1, 366)
            self.day_length_cycle = day_length(days, self.lat_radian)
            if self.remote_noon_par_monthly is not None:
                self.noon_par_cycle = interpolate_monthly_to_daily(self.remote_noon_par_monthly)
            else:
                clouds = self.station_settings["clouds"]
                e0 = self.station_settings["e0"]
                self.noon_par_cycle = noon_par(days, self.lat_radian, clouds, e0)
    
    def nutrient_limitation(self, N, kN):
        return N / (kN + N)

    def temp_limitation(self, T):
        return 1.066**T

    def get_deep_water_nitrate(self, depth, aN, bN, use_remote=False, remote_value=None):
        """Calculate deep water nitrate concentration.

        Args:
            depth: Mixed layer depth (m)
            aN: Linear coefficient for depth relationship
            bN: Baseline nutrient concentration
            use_remote: If True, use remote_value directly
            remote_value: Remote nutrient data value

        Returns:
            N0: Deep water nutrient concentration (mmol N m-3)
        """
        if use_remote and remote_value is not None:
            # Use remote data directly - it already represents mixed layer nutrients
            return remote_value
        else:
            # Fall back to linear parameterization
            return aN * depth + bN        

    def calculate_time_step(self, current_step, n_day):
        n_pday= int(1/self.settings["delta_t"])
        day_number = (current_step // n_pday) % n_day + 1
        step_in_day = current_step % n_pday

        return day_number, step_in_day

    def ode(self, X, t):
        P, N, Z, D = X

        if self.settings["flag_force"]:
            MLDnow = self.settings['force_MLD']
            N0 = self.settings['force_N0']
            SSTnow = self.settings['force_SST']
            noonparnow = self.settings['force_irradiance']
            daylnow = 12
        else:
            if self.use_time_series:
                # For time series data, use time directly as index
                time_index = min(int(t), len(self.MLD) - 1)

                MLD1 = self.MLD[time_index]
                MLD2 = self.MLD[min(time_index + 1, len(self.MLD) - 1)]
                MLDnow = MLD1 + (MLD2 - MLD1) * (t - int(t))

                if self.remote_nutrient_cycle is not None:
                    remote_nutrient_today = self.remote_nutrient_cycle[time_index]
                    N0 = self.get_deep_water_nitrate(MLDnow, self.station_settings["aN"], self.station_settings["bN"],
                                                   use_remote=True, remote_value=remote_nutrient_today)
                else:
                    N0 = self.get_deep_water_nitrate(MLDnow, self.station_settings["aN"], self.station_settings["bN"])
                SSTnow = self.SST[time_index]

                # For time series, compute day length and PAR based on time
                if self.time_coords is not None and len(self.time_coords) > 0:
                    # Convert time coordinate to day of year for light calculations
                    start_time = pd.to_datetime(self.time_coords[0])
                    current_time = start_time + pd.Timedelta(days=t)
                    day_of_year = current_time.timetuple().tm_yday
                else:
                    day_of_year = int(t) % 365 + 1

                daylnow = day_length(np.array([day_of_year]), self.lat_radian)[0]

                if self.remote_noon_par_monthly is not None and time_index < len(self.remote_noon_par_monthly):
                    noonparnow = self.remote_noon_par_monthly[time_index]
                else:
                    clouds = self.station_settings["clouds"]
                    e0 = self.station_settings["e0"]
                    noonparnow = noon_par(np.array([day_of_year]), self.lat_radian, clouds, e0)[0]
            else:
                # Original climatology logic
                day_index = int(t)
                day_in_year = day_index % 365
                hour_in_day = int((t - day_index) * 24)

                MLD1 = self.MLD[day_in_year]
                MLD2 = self.MLD[day_in_year + 1]
                MLDnow = MLD1 + (MLD2 - MLD1) * hour_in_day / 24
                if self.remote_nutrient_cycle is not None:
                    remote_nutrient_today = self.remote_nutrient_cycle[day_in_year]
                    N0 = self.get_deep_water_nitrate(MLDnow, self.station_settings["aN"], self.station_settings["bN"],
                                                   use_remote=True, remote_value=remote_nutrient_today)
                else:
                    N0 = self.get_deep_water_nitrate(MLDnow, self.station_settings["aN"], self.station_settings["bN"])
                SSTnow = self.SST[day_in_year]
                noonparnow = self.noon_par_cycle[day_in_year]
                daylnow = self.day_length_cycle[day_in_year]
        
        VpT = self.params["Vp0"] * self.temp_limitation(SSTnow)                        

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
            Nmix = self.params["wmix"] * (N0 - N) # rapid mixing for chemostat
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
        

        if self.settings['flag_force']:
            diffMLD = 0
        else:
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
        if self.use_time_series:
            # For time series, use the length of the forcing data
            ndays = len(self.MLD)
            print(f"Using time series forcing data for {ndays} days ({ndays/365:.1f} years)")
        else:
            # For climatology, use the nyears setting
            ndays = int(self.settings["nyears"] * 365)
            print(f"Using climatology forcing for {ndays} days ({ndays/365:.1f} years)")

        nhr_pday = int(24/self.settings["delta_t"])
        t = np.linspace(0, ndays, ndays * nhr_pday, endpoint=False)
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
