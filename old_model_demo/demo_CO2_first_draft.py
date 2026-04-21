import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mods.mod_co2_first_draft import tfmod 


import numpy as np
import pandas as pd

# import the data
data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\09_04_2026_co2_measurements.csv",
    sep=';'
)


data['Timestamp'] = pd.to_datetime(data['Timestamp'])

# start and end for the inlet measurment
start = "2026-04-09 10:16:41"
end   = "2026-04-09 10:21:41"

# extract all the values in the interval
inlet = data[(data["Timestamp"] >= start) & (data["Timestamp"] <= end)]

# mean inlet conc.
inlet_conc_no_cal = inlet['SCD30_CO2'].mean() # this is the diluted value, the real value is calculated using the dilution factor
inlet_conc = (inlet_conc_no_cal - 33.475)/0.9576

Q_gas_bund = 0.387  # L/min
Q_air_mix = 2.8     # L/min

# Before the dilution
inlet_conc_actual = float((Q_gas_bund+Q_air_mix) / Q_gas_bund * inlet_conc) # ppm

# oulet concentrations
start_out = "2026-04-09 10:27:41"
end_out = "2026-04-09 12:08:41"
outlet = data[(data["Timestamp"] >= start_out) & (data["Timestamp"] <= end_out)]

Q_gas_top = 0.398  # L/min
Q_air_mix = 3.0705   # L/min

# actual outlet conc.
outlet_conc_no_cal = outlet['SCD30_CO2'] * (Q_gas_top + Q_air_mix)/Q_gas_top # ppm
outlet_conc = (outlet_conc_no_cal- 33.475)/0.9576 
temp   = 22.0           # Celcius
pres   = 1.0            # bar
M_co2 = 44.009          # g/mol
R     = 8.314e-5        # m3 * bar / K-mol

outlet_conc_gm3 = pres/(R*(temp+273.15)) * outlet_conc/10**6 * M_co2


outlet_times = outlet["Timestamp"]
t0 = outlet_times.iloc[0]
outlet_t_sec = (outlet_times - t0).dt.total_seconds().to_numpy()

# pH data 
time = ["11:00", "11:05", "11:10", "11:15", "11:20", "11:25", "11:30", 
        "11:35", "11:40", "11:45", "11:50", "11:55", "12:00"]
pH   = [12.88, 12.88, 12.88, 12.91, 12.90, 12.89, 12.89, 12.90, 12.90, 12.89, 12.89, 12.89, 12.89]
pH_data = pd.DataFrame({
    "time": time,
    "pH"  : pH
})


pH_data["Timestamp"] = pd.to_datetime(
    "2026-04-09 " + pH_data["time"]
)


t0 = outlet_times.iloc[0]
pH_data["t_sec"] = (pH_data["Timestamp"] - t0).dt.total_seconds()



# experimental removal efficiency
removal_efficiency_experimental = (inlet_conc_actual - outlet_conc)/inlet_conc_actual * 100

D = 0.19
Area = (np.pi * D**2) / 4
Q_g_lmin = 10.84 
Q_l_mlmin = 505.4
v_g = (Q_g_lmin * 1e-3 / 60) / Area
v_l = (Q_l_mlmin * 1e-6 / 60) / Area
cgin = pres / ((temp + 273.15) * R) * (inlet_conc_actual / 1e6) * M_co2

L, por_g, por_l, ssa, wet_eff = 0.3, 0.91, 0.05, 260, 0.5
ssa_eff = ssa * wet_eff
TK = temp + 273.15
k_molar = 4.315 * 10**13 * np.exp(-6666/TK) 
k_mass = k_molar * (1/M_co2) * (1/1000) 


results = tfmod(
    L=L, por_g=por_g, por_l=por_l, v_g=v_g, v_l=v_l, nc=60,
    cg0=10.7, cl0=0.0, cgin=cgin, clin=0.0, 
    k=k_mass, Kga='onda', henry=[3.3e-2, 2400], pKa=6.35, pH=13.01, 
    temp=temp, dens_l=997.0, times=outlet_t_sec, 
    v_res=(20/1000/Area), pres=pres, ssa=ssa_eff, counter=True, recirc=True
)


model_gas_out = results['gas_conc'][-1, :]

# Molar flow of CO2 absorbed (mol/s) = (cgin - cgout) * FlowRate_Gas
# cgin and model_gas_out are in g/m3, so divide by M_co2 for mol/m3
moles_co2_absorbed_sec = (cgin - model_gas_out) / M_co2 * (v_g * Area)




plt.figure(figsize=(15, 5))

plt.subplot(1, 2, 1)
plt.title('Gas Outlet Concentration')
plt.plot(outlet_t_sec, model_gas_out, 'r-', label='Simple Model')
plt.plot(outlet_t_sec, outlet_conc_gm3, 'bo', label='Lab Data', markersize=2)
plt.ylabel('CO2 [g/m3]')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.title('Removal Efficiency')
eff = 100 * (cgin - model_gas_out) / cgin
plt.plot(outlet_t_sec, eff, 'r-', label='Simple Model')
plt.plot(outlet_t_sec, removal_efficiency_experimental, 'bo', label='Lab Data', markersize=2)
plt.ylim(0, 100)
plt.ylabel('% Removal')
plt.grid(True)


plt.tight_layout()
plt.show()