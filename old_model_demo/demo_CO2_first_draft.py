import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from mods.mod_co2_first_draft import tfmod 


data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\09_04_2026_co2_measurements.csv",
    sep=';'
)

data['Timestamp'] = pd.to_datetime(data['Timestamp'])

# Inlet processing
start, end = "2026-04-09 10:16:41", "2026-04-09 10:21:41"

inlet = data[(data["Timestamp"] >= start) & (data["Timestamp"] <= end)]

inlet_conc_no_cal = inlet['SCD30_CO2'].mean()
inlet_conc = (inlet_conc_no_cal - 33.475)/0.9576
Q_gas_bund, Q_air_mix_in = 0.387, 3.0705

inlet_conc_actual = float((Q_gas_bund + Q_air_mix_in) / Q_gas_bund * inlet_conc)

# Outlet processing
start_out, end_out = "2026-04-09 10:27:41", "2026-04-09 12:08:41"

outlet = data[(data["Timestamp"] >= start_out) & (data["Timestamp"] <= end_out)]

Q_gas_top, Q_air_mix_out = 0.398, 3.0705

outlet_conc_no_cal = outlet['SCD30_CO2'] * (Q_gas_top + Q_air_mix_out)/Q_gas_top
outlet_conc = (outlet_conc_no_cal - 33.475)/0.9576 

temp, pres, M_co2, R_gas = 22.0, 1.0, 44.009, 8.314e-5

outlet_conc_gm3 = pres/(R_gas*(temp+273.15)) * outlet_conc/10**6 * M_co2

outlet_times = outlet["Timestamp"]

t0 = outlet_times.iloc[0]

outlet_t_sec = (outlet_times - t0).dt.total_seconds().to_numpy()

# pH data processing
time_vals = ["11:00", "11:05", "11:10", "11:15", "11:20", "11:25", "11:30", "11:35", "11:40", "11:45", "11:50", "11:55", "12:00"]
pH_vals = [12.88, 12.88, 12.88, 12.91, 12.90, 12.89, 12.89, 12.90, 12.90, 12.89, 12.89, 12.89, 12.89]
pH_data = pd.DataFrame({"time": time_vals, "pH": pH_vals})
pH_data["Timestamp"] = pd.to_datetime("2026-04-09 " + pH_data["time"])
pH_data["t_sec"] = (pH_data["Timestamp"] - t0).dt.total_seconds()

removal_efficiency_experimental = (inlet_conc_actual - outlet_conc)/inlet_conc_actual * 100

# MODEL PARAMETERS
D_column = 0.19
Area = (np.pi * D_column**2) / 4
Q_g_lmin, Q_l_mlmin = 10.84, 505.4
v_g = (Q_g_lmin * 1e-3 / 60) / Area
v_l = (Q_l_mlmin * 1e-6 / 60) / Area
cgin = pres / ((temp + 273.15) * R_gas) * (inlet_conc_actual / 1e6) * M_co2

L, por_g, por_l, ssa, wet_eff = 0.3, 0.91, 0.05, 260, 0.5
ssa_eff = ssa * wet_eff
TK = temp + 273.15
k_molar = 4.315 * 10**13 * np.exp(-6666/TK) 
k_mass = k_molar * (1/M_co2) * (1/1000) 

# 3. MODEL EXECUTION
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

# Stoichiometry: 2 OH- consumed per 1 CO2 absorbed
moles_oh_consumed_sec = 2 * moles_co2_absorbed_sec

# OH- concentration in inlet (based on pH_lab = 13.01)
pH_inlet = 13.01
oh_inlet_molar = 10**(-(14 - pH_inlet)) # mol/L
moles_oh_in_sec = oh_inlet_molar * (Q_l_mlmin / 1000 / 60) # mol/s

# OH- concentration in outlet
moles_oh_out_sec = moles_oh_in_sec - moles_oh_consumed_sec
oh_outlet_molar = moles_oh_out_sec / (Q_l_mlmin / 1000 / 60) # mol/L

# Convert molarity to pH
# Use clip to avoid log of negative numbers if OH is exhausted
model_pH_outlet = 14 + np.log10(np.clip(oh_outlet_molar, 1e-10, None))

# ==========================================
# 5. PLOTTING
# ==========================================
plt.figure(figsize=(15, 5))

# Plot 1: Gas Outlet
plt.subplot(1, 3, 1)
plt.title('Gas Outlet Concentration')
plt.plot(outlet_t_sec, model_gas_out, 'r-', label='Simple Model')
plt.plot(outlet_t_sec, outlet_conc_gm3, 'bo', label='Lab Data', markersize=2)
plt.ylabel('CO2 [g/m3]')
plt.legend()
plt.grid(True)

# Plot 2: Removal Efficiency
plt.subplot(1, 3, 2)
plt.title('Removal Efficiency')
eff = 100 * (cgin - model_gas_out) / cgin
plt.plot(outlet_t_sec, eff, 'r-', label='Simple Model')
plt.plot(outlet_t_sec, removal_efficiency_experimental, 'bo', label='Lab Data', markersize=2)
plt.ylim(0, 100)
plt.ylabel('% Removal')
plt.grid(True)

# Plot 3: pH Outlet
plt.subplot(1, 3, 3)
plt.title('Outlet pH')
plt.plot(outlet_t_sec, model_pH_outlet, 'r-', label='Model (Calculated)')
plt.plot(pH_data["t_sec"], pH_data["pH"], 'bo', label='Lab Data', markersize=4)
plt.ylabel('pH')
plt.xlabel('Time [s]')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()