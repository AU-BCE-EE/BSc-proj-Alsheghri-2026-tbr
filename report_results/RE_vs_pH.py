'''
What i want here i to plot Removal efficiency vs pH both experimental and model values
For the experiments i only have 2 values, one at pH 13 and the other at pH 12.5 
I will use the results with Ql = 505 mL/min so we do not have the wetting efficiency problem
'''

# =========================== Experimental data from 09/04/2026 =======================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)

# import the data
data_13 = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\09_04_2026_co2_measurements.csv",
    sep=';'
)


data_13['Timestamp'] = pd.to_datetime(data_13['Timestamp'])

# start and end for the inlet measurment
start_13 = "2026-04-09 10:16:41"
end_13   = "2026-04-09 10:21:41"

# extract all the values in the interval
inlet_13 = data_13[(data_13["Timestamp"] >= start_13) & (data_13["Timestamp"] <= end_13)]

# mean inlet conc.
inlet_conc_no_cal_13 = inlet_13['SCD30_CO2'].mean() # this is the diluted value, the real value is calculated using the dilution factor
inlet_conc_13 = (inlet_conc_no_cal_13 - 33.475)/0.9576

Q_gas_bund_13 = 0.387  # L/min
Q_air_mix_13 = 3.0705     # L/min

# Before the dilution
inlet_conc_actual_13 = float((Q_gas_bund_13+Q_air_mix_13) / Q_gas_bund_13 * inlet_conc_13) # ppm



# oulet concentrations
start_out_13 = "2026-04-09 10:27:41"
end_out_13 = "2026-04-09 12:08:41"
outlet_13 = data_13[(data_13["Timestamp"] >= start_out_13) & (data_13["Timestamp"] <= end_out_13)]

Q_gas_top_13 = 0.398  # L/min
Q_air_mix_13 = 3.0705   # L/min

# actual outlet conc.
outlet_conc_no_cal_13 = outlet_13['SCD30_CO2'] * (Q_gas_top_13 + Q_air_mix_13)/Q_gas_top_13 # ppm
outlet_conc_13 = (outlet_conc_no_cal_13- 33.475)/0.9576 
temp   = 22.0           # Celcius
pres   = 1.0            # bar
M_co2 = 44.009          # g/mol
R     = 8.314e-5        # m3 * bar / K-mol

outlet_conc_gm3_13 = pres/(R*(temp+273.15)) * outlet_conc_13/10**6 * M_co2


outlet_times_13 = outlet_13["Timestamp"]
t0_13 = outlet_times_13.iloc[0]
outlet_t_sec_13 = (outlet_times_13 - t0_13).dt.total_seconds().to_numpy()


removal_efficiency_experimental_13 = (inlet_conc_actual_13 - outlet_conc_13)/inlet_conc_actual_13 * 100
# we take an average 
average_RE_13 = np.mean(removal_efficiency_experimental_13)


# =========================== Experimental data from 17/04/2026 =======================


data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\17_04_2026_co2_measurements.csv",
    sep=';'
)


data['Timestamp'] = pd.to_datetime(data['Timestamp'])

# start and end for the inlet measurment
start = "2026-04-17 09:40:38"
end   = "2026-04-17 09:54:38"

# extract all the values in the interval
inlet = data[(data["Timestamp"] >= start) & (data["Timestamp"] <= end)]

# mean inlet conc.
inlet_conc_no_cal = inlet['SCD30_CO2'].mean() # this is the diluted value, the real value is calculated using the dilution factor
inlet_conc = (inlet_conc_no_cal - 33.475)/0.9576 # here i used the devlabs calibration curve

Q_gas_bund = 0.389  # L/min
Q_air_mix = 2.736     # L/min

# Before the dilution
inlet_conc_actual = float((Q_gas_bund+Q_air_mix) / Q_gas_bund * inlet_conc) # ppm



# oulet concentrations
start_out = "2026-04-17 09:57:38"
end_out = "2026-04-17 11:25:38"
outlet = data[(data["Timestamp"] >= start_out) & (data["Timestamp"] <= end_out)]

Q_gas_top = 0.396  # L/min
Q_air_mix = 2.736   # L/min

# actual outlet conc.
outlet_conc_no_cal = outlet['SCD30_CO2'] * (Q_gas_top + Q_air_mix)/Q_gas_top # ppm
outlet_conc = (outlet_conc_no_cal- 33.475)/0.9576  # here i used the devlabs calibration curve
temp   = 22.0           # Celcius
pres   = 1.0            # bar
M_co2 = 44.009          # g/mol
R     = 8.314e-5        # m3 * bar / K-mol

outlet_conc_gm3 = pres/(R*(temp+273.15)) * outlet_conc/10**6 * M_co2


outlet_times = outlet["Timestamp"]
t0 = outlet_times.iloc[0]
outlet_t_sec = (outlet_times - t0).dt.total_seconds().to_numpy()

# experimental removal efficiency
removal_efficiency_experimental = (inlet_conc_actual - outlet_conc)/inlet_conc_actual * 100
# average RE
average_RE_12 = np.mean(removal_efficiency_experimental)

# now we want to run the model for different pH values 

def modelrun(Q_g = 10,
             Q_l = 350,
             D   = 0.19,
             pH  = 13.7,
             vres = 20,
             times = outlet_t_sec,
             frac_co2 = 0.05,
             wet_eff = 1.0,
             Kga = 'onda',
             counter=True,
             recirc=True,
             enh_method='PFO',
             constant_res_pH = True):
    """
    A function for model running
    units of arguments:
    Q_g             L/min
    Q_l             mL/min
    """
    #========= fixed parameters ==========# 
    L     = 0.3       # m
    por_g = 0.86
    por_l = 0.05
    ssa   = 260     # m2/m3
    nc    = 30

    cg0 = 28
    cl_co20   = 0.0
    cl_TOTC0  = 0.0
    clin_co2  = 0.0
    clin_TOTC = 0.0
    cr_co20   = 0.0
    cr_TOTC0  = 0.0



    henry = [3.3e-2, 2400]

    temp   = 22.0    # Celcius
    dens_l = 997.0   # kg/m3
    pres   = 1.0     # bar


    M_co2 = 44.009          # g/mol
    R     = 8.314e-5        # m3 * bar / K-mol

    # ================= Derived values =================

    ex_oh = 10**(-(14-pH)) * 1000 # mol/m3

    # cross sectional area of the reactor
    A = (np.pi*D**2)/4   # m2
    # flow velocity using flow rate and area 
    v_g = (Q_g * 1/1000 * 1/60) / A  # L/min * 1m3/1000L * 1min/60s = m3/s
    v_l = (Q_l * 1/10**6 * 1/60) / A # mL/min * 1m3/10^6mL * 1min/60s = m3/s 

    wet_eff = wet_eff
    # Kga = 0.02286
    Kga = Kga 
    v_res = vres/1000/A
    cgin = pres / ((temp + 273.15) * R) * frac_co2 * M_co2 # g/m3
    # cgin = cgin_df
    results = md.tfmod(
        L, por_g, por_l, v_g, v_l, nc,
        cg0, cl_co20, cl_TOTC0,
        cgin, ex_oh,
        clin_co2, clin_TOTC,
        cr_co20, cr_TOTC0,
        Kga, henry, temp, dens_l,
        times,
        kg='onda',
        kl='onda',
        ae='onda',
        v_res=v_res,
        pres=pres,
        ssa=ssa,
        wet_eff = wet_eff,
        typ='PR',
        counter=counter,
        recirc=recirc,
        enh_method=enh_method,
        constant_res_pH = constant_res_pH
    )

    return results,cgin

pH_span = np.array([12.0, 12.2, 12.4, 12.5, 12.6, 12.8, 12.9, 13.01, 13.2, 13.4, 13.6, 13.8, 14])
pH_span = pH_span - 0.096
removal_list =[]
for pH in pH_span:
    frac_co2 = 0.023671337124733516
    results,cgin = modelrun(Q_g = 10.84, Q_l = 505.4,
                         pH = pH, times = outlet_t_sec,frac_co2 = frac_co2,constant_res_pH=True,
                         enh_method='PFO' if pH>12.6 else 'DC', wet_eff = 1, Kga = 'onda',recirc=True
                         )
    gas = results ['gas_conc']
    gas_out = gas[-1,:]
    gas_in = float(cgin)

    removal = 100 * (gas_in - gas_out) / gas_in
    removal_list.append(removal[-1])

pH_span_plot = pH_span +0.096
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Garamond']
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 11
mpl.rcParams['ytick.labelsize'] = 11
mpl.rcParams['legend.fontsize'] = 11

plt.figure(figsize=(7,5))
plt.title('pH effect on Removal Efficiency')
plt.plot(pH_span_plot, removal_list, c = 'black', marker = 'o', linestyle = 'none', label = 'Model', fillstyle = 'none')
plt.plot(13.01, average_RE_13, c = 'blue', marker = '^', linestyle = 'none', label = 'Experimental', fillstyle = 'none')
plt.plot(12.5, average_RE_12,c = 'blue', marker = '^', linestyle = 'none', fillstyle = 'none')
plt.grid(True, linestyle = '--', linewidth = 0.5, alpha = 0.6)
plt.legend()
plt.ylabel('Removal efficiency [%]')
plt.xlabel('pH value')
plt.show()