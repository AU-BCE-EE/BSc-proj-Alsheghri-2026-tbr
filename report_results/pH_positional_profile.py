"""
Model validation with experimental data
The pH in the reservoir was 
"""


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
Q_air_mix = 3.0705     # L/min

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



import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)

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
    nc    = 60

    cg0 = 9.9
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

    TK_run  = temp + 273.15
    KW_run  = 10**(-4.2195 - 2915.16 / TK_run)   # same formula as tfmod
    c_h_run = 10**(-pH)                           # mol/L
    # Charge balance for NaOH solution: [Na+] = [OH-] - [H+]
    ex_oh   = (KW_run / c_h_run - c_h_run) * 1000  # mol/m3

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

    return results, cgin
   

frac_co2 = inlet_conc_actual/10**6
results, cgin = modelrun(Q_g = 10.84, Q_l = 505.4,
                         pH = 13.01, times = outlet_t_sec,frac_co2 = frac_co2,constant_res_pH=True,
                         enh_method='PFO', wet_eff = 1, Kga = 'onda',recirc=True
                         )


pH = results['pH_profile'] # outlet is at 0 and inlet is at -1
x = results['cell_pos']
t = results['time']
pH_plot = pH[::-1,:] # since it is counter current so we flip it, this way it makes more sense
# print(t)

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Garamond'],
    'font.size': 12,

    'mathtext.fontset': 'stix',

    'axes.linewidth': 1,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'axes.spines.top': False,
    'axes.spines.right': False,

    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.direction': 'in',
    'ytick.direction': 'in',

    'legend.fontsize': 11,
    'legend.frameon': False,

    'lines.linewidth': 1.2,

    # THIS is the key difference
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

indices = [0,1,2,-1]
for i in indices:
    plt.plot(x,pH_plot[:,i], label = f'{t[i]} s')
plt.title('pH vs position\n Ql = 505.5 mL/min, Qg = 10.8 L/min')
plt.grid(False)
plt.ylabel('pH value')
plt.xlabel('Position [m]')

plt.legend(loc = 'upper right', frameon = False, bbox_to_anchor=(1, 0.95))

plt.show()


