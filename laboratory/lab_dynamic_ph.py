import numpy as np
import pandas as pd 


# =============================== inlet data  =======================
inlet_data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\23_04_2026_inlet_co2_measurements.csv",
    sep = ';')

inlet_data['Timestamp'] = pd.to_datetime(inlet_data['Timestamp']) # convert to pandas datetime object

start_inlet = '2026-04-23 09:03:57'
end_inlet   = '2026-04-23 10:12:57'

# extract the useful data from the data file
inlet = inlet_data[(inlet_data['Timestamp'] >= start_inlet) & (inlet_data['Timestamp'] <= end_inlet)]
# convert time to second
inlet = inlet.copy()

inlet['t_sec'] = (inlet['Timestamp'] - inlet['Timestamp'].iloc[0]).dt.total_seconds()

Q_gas_bund = 0.390  # L/min
Q_air_mix  = 2.737   # L/min

# calculate conc. before dilution
inlet_conc_no_cal = inlet['SCD30_CO2'] * (Q_gas_bund + Q_air_mix)/Q_gas_bund  # ppm 

# no use the calibration curve
inlet_conc = (inlet_conc_no_cal- 33.475)/0.9576                               # ppm

frac_inlet = inlet_conc/10**6                                                 # fraction of co2 in inlet 

temp   = 22.0           # Celcius
dens_l = 997.0          # kg/m3
pres   = 1.0            # bar
M_co2 = 44.009          # g/mol
R     = 8.314e-5        # m3 * bar / K-mol

cgin = pres / ((temp + 273.15) * R) * frac_inlet * M_co2 # g/m3

# now make the data frame 

cgin_df = pd.DataFrame({
    'time': inlet['t_sec'].values,
    'cgin': cgin.values
})


# ======================== Now the outlet data.  ===================

outlet_data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\23_04_2026_outlet_co2_measurements.csv",
    sep = ';')

outlet_data['Timestamp'] = pd.to_datetime(outlet_data['Timestamp'])


start_outlet = "2026-04-23 11:52:51"
end_outlet   = "2026-04-23 13:06:51"

outlet = outlet_data[(outlet_data['Timestamp'] >= start_outlet) & (outlet_data['Timestamp'] <= end_outlet)]
outlet = outlet.copy()
outlet['t_sec'] = (outlet['Timestamp'] - outlet['Timestamp'].iloc[0]).dt.total_seconds()

Q_gas_top = 0.394   # L/min
Q_air_mix = 2.737    # L/min


# calculate conc. before dilution
outlet_conc_no_cal = outlet['SCD30_CO2'] * (Q_gas_top + Q_air_mix)/Q_gas_top  # ppm 

# no use the calibration curve
outlet_conc = (outlet_conc_no_cal- 33.475)/0.9576                             # ppm

temp   = 22.0           # Celcius
pres   = 1.0            # bar
M_co2 = 44.009          # g/mol
R     = 8.314e-5        # m3 * bar / K-mol

outlet_conc_gm3 = pres/(R*(temp+273.15)) * outlet_conc/10**6 * M_co2

# ======================================= pH =============================00
# pH data from 23/4 inlet experiment 
# time = ["9:06", "9:09", "9:12", "9:15", "9:19", "9:22", "9:25", "9:28",
#         "9:31", "9:34", "9:37", "9:40", "9:43", "9:46", "9:49", "9:52",
#         "9:55", "9:58", "10:01", "10:04", "10:07", "10:10"]

# pH   = [12.64, 12.60, 12.57, 12.54, 12.42, 12.36, 12.35, 12.34,
#         12.22, 12.25, 12.10, 11.92, 11.99, 11.92, 11.96, 11.92,
#         11.76, 11.70, 11.62, 11.82, 11.77, 11.60 ]

# # construct a dataframe
# pH_data = pd.DataFrame({
#     'time': time,
#     'pH'  : pH
# })

# # add timestamp like in inlet and convert to seconds.
# pH_data["Timestamp"] = pd.to_datetime("2026-04-23 " + pH_data["time"])
# pH_data['t_sec'] = (pH_data['Timestamp'] - inlet['Timestamp'].iloc[0]).dt.total_seconds()



# pH data from outlet experiment
time = ["12:03", "12:06", "12:09", "12:12", "12:15", "12:18", "12:21", "12:24",
        "12:27", "12:30", "12:33", "12:36", "12:39", "12:42", "12:45", "12:48",
        "12:51", "12:54", "12:57", "13:00", "13:03", "13:06"]

pH   = [12.29, 12.27, 12.22, 12.18, 12.22, 12.17, 12.02, 12.05, 
        12.05, 12.03, 12.00, 11.95, 11.82, 11.73, 11.69, 11.86, 
        11.81, 11.79, 11.62, 11.33, 11.15, 10.95]  

pH_data = pd.DataFrame({
    'time': time,
    'pH'  : pH
})

# add timestamp like in inlet and convert to seconds.
pH_data["Timestamp"] = pd.to_datetime("2026-04-23 " + pH_data["time"])
pH_data['t_sec'] = (pH_data['Timestamp'] - outlet['Timestamp'].iloc[0]).dt.total_seconds()

# ============================================= #
import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)
t_model = (inlet['Timestamp'] - inlet['Timestamp'].iloc[0]).dt.total_seconds().to_numpy()
def modelrun(Q_g = 10,
             Q_l = 350,
             D   = 0.19,
             pH  = 13.7,
             vres = 20,
             times = t_model,
             cgin = cgin_df,
             cf = 1.0,
             Kga = 'onda',
             counter=True,
             recirc=True,
             enh_method='PFO',
             constant_res_pH = True):
    L     = 0.3       # m
    por_g = 0.91
    por_l = 0.05
    ssa   = 260     # m2/m3
    nc    = 60

    cg0 = 25
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

    cf = cf
    # Kga = 0.02286
    Kga = Kga 
    v_res = vres/1000/A
    cgin = cgin
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
        cf = cf,
        typ='PR',
        counter=counter,
        recirc=recirc,
        enh_method=enh_method,
        constant_res_pH = constant_res_pH
    )

    return results, cgin





results, cgin = modelrun(Q_g = 10.84,
             Q_l = 220.645,
             D   = 0.19,
             pH  = 12.88,
             vres = 5,
             times = t_model,
             cgin = cgin_df,  
             cf = 0.65,
             Kga = 'onda',
             counter=True,
             recirc=True,
             enh_method='PFO',
             constant_res_pH = False)


gas = results['gas_conc']  # final position at all times
pH  = results['pH_profile']   
pH_plot     = results['pH_profile'][::-1,:]
x = results['cell_pos']
t = results['time']

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


plt.figure(figsize=(12, 5))
plt.suptitle('Ql = 220.6 mL/min, Qg = 10.8 L/min, cf = 0.65 \n pH = 12.88 (dynamic no adjusting)')
plt.subplot(1,2,1)
plt.title(r'Gas CO$_2$ concentration at the outlet')
plt.plot(t, gas[-1,:], 'r-', label = "Model")
plt.plot(outlet['t_sec'], outlet_conc_gm3, 'bo', label = "Experimental", markersize = 2 )
plt.ylabel(r'CO$_2$ conc. [g/m$^3$]')
plt.xlabel('Time [s]')
plt.grid(False)

plt.subplot(1,2,2)
plt.title('pH at the outlet')
plt.plot(t, pH[0,:], 'r-', label = "Model")
plt.plot(pH_data["t_sec"], pH_data["pH"],'bo', label="Experimental", markersize = 2)
plt.ylabel('pH')
plt.xlabel('Time [s]')
plt.legend(loc = 'upper right', frameon = False)
plt.grid(False)

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.show()