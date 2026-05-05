import numpy as np
import pandas as pd

# import the data
data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\16_04_2026_co2_measurements.csv",
    sep=';'
)


data['Timestamp'] = pd.to_datetime(data['Timestamp'])

# start and end for the inlet measurment
start = "2026-04-16 13:00:18"
end   = "2026-04-16 13:07:18"

# extract all the values in the interval
inlet = data[(data["Timestamp"] >= start) & (data["Timestamp"] <= end)]

# mean inlet conc.
inlet_conc_no_cal = inlet['SCD30_CO2'].mean() # this is the diluted value, the real value is calculated using the dilution factor
inlet_conc = (inlet_conc_no_cal - 33.475)/0.9576 # here i used the devlabs calibration curve

Q_gas_bund = 0.382  # L/min
Q_air_mix = 2.736     # L/min

# Before the dilution
inlet_conc_actual = float((Q_gas_bund+Q_air_mix) / Q_gas_bund * inlet_conc) # ppm



# oulet concentrations
start_out = "2026-04-16 13:15:18"
end_out = "2026-04-16 14:27:18"
outlet = data[(data["Timestamp"] >= start_out) & (data["Timestamp"] <= end_out)]

Q_gas_top = 0.388  # L/min
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




# pH data 
time = ["13:20", "13:25", "13:30", "13:35", "13:40", "13:45", "13:50", 
        "13:55", "14:00", "14:05", "14:10", "14:15", "14:20", "14:25"]
pH   = [12.75, 12.75, 12.76, 12.76, 12.76, 12.76, 12.76,
        12.76, 12.76, 12.76, 12.76, 12.76, 12.77,12.77]
pH_data = pd.DataFrame({
    "time": time,
    "pH"  : pH
})


pH_data["Timestamp"] = pd.to_datetime(
    "2026-04-16 " + pH_data["time"]
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
             Kga = 'onda',
             wet_eff = 1,
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

    cg0 = 23

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

    Kga = Kga
    # 0.023

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
        wet_eff=wet_eff,
        typ='PR',
        counter=counter,
        recirc=recirc,
        enh_method=enh_method,
        constant_res_pH = constant_res_pH
    )

    return results, cgin

frac_co2 = inlet_conc_actual/10**6
results, cgin = modelrun(Q_g = 10.84, Q_l = 220.645,
                         pH = 13.00, times = outlet_t_sec,frac_co2 = frac_co2,constant_res_pH=True,
                         enh_method='PFO',wet_eff=1,recirc = False,Kga='onda'
                         )

results_2,cgin_2 = modelrun(Q_g = 10.84, Q_l = 220.645,
                         pH = 13.00, times = outlet_t_sec,frac_co2 = frac_co2,constant_res_pH=True,
                         enh_method='PFO',wet_eff=0.65,recirc = False,Kga='onda'
                         )


gas_1 = results['gas_conc']
pH_1 = results['pH_profile']
x_1       = results['cell_pos'] 
t_1       = results['time']

pH_outlet_1  = pH_1[0,:]
gas_outlet_1 = gas_1[-1,:]

removal_eff_vs_t_1 = 100 * (cgin - gas_outlet_1) / cgin

gas_2 = results_2['gas_conc']
pH_2 = results_2['pH_profile']
x_2       = results_2['cell_pos'] 
t_2       = results_2['time']

pH_outlet_2  = pH_2[0,:]
gas_outlet_2 = gas_2[-1,:]

removal_eff_vs_t_2 = 100 * (cgin_2 - gas_outlet_2) / cgin_2




# ======================= old model ==============================
import mods.mod_co2_first_draft as old_mod
def modelrun_old(Q_g = 10,
                 Q_l = 350,
                 D   = 0.19,
                 pH  = 13.7,
                 vres = 20,
                 wet_eff = 1,
                 times = outlet_t_sec,
                 frac_co2 = 0.05,
                 Kga = 'onda',
                 counter = True,
                 recirc = True):
    L     = 0.3       # m
    por_g = 0.86
    por_l = 0.05
    ssa   = 260       # m2/m3
    ssa = ssa*wet_eff
    nc    = 60

    cg0  = 23
    cl0  = 0.0
    clin = 0.0
    ccr  = 0.0

    M_co2 = 44.009
    temp = 22.0
    TK = temp + 273.15
    k_molar = 4.315 * 10**13 * np.exp(-6666/TK) 
    k_mass = k_molar * (1/M_co2) * (1/1000) 
    k2 = 'default'
    
    henry = [3.3e-2, 2400]


    dens_l = 997.0     # kg/m3
    pres   = 1.0       # bar
    pKa    = 6.35      # pKa for CO2/HCO3- (needed by old model function signature)

    R = 8.314e-5       # m3 * bar / K-mol

    # ================= Derived values =================
    # cross sectional area of the reactor
    A = (np.pi * D**2) / 4    # m2
    # flow velocity using flow rate and area 
    v_g = (Q_g * 1/1000 * 1/60) / A   # L/min * 1m3/1000L * 1min/60s = m3/s
    v_l = (Q_l * 1/10**6 * 1/60) / A  # mL/min * 1m3/10^6mL * 1min/60s = m3/s 

    v_res = vres / 1000 / A   # Convert reservoir volume to m3 per m2 cross section
    
    # Inlet gas concentration
    cgin = pres / ((temp + 273.15) * R) * frac_co2 * M_co2  # g/m3
    
    results = old_mod.tfmod(
        L, por_g, por_l, v_g, v_l, nc,
        cg0, cl0,
        cgin, clin,
        k_mass, Kga, henry, pKa, pH,
        temp, dens_l,
        times,
        kg='onda',
        kl='onda',
        ae='onda',
        v_res=v_res,
        k2=k2,
        ccr=ccr,
        pres=pres,
        ssa=ssa,
        typ='PR',
        counter=counter,
        recirc=recirc
    )

    return results, cgin



results_old, cgin_old = modelrun_old(Q_g = 10.84, Q_l = 220.645,
                         pH = 13, times = outlet_t_sec,frac_co2 = frac_co2,
                        Kga = 'onda',wet_eff = 0.65, recirc=False, counter = True
                         )

gas_old = results_old['gas_conc']
x_old       = results_old['cell_pos'] 
t_old       = results_old['time']


gas_outlet_old = gas_old[-1,:]

removal_eff_vs_t_old = 100 * (cgin_old - gas_outlet_old) / cgin_old
# ============================= results ====================

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
plt.suptitle('pH = 13.0, Ql = 220.6 mL/min, Qg = 10.8 L/min', fontsize = 14)
plt.subplot(1,3,1)
plt.title(r'(a) Gas CO$_2$ concentration at the outlet')
plt.plot(t_1, gas_outlet_1, 'r-', label = "Model")
plt.plot(t_2,gas_outlet_2, 'k-', label = "Model, We = 0.65")
plt.plot(t_old, gas_outlet_old, color = 'grey', linestyle = '--', label = "Simple Model, We = 0.65")
plt.plot(t_1, outlet_conc_gm3, 'bo', label = "Experimental", markersize = 3 )
plt.ylabel(r'CO$_2$ conc. [g/m$^3$]')
# plt.ylim(15,25)
plt.xlabel('Time [s]')
# plt.legend()
plt.grid(False)


plt.subplot(1,3,2)
plt.title('(b) pH at the outlet')
plt.plot(t_1, pH_outlet_1, 'r-', label = "Model")
plt.plot(t_2, pH_outlet_2, 'k-', label = "Model, We = 0.65")
plt.plot(pH_data["t_sec"], pH_data["pH"],'bo', label="Experimental", markersize = 3)
plt.ylabel('pH')
plt.xlabel('Time [s]')
# plt.legend()
plt.grid(False)

plt.subplot(1,3,3)
plt.title(r'(c) CO$_2$ removal efficiency')
plt.plot(t_1, removal_eff_vs_t_1, 'r-', label = 'Model')
plt.plot(t_2, removal_eff_vs_t_2, 'k-', label = 'Model, We = 0.65')
plt.plot(t_old, removal_eff_vs_t_old, color = 'grey',linestyle = '--', label = "Simple Model, We = 0.65" )
plt.plot(t_1,removal_efficiency_experimental,'bo', label = 'Experimental', markersize = 2)
plt.ylabel('Removal efficiency [%]')
plt.xlabel('Time [s]')
plt.ylim(0, 100)
plt.grid(False)
plt.legend(loc = 'upper right', frameon = False)

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.show()
