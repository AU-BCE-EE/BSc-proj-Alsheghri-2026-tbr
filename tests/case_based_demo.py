"""
Old case based demo version, here some calculations are done manual
"""

import numpy as np
import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)

# parameters: 
L     = 2                    # m (half-bed used for evaluation)
por_g = 0.91                 # estimated (0.91 total - 0.05 liquid) # i dont know what this is maybe from void fraciton
por_l = 0.05                 # typical trickle-bed liquid holdup    # i dont know what this is
ssa   = 260                  # m2/m3 (plastic Pall rings) (from the thesis)

D_real = 0.19                   # m (from thesis)
A_real = np.pi * D_real**2 / 4  # 0.0284 m2  (from thesis)

# ---- Experimental flow case ----
# Q_g_Lmin = 40.0            # gas flow [L/min] (from theis)
# Q_l_Lmin = 0.1325          # liquid flow [L/min] (from thesis)

# # Convert to m3/s
# Q_g = Q_g_Lmin / 1000 / 60
# Q_l = Q_l_Lmin / 1000 / 60

# Convert to superficial velocities (1 m2 model expects this)
v_g = 0.03   # m/s  (from Feilberg) # before 0.5
v_l = 1e-3  # m/s                   # before 1e-4

nc = 10

# gas phase 
cg0 = 0.0
# ---- Liquid initial conditions ----
cl_co20   = 0.0
cl_TOTC0  = 0.0
clin_co2  = 0.0
clin_TOTC = 0.0
cr_co20   = 0.0
cr_TOTC0  = 0.0
# i dont know what to do for the reservoir. 
v_res = 1 / 1000 # m3

# pH
pH    = 13.7
ex_oh = 10**(-(14-pH))*1000 # mol/m3

henry = [3.3e-2, 2400]      # NIST-based


temp   = 25.0              # temp in C
dens_l = 997.0             # density  kg/m3 
pres   = 1.0               # bar


Kga = 'onda'
times = np.linspace(0, 1200, 5) # 


M_co2 = 44.009  # g/mol 
R     = 8.314 * 10**-5   # m3 * bar / K-mol



# ==================== Case 1 Biogas 30% - 50% ==========================
frac_biogas = 0.45    # fraction of co2 in the gas 
cgin_biogas = pres / ((temp + 273.15) * R) * frac_biogas * M_co2 

results_biogas = md.tfmod(
    L, por_g, por_l, v_g, v_l, nc,
    cg0, cl_co20, cl_TOTC0,
    cgin_biogas, ex_oh,
    clin_co2, clin_TOTC,
    cr_co20, cr_TOTC0,
    Kga, henry, temp, dens_l,
    times,
    kg='onda',
    kl='onda',
    ae='onda',
    v_res=v_res,          # now > 0,
    pres=pres,
    ssa=ssa,
    typ='PR',
    counter=False,         # counter-current gas–liquid
    recirc=False,          # enable recycle
    enh_method='PFO'
)



m_gin  = results_biogas['m_gin']
m_gout = results_biogas['m_gout']
m_lout = results_biogas['m_lout']
m_tout = results_biogas['m_tout']

print("=" * 10)
print(f'Gas in: {m_gin}\n Gas out: {m_gout}\n liquid out: {m_lout}\n total: {m_tout}')
print(f'Balance: {m_tout - m_gin}')

print(f'actual conc. {results_biogas['co2_liq_conc']}')
print(f'equilibrium con. {results_biogas['eq_conc']}')


TOTC_actual = results_biogas['TOTC_liq_conc']
TOTC_eq     = results_biogas['TOTC_eq']

print("\n ========================== \n")
print(f'TOTC_liquid acutal = {TOTC_actual}')
print(f'TOTC_equilirbium {TOTC_eq}')

print(f'pH profile {results_biogas['pH_profile']}')

print("\n ========================== \n")
print(f'c_oh from c_h {results_biogas['c_oh']}')


# ==================== Case 2 cement production 14% - 30% ==========================
frac_cement = 0.20    # fraction of co2 in the gas 
cgin_cement = pres / ((temp + 273.15) * R) * frac_cement * M_co2 

results_cement = md.tfmod(
    L, por_g, por_l, v_g, v_l, nc,
    cg0, cl_co20, cl_TOTC0,
    cgin_cement, ex_oh,
    clin_co2, clin_TOTC,
    cr_co20, cr_TOTC0,
    Kga, henry, temp, dens_l,
    times,
    kg='onda',
    kl='onda',
    ae='onda',
    v_res=v_res,          # now > 0,
    pres=pres,
    ssa=ssa,
    typ='PR',
    counter=False,         # counter-current gas–liquid
    recirc=False,          # enable recycle
    enh_method='PFO'
)

# ==================== Case 3 low concnetration 1% - 10% ==========================

frac_low = 0.05    # fraction of co2 in the gas 
cgin_low = pres / ((temp + 273.15) * R) * frac_low * M_co2 

results_low = md.tfmod(
    L, por_g, por_l, v_g, v_l, nc,
    cg0, cl_co20, cl_TOTC0,
    cgin_low, ex_oh,
    clin_co2, clin_TOTC,
    cr_co20, cr_TOTC0,
    Kga, henry, temp, dens_l,
    times,
    kg='onda',
    kl='onda',
    ae='onda',
    v_res=v_res,          # now > 0, activates reservoir
    pres=pres,
    ssa=ssa,
    typ='PR',
    counter=False,         # counter-current gas–liquid
    recirc=False,          # enable recycle
    enh_method='PFO'
)

# ======================= resutl processing ============================ 

def result_processing(res,cgin,label):
    gas = res[  'gas_conc'   ]
    pH  = res[  'pH_profile' ]
    x   = res[  'cell_pos'   ] 
    liq = res[ 'co2_liq_conc']
    t   = res[    'time'     ]
    
    gas_inlet       = float(cgin)
    gas_out_initial = gas[-1,0]
    gas_out_final   = gas[-1,-1]

    removal_eff = 100 * (gas_inlet - gas_out_final) / gas_inlet

    pH_initial = pH[-1, 0]
    pH_final   = pH[-1, -1]


    print(f'\n===== {label} ======\n')

    print(f"\nGas inlet CO2         : {gas_inlet:.1f} g/m3")
    print(f"Gas outlet initial      : {gas_out_initial:.2f} g/m3")
    print(f"Gas outlet final        : {gas_out_final:.2f} g/m3")
    print(f"CO2 removal efficiency  : {removal_eff:.2f} %")

    print("\n--- pH ---")
    print(f"Initial pH      : {pH_initial:.3f}")
    print(f"Final pH        : {pH_final:.3f}")
    print("\n======================================\n")


    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = ['Garamond']
    mpl.rcParams['axes.linewidth'] = 1
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['xtick.labelsize'] = 11
    mpl.rcParams['ytick.labelsize'] = 11
    mpl.rcParams['legend.fontsize'] = 11    

    plt.figure(figsize=(12, 5))
    plt.subplot(1,3,1)
    plt.title('gas conc. vs position')
    plt.plot(x, gas[:,-1], 'b-', linewidth = 1.8) 
    plt.ylabel('gas conc g/m3')
    plt.ticklabel_format(style='plain', axis='y', useOffset=False)
    plt.xlabel('position')
    plt.grid()

    plt.subplot(1,3,2)
    plt.title('liq conc. vs position')
    plt.plot(x, liq[:,-1], 'g-', linewidth = 1.8)
    plt.ticklabel_format(style='plain', axis='y', useOffset=False)
    plt.ylabel('liq conc g/m3')
    plt.xlabel('position')
    plt.grid()

    plt.subplot(1,3,3)
    plt.title('pH vs position')
    plt.plot(x, pH[:,-1], 'r-', linewidth = 1.8)
    plt.ticklabel_format(style='plain', axis='y', useOffset=False)
    plt.ylabel('pH value')
    plt.xlabel('position')
    plt.grid()

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(12, 5))
    plt.subplot(1,3,1)
    plt.title('Gas concentraion at the outlet vs time')
    plt.plot(t, gas[-1,:])
    plt.ylabel('conc[g/m3]')
    plt.xlabel('t[s]')
    
    plt.subplot(1,3,2)
    plt.title('liquid concentraion at the outlet vs time')
    plt.plot(t,liq[-1,:])
    plt.ylabel('conc[g/m3]')
    plt.xlabel('t[s]')

    plt.subplot(1,3,3)
    plt.title('pH at the outlet vs time')
    plt.plot(t, pH[-1,:])
    plt.ylabel('pH')
    plt.xlabel('time[s]')

    plt.tight_layout()
    plt.show()

# biogas
result_processing(results_biogas, cgin_biogas, 'Biogas')

# cement 
result_processing(results_cement,cgin_cement, 'Cement')

# low conc.
result_processing(results_low, cgin_low, 'Low conc.')
