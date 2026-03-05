import numpy as np
import matplotlib.pyplot as plt
import mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)

# parameters: 
L     = 2                    # m (half-bed used for evaluation)
por_g = 0.86                 # estimated (0.91 total - 0.05 liquid) # i dont know what this is maybe from void fraciton
por_l = 0.05                 # typical trickle-bed liquid holdup    # i dont know what this is
ssa   = 260                  # m2/m3 (plastic Pall rings) (from the thesis)

D_real = 0.19                   # m (from thesis)
A_real = np.pi * D_real**2 / 4  # 0.0284 m2  (from thesis)

# ---- Experimental flow case ----
Q_g_Lmin = 40.0            # gas flow [L/min] (from theis)
Q_l_Lmin = 0.1325          # liquid flow [L/min] (from thesis)

# Convert to m3/s
Q_g = Q_g_Lmin / 1000 / 60
Q_l = Q_l_Lmin / 1000 / 60

# Convert to superficial velocities (1 m2 model expects this)
v_g = 0.5   # m/s  (from Feilberg)
v_l = 1e-4  # m/s

nc = 20

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
pH    = 11
ex_oh = 10**(-(14-pH))*1000 # mol/m3

henry = [3.3e-2, 2400]      # NIST-based


temp   = 25.0              # temp in C
dens_l = 997.0             # density  kg/m3 
pres   = 1.0               # bar


Kga = 'onda'
times = np.linspace(0, 600, 200) # 10 min 


M_co2 = 44.009  # g/mol 
R     = 8.314 * 10**-5   # m3 * bar / K-mol



# ==================== Case 1 Biogas 30% - 50% ==========================
frac_biogas = 0.45    # fraction of co2 in the gas 
cgin_biogas = pres / ((temp + 273.15) * R) * frac_biogas * M_co2 

# results_biogas = md.tfmod(
#     L, por_g, por_l, v_g, v_l, nc,
#     cg0, cl_co20, cl_TOTC0,
#     cgin_biogas, ex_oh,
#     clin_co2, clin_TOTC,
#     cr_co20, cr_TOTC0,
#     Kga, henry, temp, dens_l,
#     times,
#     kg='onda',
#     kl='onda',
#     ae='onda',
#     v_res=v_res,          # now > 0,
#     pres=pres,
#     ssa=ssa,
#     typ='PR',
#     counter=False,         # counter-current gas–liquid
#     recirc=False,          # enable recycle
#     enh_method='PFO'
# )


# ==================== Case 2 cement production 14% - 30% ==========================
frac_cement = 0.20    # fraction of co2 in the gas 
cgin_cement = pres / ((temp + 273.15) * R) * frac_cement * M_co2 

# results_cement = md.tfmod(
#     L, por_g, por_l, v_g, v_l, nc,
#     cg0, cl_co20, cl_TOTC0,
#     cgin_cement, ex_oh,
#     clin_co2, clin_TOTC,
#     cr_co20, cr_TOTC0,
#     Kga, henry, temp, dens_l,
#     times,
#     kg='onda',
#     kl='onda',
#     ae='onda',
#     v_res=v_res,          # now > 0,
#     pres=pres,
#     ssa=ssa,
#     typ='PR',
#     counter=False,         # counter-current gas–liquid
#     recirc=False,          # enable recycle
#     enh_method='PFO'
# )

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

    plt.figure(figsize=(15, 5))
    plt.subplot(1,3,1)
    plt.title('gas conc. vs position')
    plt.plot(x, gas[:,-1], 'b-', linewidth = 1.8) 
    plt.ylabel('gas conc g/m3')
    plt.xlabel('position')
    plt.grid()

    plt.subplot(1,3,2)
    plt.title('liq conc. vs position')
    plt.plot(x, liq[:,-1], 'g-', linewidth = 1.8)
    plt.ylabel('liq conc g/m3')
    plt.xlabel('position')
    plt.grid()

    plt.subplot(1,3,3)
    plt.title('pH vs position')
    plt.plot(x, pH[:,-1], 'r-', linewidth = 1.8)
    plt.ylabel('pH value')
    plt.xlabel('position')
    plt.grid()

    plt.tight_layout()
    plt.show()

# biogas
# result_processing(results_biogas, cgin_biogas, 'Biogas')

# # cement 
# result_processing(results_cement,cgin_cement, 'Cement')

# low conc.
result_processing(results_low, cgin_low, 'Low conc.')
