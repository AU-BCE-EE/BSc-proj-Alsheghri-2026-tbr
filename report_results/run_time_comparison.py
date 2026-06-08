import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import mods.mod_co2_first_draft as old_mod
import importlib
import matplotlib as mpl
importlib.reload(md)
from time import perf_counter
import numpy as np


times = np.arange(0,6060+60,60)
def modelrun(Q_g = 10.8,
             Q_l = 505,
             D   = 0.19,
             pH  = 13.0,
             vres = 20,
             nc = 10,
             times = times,
             frac_co2 = 0.025,
             cf = 1.0,
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
    nc = nc

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

    # ex_oh = 10**(-(14-pH)) * 1000 # mol/m3

    # cross sectional area of the reactor
    A = (np.pi*D**2)/4   # m2
    # flow velocity using flow rate and area 
    v_g = (Q_g * 1/1000 * 1/60) / A  # L/min * 1m3/1000L * 1min/60s = m3/s
    v_l = (Q_l * 1/10**6 * 1/60) / A # mL/min * 1m3/10^6mL * 1min/60s = m3/s 

    cf = cf
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
        cf = cf,
        typ='PR',
        counter=counter,
        recirc=recirc,
        enh_method=enh_method,
        constant_res_pH = constant_res_pH
    )

    return results, cgin





def modelrun_old(Q_g = 10.8,
                 Q_l = 505,
                 D   = 0.19,
                 pH  = 13.0,
                 vres = 20,
                 cf = 1,
                 nc = 10,
                 times = times,
                 frac_co2 = 0.025,
                 Kga = 'onda',
                 counter = True,
                 recirc = True):
    L     = 0.3       # m
    por_g = 0.86
    por_l = 0.05
    ssa   = 260       # m2/m3
    ssa = ssa
    nc    = nc

    cg0  = 9.9
    cl0  = 0.0
    clin = 0.0
    ccr  = 0.0

    M_co2 = 44.009
    temp = 22.0
    TK = temp + 273.15
    k_molar = 4.315 * 10**13 * np.exp(-6666/TK) 
    k_mass = k_molar * (1/M_co2) * (1/1000) 
    k2 = 'default'
    cf = cf
    
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
        L, por_g, por_l,cf, v_g, v_l, nc,
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



cells = np.arange(10, 90+10, 10)


t_old = []
t_new = []
for n in cells:
    start_new = perf_counter()
    res,cgin = modelrun(nc=n)
    t_new.append((perf_counter() - start_new))

    start_old = perf_counter()
    res_old,cgin_old = modelrun_old(nc=n)
    t_old.append((perf_counter() - start_old))






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

plt.figure(figsize=(6,4))
plt.semilogy(cells, t_old, '-', color='grey', label='M2')
plt.semilogy(cells, t_new, '-', color='red', label='M1')
plt.xlabel('Number of cells')
plt.ylabel('Runtime [s] (log scale)')
plt.title('Runtime comparison')
plt.grid(False)
plt.legend(frameon = False)
plt.tight_layout()
plt.show()