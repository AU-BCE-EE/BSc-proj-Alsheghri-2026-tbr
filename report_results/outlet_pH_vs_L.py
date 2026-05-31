import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)
import numpy as np

def modelrun(Q_g = 10,
             Q_l = 350,
             D   = 0.19,
             pH  = 13.7,
             vres = 20,
             times = np.arange(0,6060+60,60),
             frac_co2 = 0.05,
             L = 0.30,
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
    por_g = 0.86
    por_l = 0.05
    ssa   = 260     # m2/m3
    nc    = 10

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

Length = np.arange(0.10,1.00+0.05,0.05)
Ql = [220.6, 505.5]

pH_diff_pH13_Q505 = []
pH_diff_pH13_Q220 = []
for L in Length:
    for Q in Ql:
        results, cgin = modelrun(Q_g = 10.8, Q_l = Q , pH = 13,
                                frac_co2=0.025, L = L)
        pH_inlet = results['pH_profile'][-1,-1]
        pH_outlet = results['pH_profile'][0,-1]
        pH_diff = pH_inlet - pH_outlet
        if Q == 505.5:
            pH_diff_pH13_Q505.append(pH_diff)
        else: 
            pH_diff_pH13_Q220.append(pH_diff)


pH_diff_pH12_5_Q505 = []
pH_diff_pH12_5_Q220 = []
for L in Length:
    for Q in Ql:
        results, cgin = modelrun(Q_g = 10.8, Q_l = Q , pH = 12.5,
                                frac_co2=0.025, L = L)
        pH_inlet = results['pH_profile'][-1,-1]
        pH_outlet = results['pH_profile'][0,-1]
        pH_diff = pH_inlet - pH_outlet
        if Q == 505.5:
            pH_diff_pH12_5_Q505.append(pH_diff)
        else: 
            pH_diff_pH12_5_Q220.append(pH_diff)

from matplotlib.lines import Line2D
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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), sharey=True) 
fig.suptitle('Steady-state pH difference between reactor inlet and outlet')


ax1.set_title('(a) Inlet pH = 13')
ax1.plot(Length, pH_diff_pH13_Q505, color='blue',
          label =r'Q$_l$ = 505 mL/min')
ax1.plot(Length, pH_diff_pH13_Q220, color='red',
             label = r'Q$_l$ = 220 mL/min')
ax1.legend(loc = 'upper left', fontsize='large', frameon = False, bbox_to_anchor=(0, 1))
ax1.set_xlabel('Packed bed length [m]')
ax1.set_ylabel(r'$pH_{in}$ - $pH_{out}$ ')
ax1.grid(False)


ax2.set_title('(b) Inlet pH = 12.5')
ax2.plot(Length, pH_diff_pH12_5_Q505, color='blue', 
            label = r'Q$_l$ = 505 mL/min ')
ax2.plot(Length, pH_diff_pH12_5_Q220, color='red',
          label = r'Q$_l$ = 220 mL/min')

ax2.set_xlabel('Packed bed length [m]')
ax2.grid(False)
ymin = 0.25
ymax = 2.25
yticks = np.arange(ymin, ymax+0.25, 0.25)
ax1.set_yticks(yticks)
# plt.legend(loc = 'upper left', fontsize='small', frameon = False, bbox_to_anchor=(0, 1))
plt.tight_layout()
plt.show()