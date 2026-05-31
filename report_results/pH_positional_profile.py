import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
import matplotlib as mpl
importlib.reload(md)
import numpy as np

t = np.arange(0,6060+30,30)

def modelrun(Q_g = 10,
             Q_l = 350,
             D   = 0.19,
             pH  = 13.7,
             vres = 20,
             times = t,
             frac_co2 = 0.05,
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
   

frac_co2 = 0.025
results, cgin = modelrun(Q_g = 10.84, Q_l = 505.4,
                         pH = 13.01, times = t,frac_co2 = frac_co2,constant_res_pH=True,
                         enh_method='PFO', cf = 1, Kga = 'onda',recirc=True
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

indices = [0,1,2,3,4,-1]
plt.figure(figsize=(7,4))
for i in indices:
    plt.plot(x,pH_plot[:,i], label = f'{t[i]} s')
plt.title('Ql = 505.5 mL/min, Qg = 10.8 L/min')
plt.grid(False)
plt.ylabel('pH value')
plt.xlabel('Position [m]')

plt.legend(loc = 'best', frameon = False, bbox_to_anchor=(0.3, 0.75))

plt.show()


