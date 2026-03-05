import numpy as np
import matplotlib.pyplot as plt
import mods.mod_co2_main as md
import importlib
importlib.reload(md)


L = 2                    # m (half-bed used for evaluation)
por_g = 0.86             # estimated (0.91 total - 0.05 liquid) # i dont know what this is
por_l = 0.05             # typical trickle-bed liquid holdup
ssa = 260                # m2/m3 (plastic Pall rings)

D_real = 0.19                   # m
A_real = np.pi * D_real**2 / 4  # 0.0284 m2

# ---- Experimental flow case ----
Q_g_Lmin = 40.0            # gas flow [L/min]
Q_l_Lmin = 0.1325          # liquid flow [L/min]

# Convert to m3/s
Q_g = Q_g_Lmin / 1000 / 60
Q_l = Q_l_Lmin / 1000 / 60

# Convert to superficial velocities (1 m2 model expects this)
v_g = 0.5  # m/s
v_l = Q_l / A_real  # m/s

nc = 30

# Low conc. case
# 5 vol% CO2 at 1 bar, 25C ≈ 90 g/m3 (using ideal gas law)
cgin = 90.0              # g/m3
cg0 = 0.0

# ---- Liquid initial conditions ----
cl_co20 = 0.0
cl_TOTC0 = 0.0
clin_co2 = 0.0
clin_TOTC = 0.0
cr_co20 = 0.0
cr_TOTC0 = 0.0

# i dont know what to do for the reservoir. 
v_res = 0.0


ex_oh = 100              # mol/m3 (0.1 M NaOH, pH ≈ 13)
henry = [3.3e-2, 2400]   # NIST-based


temp = 25.0              # temp in C
dens_l = 997.0           # density 
pres = 1.0


Kga = 'onda'

times = np.linspace(0, 600, 200) # 10 min 

results = md.tfmod(
    L, por_g, por_l, v_g, v_l, nc,
    cg0, cl_co20, cl_TOTC0,
    cgin, ex_oh,
    clin_co2, clin_TOTC,
    cr_co20, cr_TOTC0,
    Kga, henry,temp, dens_l,
    times,
    kg='onda',
    kl='onda',
    ae='onda',
    v_res=v_res,
    pres=pres,
    ssa=ssa,
    typ='PR',
    counter=True,
    recirc=False,
    enh_method='PFO'
)

gas = results['gas_conc']
pH = results['pH_profile']
TOTC = results['TOTC_liq_conc']
z = results['cell_pos']
t = results['time']

gas_inlet = cgin
gas_out_initial = gas[-1, 0]
gas_out_final = gas[-1, -1]

removal_eff = 100 * (gas_inlet - gas_out_final) / gas_inlet

# pH
pH_initial = pH[0, 0]
pH_final = pH[-1, -1]



print("\n===== MODEL OUTPUT =====\n")

print(f"Gas superficial velocity : {v_g:.2f} m/s")
print(f"Liquid superficial velocity: {v_l:.6f} m/s")

print(f"\nGas inlet CO2         : {gas_inlet:.1f} g/m3")
print(f"Gas outlet initial      : {gas_out_initial:.2f} g/m3")
print(f"Gas outlet final        : {gas_out_final:.2f} g/m3")
print(f"CO2 removal efficiency  : {removal_eff:.2f} %")

print("\n--- pH ---")
print(f"Initial pH      : {pH_initial:.3f}")
print(f"Final pH        : {pH_final:.3f}")

print("\n======================================\n")