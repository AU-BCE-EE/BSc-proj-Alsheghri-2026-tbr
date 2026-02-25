import mods.enhancement_mod as em
import numpy as np

# This is only for testing the consistency of the equations
# Inputs like Anders's excel file to check if i get similar results

pH = 11
c_oh = 10**-(14-pH) * 1000        # OH^- concentration [mol/m3]
c_co2 = 0.034                     # CO2 concentration [mol/m3]
k = 4*10**3 / 1000                # Second-order rate constant [m3/(mol*s)]
K = 6.60E7 / 1000                 # Equilibrium constant for CO2 + OH^- ⇌ HCO3^-
kl = 3.48e-6                      # Liquid-side mass transfer coefficient [m/s]

# Test PFO
E_pfo = em.enh_fac(c_oh, c_co2, k, K, kl, method='PFO')
print("Enhancement factor (PFO):", E_pfo)

# Test RSO
E_rso = em.enh_fac(c_oh, c_co2, k, K, kl, method='RSO')
print("Enhancement factor (RSO):", E_rso)