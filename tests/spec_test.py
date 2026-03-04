"""

Description:
    Test of speciation model functions
"""

import mods.speciation_model as sm 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import importlib

importlib.reload(sm)

TOTC = np.arange(1e-5, 1e-4 + 1e-5, 1e-5) * 1000         # mol/L * 1000 L/m3 = mol/m3 
TK = 298                                       # K
ex_oh = 3e-2

# equilibrium constants from https://github.com/sashahafner/NH3-RTM kinSpec()
K_hco3 = 10**(-353.5305 - 0.06092*TK + 21834.37/TK + 126.8339 * np.log10(TK) - 1684915/TK**2)
K_co3 = 10 **(-461.4176 - 0.093448*TK + 26986.16/TK + 165.7595*np.log10(TK) - 2248629/TK**2)
KW = 10**(-4.2195 - 2915.16/TK)

res = [] 
for i in range(len(TOTC)):
    res.append (sm.spec2_matrix(TOTC[i], K_hco3, K_co3, KW, ex_oh))

# save the result in a dataframe
df = pd.DataFrame(res)
df['TOTC'] = TOTC
df['ex_oh'] = ex_oh

# mass balance
df['mass_balance_residual'] = df.TOTC - df.c_h2co3 - df.c_hco3 - df.c_co3

# charge balance
df['charge_balance_residual'] = df.c_h - df.c_hco3 - 2*df.c_co3 - df.c_oh

print(df) # to see the results in a table
# df.to_csv('spec_model_results.csv', index=False)

print()
if np.isclose(max(abs(df.mass_balance_residual)), 0, atol=1e-9):
    print(f"Mass balance is effectively zero")
else: 
    print(f'Mass balance is {max(abs(df.mass_balance_residual))} which differs by {max(abs(df.mass_balance_residual)) - 1e-9}')
print()


if np.isclose(max(abs(df.charge_balance_residual)), 0, atol=1e-9):
    print(f"Charge balance is effectively zero")
else: 
    print(f'Charge balance is {max(abs(df.charge_balance_residual))} which differs by {max(abs(df.charge_balance_residual)) - 1e-9}')

# Check the equlibrium constants values
# Should be, around 6.4, 16.5 and 14 respectively... 
print()
print(f"K_hco3 = {K_hco3:.2e} (pKa1 = {-np.log10(K_hco3):.2f})")
print(f"K_co3  = {K_co3:.2e} (pKa2 = {-np.log10(K_co3):.2f})")
print(f"KW     = {KW:.2e} (pKW  = {-np.log10(KW):.2f})")


plt.figure(figsize=(6, 4))

plt.plot(df.TOTC, df.pH , 'ro-')
plt.xlabel('TOTC [mol/m3]')
plt.ylabel('pH value')
plt.grid()
# plt.savefig('pH_vs_TOTC.png')
plt.show()

plt.figure(figsize=(12, 4))
plt.plot(df.TOTC, df.c_h2co3/df.TOTC * 100, 'o-', label = 'H2CO3')
plt.plot(df.TOTC, df.c_hco3/df.TOTC * 100, 'o-', label = 'HCO3')
plt.plot(df.TOTC, df.c_co3/df.TOTC * 100, 'o-', label = 'CO3')
plt.ylim(-10,110)
plt.xlabel('Total carbonate [mol/m3]')
plt.ylabel('Species distribution')
plt.legend()
plt.grid()
# plt.savefig('Species_distribution.png')
plt.show()

# ==================================== test the new vectorized function 
res2 = sm.spec2_vector(TOTC, K_hco3, K_co3, KW, ex_oh)
df2 = pd.DataFrame(res2) 
df2['TOTC'] = TOTC 
df2['ex_oh'] = ex_oh 
df2['mass_balance_residual'] = df2.TOTC - df2.c_h2co3 - df2.c_hco3 - df2.c_co3 
# charge balance 
df2['charge_balance_residual'] = df2.c_h - df2.c_hco3 - 2*df2.c_co3 - df2.c_oh 
print(df2) # to see the results in a table 
print() 
if np.isclose(max(abs(df2.mass_balance_residual)), 0, atol=1e-9): 
    print(f"Mass balance is effectively zero") 
else: 
    print(f'Mass balance is {max(abs(df2.mass_balance_residual))} which differs by {max(abs(df2.mass_balance_residual)) - 1e-9}') 

print() 
if np.isclose(max(abs(df2.charge_balance_residual)), 0, atol=1e-9): 
    print(f"Charge balance is effectively zero") 
else: 
    print(f'Charge balance is {max(abs(df2.charge_balance_residual))} which differs by {max(abs(df2.charge_balance_residual)) - 1e-9}') 


plt.figure(figsize=(6, 4)) 
plt.plot(df2.TOTC, df2.pH , 'ro-') 
plt.xlabel('TOTC [mol/m3]') 
plt.ylabel('pH value') 
plt.grid() 
# plt.savefig('pH_vs_TOTC.png') 
plt.show() 

plt.figure(figsize=(12, 4)) 
plt.plot(df2.TOTC, df2.c_h2co3/df2.TOTC * 100, 'o-', label = 'H2CO3') 
plt.plot(df2.TOTC, df2.c_hco3/df2.TOTC * 100, 'o-', label = 'HCO3') 
plt.plot(df2.TOTC, df2.c_co3/df2.TOTC * 100, 'o-', label = 'CO3') 
plt.ylim(-10,110) 
plt.xlabel('Total carbonate [mol/m3]') 
plt.ylabel('Species distribution') 
plt.legend() 
plt.grid() 
# plt.savefig('Species_distribution.png') 
plt.show()