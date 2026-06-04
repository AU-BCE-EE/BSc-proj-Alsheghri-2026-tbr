import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

temp = 25
TK = temp+273.15 
k3 = 4.315 * 10**13 * np.exp(- 6666/TK) / 1000               # Second-order rate constant for CO2 + OH⁻ → HCO3⁻  [m3/mol-s]
k1 = 6.672 * 10**12 * np.exp(- 9.724 * 10**3 / TK)           # First-order rate constant for CO2 + H2O → H2CO3   [s^-1]



c_co2 = 1

pH = np.arange(7,10+0.1,0.1)
pOH = 14-pH
c_oh = 10**(-pOH) * 1000


c_co2 = np.ones_like(c_oh)
r1 = k1 * c_co2 
r2 = k3 * c_co2 * c_oh

ratio = r2/r1




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


cross_pH = pH[np.argmin(np.abs(ratio - 1))]

plt.figure(figsize=(8,5))
plt.title(r'Reaction path dominance vs pH')
plt.plot(pH, ratio, label = r'$r_2/r_1$', color ='red')
plt.xlabel('pH value')
plt.ylabel(r'$r_2/r_1$')
plt.axhline(1, linestyle ='--', label = r'$r_2 = r_1$')
plt.axvline(cross_pH, linestyle='--', color='gray', label=f'Crossover pH ≈ {cross_pH:.1f}')
plt.axvspan(pH[0], cross_pH, alpha=0.05, color='blue', label=r'$r_1$ dominant')
plt.axvspan(cross_pH, pH[-1], alpha=0.05, color='red', label=r'$r_2$ dominant')
plt.legend(frameon = False, bbox_to_anchor=(0.3, 0.9))
plt.tight_layout()
plt.grid(False)
plt.show()