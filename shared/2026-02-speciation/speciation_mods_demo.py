"""
File name: speciation_mods_demo.py
Author: S. Hafner

Description:
    Demo of the example inorganic carbon speciation functions in 
    the `CO2_spec_mods.py` module.
"""

import numpy as np
import CO2_spec_mods as csm
#from importlib import reload

#reload(csm)   

# Equilibrium with an atmosphere that contains CO2 at a fixed partial pressure
# CO2 partial pressure in atm
pCO2 = 400E-6
csm.CO2eqgas(pCO2)
# So the pH is:
res = csm.CO2eqgas(pCO2)
-np.log10(res['H+'])

# Equilibrium with a fixed total dissolved inorganic C concentration
# Concentration in mol/kg
tot = 1E-5
csm.CO2eqtot(tot)
