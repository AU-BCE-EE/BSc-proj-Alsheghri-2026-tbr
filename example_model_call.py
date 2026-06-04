import numpy as np
import matplotlib.pyplot as plt
# import the main model
import mods.mod_co2_main as md

# Call the main model function 

# lets take some inputs 

L = 0.3            # the packed bed length                          units: m
por_g = 0.86       # gas phase porosity 
por_l = 0.05       # liquid phase porosity
v_g = 0.006        # gas superficial velocity                       units: m/s
v_l = 0.00013      # liquid superficial velocity                    units : m/s
nc = 10            # number of grid cells 
cg0 = 0            # inital gas concentration in the column         units: g/m3
cl_co20 = 0        # initial CO2 concentration in liquid phase      units: g/m3
cl_TOTC0 = 0       # initial TOTC concetration in the liquid phase  units: mol/m3 
cgin = 40          # inlet gas CO2 concentration g/m3
ex_oh = 80         # excess OH⁻/Na+ concentration                   units: mol/m3
# ex_oh is the variable that sets the initial pH of the liquid. 
clin_co2 = 0       # CO2 concentration in the liquid in flow        units: g/m3
clin_TOTC = 0      # TOTC concnetration in the liquid in flow       units: mol/m3
cr_co20 = 0        # initial CO2 concentration in reservoir         units: g/m3
cr_TOTC0 = 0       # initial CO2 concentration in reservoir         units: g/m3
Kga = 'onda'         # overall mass transfer coefficient
# note that by giving Kga as a numeric value, the model skips onda correlations 
# and enhancement factor calculation. 

# henry law cosntant is given in the form of [k_H at 25°C, d(ln(kH))/d(1/T)] 
henry = [3.3e-2, 2400]
temp = 22       # temperture in                                     units: Celcius
dens_l = 997    # Liquid density                                    units: kg/m3
times = np.arange(0, 4500+60, 60) # Time points for which output is desired units: s

# these are not changed, they are like the old model
kg ='onda' # gas mass transfer coefficient
kl ='onda' # liquid mass transfer coefficient
ae ='onda' # effective interfacial area estimation 

# for the reserovir volume you need to deviede by the cross-sectional are of the column
D = 0.19 # m 
A_cross = D**2/4 * np.pi  # cross sectional area of the column m2
v_res_m3 = 20 * 1/1000  # reservoir volume in L * 1 m3/1000 L = m3
v_res = v_res_m3 / A_cross #                                        units:m3/m2

pres = 1.  # pressure                                               untis: bar
ssa = 260 # actual surface area of the packing material             units: m2/m3

cf = 1.0 # correction factor for over and underestimation of onda, this was mainly 
# used for the report (Osman's bachelor), keep it to 1. 

# this parameter was kept as in the original model, ignore it...
typ = 'PR'

# Counter or co current operation. If True you have counter current, if False you have co-current
counter = True,

# If True, liquid is recirculated. If False no recirculation
recirc = True
# here you can choose the enhancement factor approximation you want to use
# Look at the enhancement_mod.py for more info. 
enh_method = 'PFO'

# This was also for some report purposes. But if you set this parameter to True
# this would set the alkalinity in the reservoir to 0, meaning that you would have
# constant pH in the reservoir. This should be equal to setting recirc = False
# if you set it to False then pH would decrease in the reservoir as CO2 is absorbed.  
constant_res_pH = False


results = md.tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl_co20, cl_TOTC0, cgin, ex_oh, 
          clin_co2, clin_TOTC, cr_co20, cr_TOTC0, Kga, henry, temp, dens_l, 
          times, kg, kl, ae, v_res, 
          pres, ssa, cf, typ, counter, recirc, enh_method,
          constant_res_pH)

# the output is a dictionary. You can extract the results:


#[position, time]
gas_co2_conc = results['gas_conc']
liquid_co2_conc = results['co2_liq_conc']
liquid_pH = results['pH_profile']


cell_positions = results['cell_pos']
time = results['time']

# now you can make some plots if you like. 

# Ignore that the plots are ugly with no units and no labels. 

plt.plot(time, gas_co2_conc[-1,:]) # here we plot the outlet gas concentration
plt.show()

plt.plot(time, liquid_pH[0,:]) # here we plot the liquid pH at the reactor outlet, 
# not that since we have counter current, the liquid outlet is at position 0
plt.ylim(12.0,13.0)
plt.show()

# we could also plot the liquid phase CO2 concentration 
plt.plot(time,liquid_co2_conc[0,:])
plt.show()

# We could plot the pH profile in the column
# For i feel it makes more sense to flip the pH results when plotting the positional
# profile, in this case 0 is the inlet and 0.3 m is the outlet. 
plt.plot(cell_positions,liquid_pH[::-1][:,-1])
plt.show()




