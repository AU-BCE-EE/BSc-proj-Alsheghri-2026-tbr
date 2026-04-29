
import numpy as np
import math
import pandas as pd
import sys
from scipy.integrate import solve_ivp 
           

# Rates function
# Arg order: time, state variable, then arguments
# For CO2 absoption in aqueous solution of NaOH at high pH (>10)
def rates(t, mc, v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, pH, Kga, v_res, k2, temp, henry, pres, pKa, ssa, dens_l, por_l, por_g, counter = True, recirc = False):

  # If time-variable concentrations coming in are given, get interpolated values (Here we handle the values that may vary with time)
  # (then we interpolate to get current values at time t)
  if type(clin) is pd.core.frame.DataFrame:
    clin = np.interp(t, clin.iloc[:, 0], clin.iloc[:, 1])
  elif type(clin) is np.ndarray:
    clin = np.interp(t, clin[0], clin[1])
  elif not (type(clin) is int or type(clin) is float):
    sys.exit('Error: clin input must be float, integer, numpy array, or a pandas data frame, but is none of these')

  if type(cgin) is pd.core.frame.DataFrame:
    cgin = np.interp(t, cgin.iloc[:, 0], cgin.iloc[:, 1])
  elif type(cgin) is np.ndarray:
    cgin = np.interp(t, cgin[0], cgin[1])
  elif not (type(cgin) is int or type(cgin) is float):
    sys.exit('Error: cgin input must be float, integer, numpy array, or a pandas data frame, but is none of these')
    
  if type(pH) is pd.core.frame.DataFrame:
       pH = np.interp(t, pH.iloc[:, 0], pH.iloc[:, 1])
  elif type(pH) is np.ndarray:
       pH = np.interp(t, pH[0], pH[1])
  elif not (type(pH) is int or type(pH) is float):
       sys.exit('Error: pH input must be float, integer, numpy array, or a pandas data frame, but is none of these')

# (reaction rate, so if k2 is not specified then k2 = k which means that the same constant for both ionized and neutral species)
# For CO2 absorption: k is second-order rate constant (m3/(g*s)) for CO2 + OH- reaction
# k2 is not used for CO2 but kept for compatibility
  if type(k2) is str and k2.lower == 'default':
        k2=k
  
  # Calculate OH conentration from pH for CO2 kinetics
  # ___________________________ NEW _______
  pOH = 14-pH
  M_OH = 17.007                         # g/mol  molar mass of OH
  c_OH = 10**(-pOH) * M_OH * 1000       # g/m3   Concentration of OH ions 
  k_eff = k*c_OH                        # 1/s    The pseudo reaction rate constnat  
  # _______________________________________               
  
   # Caclulate Kga if requested (Mass transfer coefficient enhancement)
   #For ordinary Onda
  if type(Kga) is str and Kga.lower() == 'onda':
      # Hard-wired constants
      g = 9.81        # m / sec^2
      Dg = 1.16E-5    # gas diffusion coefficient in m2 / sec; compound specific      # changed
      Dliq = 1.89E-9   # liquid diffusion coefficient                                  # changed
      sigm_c = 0.029   # critical surface tension
      sigm_l = 0.072 # surface tension
      R = 0.083144    # Gas constant (L bar / K-mol)
      mw_g = 28.97    # Air (gas mix) molecular weight (molar mass) (g/mol)

      dp = 6 * (1 - por_g) / ssa  # characteristic packing length
      # empirical correlation for particle diameter
      if dp < 15:
          dp_emp = 2.0
      else:
          dp_emp = 5.23

      # temp in kelvin
      TK = temp + 273.15

      dens_g = pres * mw_g / (R * TK) # g/L = kg/m3     # density of gas 
      visc_g = 9.1E-8 * TK - 1.16E-5                    # empirical relation for gas viscosity vs TK
      visc_l = -2.55E-5 * TK + 8.51E-3                  # liquid viscosity 
      # Dimensionless numbers 
      Re = dens_l * v_l / (ssa * visc_l)                # reynold
      Fr = v_l * v_l * ssa / g                          # Froude number
      We = v_l * v_l * dens_l / (sigm_l * ssa)          # Weber number
      
      # some effective area for mass transfer? 
      ae = ssa * (1.0-2.71828**(-1.45 * (sigm_c / sigm_l)**0.75 *
                              Re**0.1 * Fr**-0.05 * We**0.2))
     
      #gas phase resistance
      kg = dp_emp * (v_g * dens_g / (ssa * visc_g))**0.7 \
          * (visc_g / (dens_g * Dg))**(1 / 3) * (ssa * dp)**-2 * ssa * Dg
      
      #liquid phase resistance
      kl = 0.0051 * (v_l * dens_l / (ae * visc_l))**(2/3) * (visc_l / (dens_l * Dliq))**(-0.5) * (ssa * dp)**0.4 * (dens_l / (visc_l * g))**(-1/3)

     
      kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas   # vant hoff equation for calculating new henry constant
      kh = kh * dens_l / 1000                                # mol/L-bar
      Kaw = 1 / (kh * R * TK)                                # Neutral air-water distribution
    #   alpha0 = 1 / (1 + 10**(pH - pKa))                                                                   # changed
    #   Daw = alpha0 * Kaw                                                                                  # changed
    #   For CO2 no acid base eq. needed at high pH all the CO2 reacts with OH-
      Daw = Kaw                      # air water distribution (simple at high pH)                           # changed
    
      #Enhancement factor for chemical reaction
      #Allows chemical reaction to happen in the liquid phase. 
      Ha = (Dliq*k_eff)**0.5/kl   # Hatta number, assumes fast reaction (Kasper and Feilberg, 2019)
      E = Ha/np.tanh(Ha)          # Enhancment factor, assumes fast reaction (Kasper and Feilberg, 2019)
      # breakpoint()
      Rtot = 1 / (kg * ae) + Daw / (kl * ae * E) + Daw / (k_eff * por_l) # reference p181 in Seader et al book (resistance in serie two film)
      Kga = 1 / Rtot  # overall mass transfer coefficient



  
  # Number of cells (layers) (note integer division)
  nc = mc.shape[0] // 2

  # Separate gas, liquid and reservoir state variables (g) (g/m2)
  mcg = mc[0:nc]                      # CO2 in the gas phase
  mcl = mc[nc:(2 * nc)]               # CO2 in the liquid phase
  mcr = mc[ (2*nc) : (2*nc)+1 ]       # CO2 in the reservoir if circulation is on

  # Concentrations (g/m3)
  ccg = mcg / vol_gas
  ccl = mcl / vol_liq
  
  dmcr = 0.0

  # Get reservoir liquid phase concentration to use at inlet if recirc = True. 
  #If recirc = false, reservoir concentration is still calculated but not used
  if recirc:
      if counter:
           oi = 0  # index depends in flow direction
      else:
           oi = nc - 1
      if v_res > 0:
        clin = mcr / v_res
        # rxnr = k * mcr * alpha0 + k2 * mcr * (1-alpha0)                   # changed
        rxnr = k_eff * mcr                                                  # changed
        #1/s * g/m2 * [] = g/(s*m2)
        #reservoir derivative (g/(s*m2)) 
        dmcr = (ccl[oi] - (mcr / v_res)) * abs(v_l) - rxnr
        #(g/m3 - g/m3) * m/s - g/(s*m2) = g/(s*m2) 
      elif v_res == 0:
          clin = ccl[oi]
      else:
          sys.exit('v_res is negative, must be a positive float or 0')
          


  # Derivatives
  # Set up empty arrays
  dmg = dml = g2l = np.zeros(nc)


  # Common term, mass transfer into liquid phase (g/s)
  #g/s  1/s    m3(t)     ----g/m3(g)-----
  g2l = Kga * vol_tot * (ccg - ccl * Daw)   # (here we calculate the mass transfer rate from gas to liquid)

  # Gas phase derivatives (g/s)
  # No reaction in gas phase
  # cddiff = concentration double difference (g/m3)
  # cvec = array of cell concentrations with inlet air added
  # rxn = 0 for as phase
  cvec = np.insert(ccg, 0, cgin) # Add the inlet CO2 concentration 
  advec = - v_g * np.diff(cvec)  # (advection) 
  dmg = advec - g2l              # Rate of change of CO2 mass in the gas phase 

  # Liquid phase derivatives (g/s)
  # Includes transport and reaction

#   rxn = k * mcl * alpha0        # reaction for neutral species      # change 
#   rxn2 = k2 * mcl * (1-alpha0)  # reaction for ionized species      # change
  rxn = k_eff * mcl               # Reaction for CO2 with OH          # change

  if not counter:
     cvec = np.insert(ccl, 0, clin)
  else:
     v_l = - v_l
     cvec = np.insert(ccl, nc, clin)

  advec = - v_l * np.diff(cvec)      # advection for CO2 in the liquid
  # dml = advec + g2l - rxn - rxn2 
  dml = advec + g2l - rxn            # Rate of change of CO2 mass in the liquid  # changed
  
 

 # Combine gas and liquid and reservoir
  dm = np.concatenate([dmg, dml])
  dm = np.append (dm, dmcr)

  #if t / 3600 > 0.05:
  #    breakpoint()

  return dm


def tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl0, cgin, clin, k, Kga, henry, pKa, pH, temp, dens_l, times, kg='onda', kl='onda', ae='onda', v_res = 0, k2 = 'default', ccr = 0, pres = 1., ssa = 1100, typ = 'TBD', counter = True, recirc = False):

   ## Note that units are defined per 1 m2 filter cross-sectional (total) area 
   ## Below, where 2 sets of units are given this applies to the first case
   ## For the second one, the cross-sectional area is used to normalize the unit
   ## The two are mathematically equivalent
   # L      = total longitudinal length/height of reactor/filter (m)
   # por_g  = gas phase porosity (m3/m3 = m3(g)/m3(t) where g = gas and t = total)
   # por_l  = liquid phase content (m3/m3 = m3(l)/m3(t) where l = liquid)
   # v_g    = superficial velocity in m/s
   # v_l    = superficial velocity in m/s
   # nc     = number of cells (layers)
   # cg0    = initial compound concentration in gas phase (g/m3 = g(compound)/m3(g))
   # cl0    = initial compound concentration in liquid phase (g/m3)
   # cgin   = compound concentration in gas inflow (g/m3)
   # clin   = compound concentration in liquid inflow (g/m3) (ignored if recirc = True)
   # Kga    = mass transfer coefficient for gas to liquid in gas phase units (1/s = g/s-m3(t) / g/m3(g))
   # k      = second-order liquid phase reaction rate constant (m3/g-s) for CO2+OH-
   # k2     = Not used but kept for compatibility
   # henry  = Henry's law constant coefficients as [k_H at 25 C, d(ln(kH)) / d(1/T)] as in NIST Chemistry Web Book
   # temp   = temperature (degrees C)
   # dens_l = solution (liquid) density (kg/m3)
   # pres   = total pressure (bar?)
   # ssa    = particle specific surface area (m2 surface / m3 bulk volume)
   # count  = Boolean for countercurrent flow
   # v_res  = volume of a reservoir for the liquid phase (m3 pr m2 cross sectional area), ignored if recirc = False
   # ccr    = initial concentration in the reservoir (g/m3)
   # typ    = type of filtermaterial if Kim and Deshusses is used for mass transfer estimation (options:LR, PUF, PR, PCB or PCR for lava rock, polyurethane foam, pall ring, porous ceramic bead and porous ceramic ring)

   # Save input arguments for echoing in output
   args_in = locals()

   # Constants
   # Ideal gas constant (L bar / K-mol)
   R = 0.083144 
   
   
   # Retention time (s)
   rt_gas = L * por_g / v_g
   rt_liq = L * por_l / v_l
   
   # Make sure some inputs are numeric to avoid integer math bug
   cg0 = float(cg0)
   cl0 = float(cl0)
   k = float(k)
   pKa = float(pKa)
   temp = float(temp)
   dens_l = float(dens_l)
   times = np.array(times).astype(float)
   henry = np.array(henry).astype(float)

   # Temperature and Henry's law constant
   TK = temp + 273.15
   kh = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas  (vant hoff)
   kh = kh * dens_l / 1000                                # mol/L-bar
   Kaw = 1 / (kh * R * TK)                                # dimensionless, gas:liq, e.g., g/L / g/L or g/m3 per g/m3
   
   
   
   # Create cells
   # we create nc equally spaced cells over the lenght L
   x = np.linspace(0, L, nc + 1)  # nc + 1 values
   dx = np.diff(x)[0]             # dx, single value, same for all cells  # cell lenght
   x = x[1:(nc + 1)] - dx / 2     # Center position in m   # cell position 
   
   # Cell volumes (m3) (m3(g/l/t)/m2(t))
   vol_tot = dx * 1  # cell volume dx * A_cross sectional , the cross sectional area is normalized to 1 m^2
   # Gas
   vol_gas = vol_tot * por_g 
   # Liquid
   vol_liq = vol_tot * por_l 
   
   # Compound conc (g/m3) and mass (g) (g/m2)
   # cc = concentration, mc = mass [position]
   # g = gas, l = liquid
   # Uniform initial concentrations throughout reactor
   ccg = np.full((nc), cg0)
   ccl = np.full((nc), cl0)
   # convert concentration to mass
   mcg = ccg * vol_gas
   mcl = ccl * vol_liq
   mcr = ccr * v_res

  # Initial state variable array
   y0 = np.concatenate([mcg, mcl])
   y0 = np.append (y0, mcr)


   
   # Solve/integrate
   out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                   t_eval = times, 
                   args = (v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, pH, Kga, v_res, k2, temp, henry, pres, pKa, ssa, dens_l, por_l, por_g, counter, recirc),
                   method = 'Radau')
   
   # Extract mass of compound [position, time]
   mcgt = out.y[0:nc]
   mclt = out.y[nc:(2 * nc)]
   mcrt = out.y[(2 * nc):(2 * nc)+1]
   mctot = np.sum(mcgt,0) + np.sum(mclt,0) #total mass of compuond in the entire column as a function of time. 

   # Get concentrations vs. time
   # Gas
   ccgt = mcgt / np.transpose(np.tile(vol_gas, (mcgt.shape[1], 1)))
   # Liquid
   cclt = mclt / np.transpose(np.tile(vol_liq, (mclt.shape[1], 1)))
   # Total
   cctt = mclt / np.transpose(np.tile(vol_tot, (mclt.shape[1], 1)))

   mct = np.concatenate([mcgt, mclt])
   
   # Return results as a dictionary
   return {'gas_conc': ccgt, 'liq_conc': cclt, 'gas_mass': mcgt, 'liq_mass': mclt, 
           'cell_pos': x, 'time': times, 'tot_mass' : mctot,
           'inputs': args_in, 
           'pars': {'gas_rt': rt_gas, 'liq_rt': rt_liq, 'Kga': Kga, 'Kaw': Kaw, 'ae':ae,'kg':kg,'kl':kl}}

