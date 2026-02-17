
import numpy as np
import math
import pandas as pd
import sys
from scipy.integrate import solve_ivp 
           
# molar masses for all species 
M_co2   = 44.01                 # g/mol
M_h2co3 = 62.03                 # g/mol
M_co3   = 60.01                 # g/mol
M_hco3  = 61.02                 # g/mol
M_oh    = 17.01                 # g/mol

# unit conversions to fix the user input and the output.
def gpm3_to_molpm3(c, M):
    return c / M

def molpm3_to_gpm3(c, M):
    return c * M

# Rates function
# For CO2 absoption in aqueous solution of NaOH.
def rates(t, n, 
          v_g, v_l, 
          cgin, clin_co2,clin_TOTC, 
          vol_gas, vol_liq, vol_tot,
          k1, k3, K1, K4, Kga, v_res, 
          temp, henry, pres, 
          ssa, dens_l, por_l, por_g, counter = True, recirc = False):  # Change the inputs(Rates and eq. cosntants)
  
    # interpolation function
  def interpolation(t, c):                                             
    if type(c) is pd.core.frame.DataFrame:
        c = np.interp(t, c.iloc[:, 0], c.iloc[:, 1])
    elif type(c) is np.ndarray:
        c = np.interp(t, c[0], c[1])
    elif not (type(c) is int or type(c) is float):
        sys.exit(f'Error: input must be float, integer, numpy array, or a pandas data frame, but is {type(c)}')
    return c
  # If time-variable concentrations coming in are given, get interpolated values
  # (then we interpolate to get current values at time t)
  clin_co2  = interpolation(t, clin_co2)
  clin_TOTC = interpolation(t, clin_TOTC)
  cgin      = interpolation(t, cgin)

  # reaction rates for the reverse reactions given by the equlibrium constant.
  k2 = k1/K1
  k4 = k3/K4 

  # Number of cells (layers) (note integer division)
  nc = (len(n)-2)//3                                                   # double check this
  # Separate gas, liquid and reservoir state variables (mol) (mol/m2)
  ncg     = n[0:nc]                              # CO2 in the gas phase
  n_co2   = n[nc:(2 * nc)]                       # CO2 in the liquid phase
  n_TOTC  = n[(2*nc):(3*nc)]                     # TOTC in the liquid phase
  ncr     = n[(3*nc) : (3*nc)+2]                 # reservoir

  # the reservoir now should handle 2 species 
  # CO2 
  ncr_co2  = ncr[0]
  # TOTC
  ncr_TOTC = ncr[1]

  # Concentrations (mol/m3)
  ccg     = ncg   / vol_gas                  # CO2 in the gas phase mol/m3                                                                 ## CHANGE
  c_co2   = n_co2 / vol_liq                  # CO2 in the liquid phase mol/m3   
  c_TOTC  = n_TOTC/ vol_liq                  # TOTC in the liquid phase mol/m3

  # note that the concentrations of H2CO3, HCO3, CO3 and OH would be calculated
  # from the speciation model, and therefore for now i would e.g. write c_h2co3 which does 
  # not excist yet, but when i have made the speciation model this would be more clear. 
  # Then it would be somethign like c_h2co3 = spec['H2CO3']
     

   # -----start----------- THis part needs to be fixed after talk with feilberg --------------
   # Caclulate Kga if requested (Mass transfer coefficient enhancement)
   # For ordinary Onda
  if type(Kga) is str and Kga.lower() == 'onda':
      # Hard-wired constants
      g      = 9.81        # m / sec^2
      Dg     = 1.16E-5     # gas diffusion coefficient in m2 / sec; compound specific     
      Dliq   = 1.9E-9      # liquid diffusion coefficient                                  
      sigm_c = 0.75        # critical surface tension
      sigm_l = 0.0073      # surface tension
      R      = 0.083144    # Gas constant (L bar / K-mol)
      mw_g   = 28.97       # Air (gas mix) molecular weight (molar mass) (g/mol)

      dp     = 6 * (1 - por_g) / ssa  # characteristic packing length
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

     
      kh  = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15))  # mol/kg-bar as liq:gas   # vant hoff equation for calculating new henry constant
      kh  = kh * dens_l / 1000                                 # mol/L-bar
      Kaw = 1 / (kh * R * TK)                                  # Neutral air-water distribution
    #   alpha0 = 1 / (1 + 10**(pH - pKa))                                                                   
    #   Daw = alpha0 * Kaw                                                                                  
    #   For CO2 no acid base eq. needed at high pH all the CO2 reacts with OH-
      Daw = Kaw                      # air water distribution (simple at high pH)                           
    
      #Enhancement factor for chemical reaction
      #Allows chemical reaction to happen in the liquid phase. 
      Ha = (Dliq* k1 * c_oh)**0.5/kl   # Hatta number, assumes fast reaction (Kasper and Feilberg, 2019)
      E  = Ha/np.tanh(Ha)              # Enhancment factor, assumes fast reaction (Kasper and Feilberg, 2019)
      # breakpoint()
      Rtot = 1 / (kg * ae) + Daw / (kl * ae * E) + Daw / (k1 * por_l) # reference p181 in Seader et al book (resistance in serie two film)
      Kga  = 1 / Rtot # overall mass transfer coefficient

##### Check units of Kga !!!! ####### bacause the rate function is now in moles...
#------------------------------------ end ----------------------------------------------------
  dmcr = np.zeros(2)

  # Get reservoir liquid phase concentration to use at inlet if recirc = True.
  #If recirc = false, reservoir concentration is still calculated but not used
  if recirc:
      if counter:
           oi = 0  # index depends in flow direction
      else:
           oi = nc - 1
      if v_res > 0:
        cr_co2  = ncr_co2  / v_res                      # CO2 in the reservoir phase mol/m3   
        cr_TOTC = ncr_TOTC / v_res                      # TOTC in the reservoir phase mol/m3

        r1_res = k1 * cr_co2 - k2 * c_h2co3
        r2_res = k3 * cr_co2 * cr_oh - k4 * cr_hco3    # reaction 2 in the reservoir
        # co2 in the reservoir
        clin_co2 = ncr_co2 / v_res                  
        rxnr_co2 = (r1_res + r2_res)

        # TOTC in the reservoir
        clin_TOTC = ncr_TOTC/v_res
        rxnr_TOTC = (r1_res + r2_res)

        # reservoir derivatives
        # CO2
        dmcr_co2  = (c_co2[oi]-cr_co2) * abs(v_l) - rxnr_co2 * v_res
        # TOTC
        dmcr_TOTC = (c_TOTC[oi]-cr_TOTC)*abs(v_l) + rxnr_TOTC * v_res

        # now we combine them together in one array
        dmcr = np.array([dmcr_co2,dmcr_TOTC])
      elif v_res == 0:
          clin_co2  = c_co2[oi]
          clin_TOTC = c_TOTC[oi]
      else:
          sys.exit('v_res is negative, must be a positive float or 0')          

  # Derivatives
  # Set up empty arrays
  dmg = dm_co2 = dm_TOTC = g2l = np.zeros(nc)


  # Common term, mass transfer into liquid phase (mol/s)
  #mol/s  1/s    m3(t)     ----mol/m3(g)-----
  g2l = Kga * vol_tot * (ccg - c_co2 * Daw)   # (here we calculate the mass transfer rate from gas to liquid)
  # Gas phase derivatives (mol/s)
  # No reaction in gas phase
  # cddiff = concentration double difference (mol/m3)
  # cvec = array of cell concentrations with inlet air added
  # rxn = 0 for gas phase
  cvec_gas  = np.insert(ccg, 0, cgin)       # Add the inlet CO2 concentration
  advec_gas = - v_g * np.diff(cvec_gas)     # (advection)  m/s * mol/m3 = mol/s-m2 but since A = 1m2 
  dmg       = advec_gas - g2l               # Rate of change of CO2 mass in the gas phase       [mol/s]


  # CO2 in the liquid phase 
  # Liquid phase derivatives (mol/s)
  # Includes transport and reaction  
  v_l_eff = v_l 

  r1 = k1 * c_co2 - k2 * c_h2co3
  r2 = k3 * c_co2 * c_oh - k4 * c_hco3

  rxn_co2   = r1 + r2
  rxn_TOTC  = r1 + r2      

  if not counter:
     cvec_co2 = np.insert(c_co2, 0, clin_co2)
  else:
     v_l_eff  = - v_l
     cvec_co2 = np.insert(c_co2, nc, clin_co2)

  advec_co2 = - v_l_eff * np.diff(cvec_co2)             # advection for CO2 in the liquid
  # mol/s   mol/s      mol/s  mol/m3-s * m3
  dm_co2 = advec_co2 + g2l + rxn_co2 * vol_liq          # Rate of change of CO2 in the liquid 
  
  # TOTC in the liquid phase
  if not counter:
     cvec_TOTC = np.insert(c_TOTC, 0, clin_TOTC)
  else:
     v_l_eff   = - v_l
     cvec_TOTC = np.insert(c_TOTC, nc, clin_TOTC)

  advec_TOTC = - v_l_eff * np.diff(cvec_TOTC)           # advection for TOTC in the liquid
  dm_TOTC    = advec_TOTC  + rxn_TOTC * vol_liq         # Rate of change of TOTC in the liquid 

 # Combine gas and liquid and reservoir
  dm = np.concatenate([dmg, dm_co2,dm_TOTC,dmcr])               
  # dm = np.append (dm, dmcr)

  #if t / 3600 > 0.05:
  #    breakpoint()

  return dm

# Model function
# TBR
def tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl_co20, cl_TOTC0, cgin, 
          clin_co2,clin_TOTC, cr_co20,cr_TOTC0, k1, K1, k3, Kga, henry, pKa,temp, dens_l, 
          times, kg='onda', kl='onda', ae='onda', v_res = 0, 
          pres = 1., ssa = 1100, typ = 'TBD', counter = True, recirc = False):

   ## Note that units are defined per 1 m2 filter cross-sectional (total) area 
   ## Below, where 2 sets of units are given this applies to the first case
   ## For the second one, the cross-sectional area is used to normalize the unit
   ## The two are mathematically equivalent
   # L         = total longitudinal length/height of reactor/filter (m)
   # por_g     = gas phase porosity (m3/m3 = m3(g)/m3(t) where g = gas and t = total)
   # por_l     = liquid phase content (m3/m3 = m3(l)/m3(t) where l = liquid)
   # v_g       = superficial velocity in m/s
   # v_l       = superficial velocity in m/s
   # nc        = number of cells (layers)
   # -------------------------initial concentrations --------------------------
   # cg0       = initial compound concentration in gas phase (g/m3 = g(compound)/m3(g))
   # cl_co20   = initial CO2 concentration in liquid phase (g/m3)
   # cl_TOTC0  = initial TOTC concetraion in the liquid phase(mol/m3) 
   # -------------------------- inflow concentrations --------------------
   # cgin      = CO2 concentration in gas inflow (g/m3(g))
   # clin_co2  = CO2 concentration in liquid inflow (g/m3(l)) (ignored if recirc = True)
   # clin_TOTC = TOTC concentration in liquid inflow (mol/m3(l)) (ignored if recirc = True)
   # --------------------- Reservoir Initial Concentrations ----------------
   # cr_co20   = initial CO2 concentration in reservoir (g/m3)
   # cr_TOTC0  = initial TOTC concentration in reservoir (mol/m3)
   # -------------------- physical parameters ----------------------
   # Kga       = mass transfer coefficient for gas to liquid in gas phase units (1/s = g/s-m3(t) / g/m3(g))
   # k1        = first-order rate constant for CO2 + H2O -> H2CO3 (s^-1)
   # k3        = second-order rate constant for CO2 + OH^- -> HCO3^- (m3/mol-s)
   # K1        = equilibrium constant for CO2 + H2O -> H2CO3 (dimensionless)
   # henry     = Henry's law constant coefficients as [k_H at 25 C, d(ln(kH)) / d(1/T)] as in NIST Chemistry Web Book
   # temp      = temperature (degrees C)
   # dens_l    = solution (liquid) density (kg/m3)
   # pres      = total pressure (bar?)
   # ssa       = particle specific surface area (m2 surface / m3 bulk volume)
   # count     = Boolean for countercurrent flow
   # v_res     = volume of a reservoir for the liquid phase (m3 pr m2 cross sectional area), ignored if recirc = False
   # typ       = type of filtermaterial if Kim and Deshusses is used for mass transfer estimation (options:LR, PUF, PR, PCB or PCR for lava rock, polyurethane foam, pall ring, porous ceramic bead and porous ceramic ring)

   # Save input arguments for echoing in output
   args_in = locals()

   # Constants
   # Ideal gas constant (L bar / K-mol)
   R = 0.083144 
   
   
   # Retention time (s)
   rt_gas = L * por_g / v_g
   rt_liq = L * por_l / v_l
   
   # Make sure some inputs are numeric to avoid integer math bug
   cg0      = float(cg0)
   cl_co20  = float(cl_co20)
   cl_TOTC0 = float(cl_TOTC0)
   cr_co20  = float(cr_co20)
   cr_TOTC0 = float(cr_TOTC0)
   k1       = float(k1)
   k3       = float(k3)
   K1       = float(K1)
   pKa      = float(pKa)
   temp     = float(temp)
   dens_l   = float(dens_l)
   times    = np.array(times).astype(float)
   henry    = np.array(henry).astype(float)

   # Temperature and Henry's law constant
   TK   = temp + 273.15
   kh   = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15))     # mol/kg-bar as liq:gas  (vant hoff)
   kh   = kh * dens_l / 1000                                    # mol/L-bar
   Kaw  = 1 / (kh * R * TK)                                     # dimensionless, gas:liq, e.g., g/L / g/L or g/m3 per g/m3
   
   # Temperature dependent constants
   Kw = 10**(-4.2195 - 2915.16/TK) # KW = 10^-pKw
   rho_water = 1000 * (1 - ((temp + 288.9414) / (508929.2 * (temp + 68.12963))) * (temp - 3.9863)**2)
   K2 = math.exp(-12092.1/TK -36.786 * math.log(TK) + 235.482) * rho_water  # equlibrium constant for CO2 + H2O -> HCO3^- + H^+
   K4 = K2/Kw                                                               # equilibrium constant for CO2 + OH^- -> HCO3^- (dimensionless)
   K_hco3 = 10**(-353.5305 - 0.06092*TK + 21834.37/TK + 126.8339 * np.log10(TK) - 1684915/TK**2)
   K_co3 = 10 **(-461.4176 - 0.093448*TK + 26986.16/TK + 165.7595*np.log10(TK) - 2248629/TK**2)

   
   # Create cells
   # we create nc equally spaced cells over the lenght L
   x    = np.linspace(0, L, nc + 1)  # nc + 1 values
   dx   = np.diff(x)[0]              # dx, single value, same for all cells  # cell lenght
   x    = x[1:(nc + 1)] - dx / 2     # Center position in m   # cell position 
   
   # Cell volumes (m3) (m3(g/l/t)/m2(t))
   vol_tot = dx * 1  # cell volume dx * A_cross sectional , the cross sectional area is normalized to 1 m^2
   # Gas
   vol_gas = vol_tot * por_g 
   # Liquid
   vol_liq = vol_tot * por_l 

   # convert the users input from g/m3 to mol/m3 so it fits the units of the model
   cg0      = gpm3_to_molpm3(cg0, M_co2)
   cgin     = gpm3_to_molpm3(cgin, M_co2)
   cl_co20  = gpm3_to_molpm3(cl_co20,  M_co2)
   clin_co2  = gpm3_to_molpm3(clin_co2,  M_co2)



   cr = np.array([
        gpm3_to_molpm3(cr_co20,  M_co2),
        cr_TOTC0,
    ]) * v_res  # mol

   # Uniform initial concentrations throughout reactor
   ccg    = np.full((nc), cg0)
   c_co2  = np.full((nc), cl_co20)   
   c_TOTC = np.full((nc), cl_TOTC0) 

   # convert concentration to moles
   ncg      = ccg * vol_gas            # mol/m3 * m3 = mol
   ncl_co2  = c_co2 * vol_liq
   ncl_TOTC = c_TOTC * vol_liq


   # Initial state variable array
   y0 = np.concatenate([ncg, ncl_co2,ncl_TOTC,cr])
   #  y0 = np.append (y0, ncr)


   
   # Solve/integrate
   out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                   t_eval = times, 
                   args = (v_g, v_l,cgin, clin_co2,clin_TOTC,
                           vol_gas, vol_liq, vol_tot, k1, k3, K1, K4, Kga, v_res,
                           temp, henry, pres,ssa, dens_l, por_l, por_g, counter, recirc
                           ),
                   method = 'Radau')
   
   # Extract moles of compound [position, time]
   ncgt      = out.y[0:nc]       
   ncl_co2t  = out.y[nc:(2 * nc)]
   ncl_TOTCt = out.y[(2*nc):(3*nc)]
   ncr       = out.y[(3*nc):(3*nc+2)]
   nctot     = np.sum(ncgt,0) + np.sum(ncl_co2t,0) + np.sum(ncl_TOTCt,0)  #total mass of compuond in the entire column as a function of time. why? 

   # Get concentrations vs. time
   # Gas
   ccgt      = ncgt / np.transpose(np.tile(vol_gas, (ncgt.shape[1], 1)))         # why not just ncgt/vol_gas
   # co2
   ccl_co2t  = ncl_co2t / np.transpose(np.tile(vol_liq, (ncl_co2t.shape[1], 1))) # why not just ncl_co2t/vol_liq
   # TOTC
   ccl_TOTCt = ncl_TOTCt / np.transpose(np.tile(vol_liq, (ncl_TOTCt.shape[1], 1)))

   # Total
  #  cctt = nctot / np.transpose(np.tile(vol_tot, (nctot.shape[1], 1)))

   nct = np.concatenate([ncgt, ncl_co2t,ncl_TOTCt])

   
   # Return results as a dictionary
   return {'gas_conc'    : molpm3_to_gpm3(ccgt,M_co2),
          'co2_liq_conc' : molpm3_to_gpm3(ccl_co2t,M_co2),
          'TOTC_liq_conc': ccl_TOTCt,
          'cell_pos'     : x,
          'time'         : times,
          'inputs'       : args_in, 
          'pars'         : {'gas_rt': rt_gas, 'liq_rt': rt_liq, 'Kga': Kga, 'Kaw': Kaw, 'ae':ae,'kg':kg,'kl':kl}}