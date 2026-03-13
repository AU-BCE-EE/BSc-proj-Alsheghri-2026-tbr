
import numpy as np
import pandas as pd
import sys
from scipy.integrate import solve_ivp 
from . import speciation_model as sm
from . import enhancement_mod as em
from importlib import reload
          
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
          cgin, clin_co2, clin_TOTC, 
          vol_gas, vol_liq, vol_tot,
          k1, k2, k3, K4, K_hco3, K_co3, KW, ex_oh, Kga, v_res, 
          temp, henry, pres, 
          ssa, dens_l, por_l, por_g, counter = True, recirc = False, enh_method = 'PFO'):
  
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
  k4 = k3 / K4 

  # Number of cells (layers) (note integer division)
  nc = (len(n)-2)//3                                                   # double check this
  # Separate gas, liquid and reservoir state variables (mol) (mol/m2)
  ncg     = n[0        : nc]                       # CO2 in the gas phase
  n_co2   = n[nc       : (2 * nc)]                 # CO2 in the liquid phase
  n_TOTC  = n[(2 * nc) : (3 * nc)]                 # TOTC in the liquid phase
  ncr     = n[(3 * nc) : (3 * nc) + 2]             # reservoir

  # the reservoir now should handle 2 species 
  # CO2 
  ncr_co2  = ncr[0]
  # TOTC
  ncr_TOTC = ncr[1]

  # Concentrations (mol/m3)
  ccg     = ncg    / vol_gas                  # CO2 in the gas phase mol/m3                                                                 ## CHANGE
  c_co2   = n_co2  / vol_liq                  # CO2 in the liquid phase mol/m3   
  c_TOTC  = n_TOTC / vol_liq                  # TOTC in the liquid phase mol/m3

  # integrating the speciation model: 
  c_h2co3 = np.zeros(nc)
  c_hco3  = np.zeros(nc)
  c_co3   = np.zeros(nc)
  c_oh    = np.zeros(nc)
  c_h     = np.zeros(nc)
  pH      = np.zeros(nc)

  for i in range(nc):
      res = sm.spec2_matrix(c_TOTC[i], K_hco3, K_co3, KW, ex_oh)
      c_h2co3[i] = res['c_h2co3']
      c_hco3[i]  = res['c_hco3']
      c_co3[i]   = res['c_co3']
      c_oh[i]    = res['c_oh']
      c_h[i]     = res['c_h']
      pH[i]      = res['pH']

#   breakpoint()
   # -----start----------- THis part needs to be fixed after talk with feilberg --------------
   # here i calculate Daw outside the Kga function, because before this 
   # Daw was only calculted of Kga == 'onda' but sometimes we want to test with Kga = float

  R        = 0.083144        # # Gas constant (L bar / K-mol)
  TK       = temp + 273.15   # Temp in Kelvin

  kh       = henry[0] * np.exp(henry[1] * (1/TK - 1/298.15)) # mol/kg-bar as liq:gas   # vant hoff equation for calculating new henry constant
  kh       = kh * dens_l / 1000                              # mol/L-bar

  Kaw      = 1 / (kh * R * TK)                                # Neutral air-water distribution 
  Daw      = Kaw                                              # air water distribution 

   # Caclulate Kga if requested (Mass transfer coefficient enhancement)
   # For ordinary Onda  
  if isinstance(Kga,str) and Kga.lower() == 'onda':
      # Hard-wired constants
      g      = 9.81                   # m / sec^2
      Dg     = 1.16E-5                # gas diffusion coefficient in m2 / sec; compound specific     
      Dliq   = 1.89E-9                # liquid diffusion coefficient                                  
      sigm_c = 0.75                   # critical surface tension
      sigm_l = 0.0073                 # surface tension
      mw_g   = 28.97                  # Air (gas mix) molecular weight (molar mass) (g/mol)
      dp     = 6 * (1 - por_g) / ssa  # characteristic packing length

      # empirical correlation for particle diameter
      if dp < 15:
          dp_emp = 2.0
      else:
          dp_emp = 5.23

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
    
      #Enhancement factor for chemical reaction
      # interfacial co2 concentration
      # note we cant just use c_co2 because c_co2 is the bulk concentration.
      c_co2i = ccg / Kaw            # Anders also did this, but he said maybe there is a 
                                    # better way to do it here in the model... 
                                    # maybe we have to use something like two film approach
                                    # c_co2i = c_co2 + N_A /kl 
      if enh_method is not None:
          E = em.enh_fac(
              c_oh   = c_oh,
              c_co2  = c_co2i,
              k      = k3,
              K      = K4,
              kl     = kl,
              method = enh_method
              )
      else:
          E = 1

      Rtot = 1 / (kg * ae) + Daw / (kl * ae * E) + Daw / (k3 * np.maximum(c_oh,1e-10) * por_l) # reference p181 in Seader et al book (resistance in serie two film)
      Kga  = 1 / Rtot # overall mass transfer coefficient
    #   breakpoint()

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

        # integrating the speciation model: 
        res_reservoir = sm.spec2_matrix(cr_TOTC, K_hco3, K_co3, KW, ex_oh)
        cr_h2co3 = res_reservoir['c_h2co3']
        cr_hco3  = res_reservoir['c_hco3']
        cr_oh    = res_reservoir['c_oh'] 

        r1_res = k1 * cr_co2 - k2 * cr_h2co3
        r2_res = k3 * cr_co2 * cr_oh - k4 * cr_hco3
        # co2 in the reservoir
        clin_co2 = ncr_co2 / v_res                  
        rxnr_co2 = (r1_res + r2_res)

        # TOTC in the reservoir
        clin_TOTC = ncr_TOTC/v_res
        rxnr_TOTC = (r1_res + r2_res)

        # reservoir derivatives
        # CO2
        dmcr_co2  = (c_co2[oi] - cr_co2) * abs(v_l) - rxnr_co2 * v_res
        # TOTC
        dmcr_TOTC = (c_TOTC[oi] - cr_TOTC)*abs(v_l) + rxnr_TOTC * v_res

        # now we combine them together in one array
        dmcr = np.array([dmcr_co2, dmcr_TOTC])
      elif v_res == 0:
          clin_co2  = c_co2[oi]
          clin_TOTC = c_TOTC[oi]
      else:
          sys.exit('v_res is negative, must be a positive float or 0')          

  # Derivatives
  # Set up empty arrays
  dmg     = np.zeros(nc)
  dm_co2  = np.zeros(nc)
  dm_TOTC = np.zeros(nc) 
  g2l     = np.zeros(nc)


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

  if not counter:
     cvec_co2 = np.insert(c_co2, 0, clin_co2)
  else:
     v_l_eff  = - v_l
     cvec_co2 = np.insert(c_co2, nc, clin_co2)

  advec_co2 = - v_l_eff * np.diff(cvec_co2)             # advection for CO2 in the liquid
  # mol/s   mol/s      mol/s  mol/m3-s * m3
  dm_co2 = advec_co2 + g2l - rxn_co2 * vol_liq          # Rate of change of CO2 in the liquid 
  
  # TOTC in the liquid phase
  if not counter:
     cvec_TOTC = np.insert(c_TOTC, 0, clin_TOTC)
  else:
     v_l_eff   = - v_l
     cvec_TOTC = np.insert(c_TOTC, nc, clin_TOTC)

  advec_TOTC = - v_l_eff * np.diff(cvec_TOTC)           # advection for TOTC in the liquid
  dm_TOTC    = advec_TOTC  + rxn_co2 * vol_liq         # Rate of change of TOTC in the liquid 

 # Combine gas and liquid and reservoir
  dm = np.concatenate([dmg, dm_co2, dm_TOTC, dmcr])               
  # dm = np.append (dm, dmcr)

  #if t / 3600 > 0.05:
  #    breakpoint()

  return dm

# Model function
# TBR
def tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl_co20, cl_TOTC0, cgin, ex_oh, 
          clin_co2, clin_TOTC, cr_co20, cr_TOTC0, Kga, henry, temp, dens_l, 
          times, kg='onda', kl='onda', ae='onda', v_res = 0, 
          pres = 1., ssa = 1100, typ = 'TBD', counter = True, recirc = False, enh_method = 'PFO'):
   """"
    Simulate CO2 absorption in a trickle bed reaction/film with chemical reaction


    Parameters
    ----------
    L : float
        Total longitudinal length/height of reactor/filter (m)
    por_g : float
        Gas phase porosity (m³ gas / m³ total volume)
    por_l : float
        Liquid phase content (m³ liquid / m³ total volume)
    v_g : float
        Superficial gas velocity (m/s)
    v_l : float
        Superficial liquid velocity (m/s)
    nc : int
        Number of cells (layers)
        
    Initial concentrations
    ----------------------
    cg0 : float
        Initial CO2 concentration in gas phase (g/m³ gas)
    cl_co20 : float
        Initial CO2 concentration in liquid phase (g/m³ liquid)
    cl_TOTC0 : float
        Initial total inorganic carbon (TOTC) concentration in liquid phase (mol/m³ liquid)
    
    Inflow concentrations
    ---------------------
    cgin : float
        CO2 concentration in gas inflow (g/m³ gas)
    ex_oh : float
        Excess OH⁻ concentration (mol/m³ liquid)
    clin_co2 : float
        CO2 concentration in liquid inflow (g/m³ liquid) (ignored if recirc = True)
    clin_TOTC : float
        TOTC concentration in liquid inflow (mol/m³ liquid) (ignored if recirc = True)
    
    Reservoir initial concentrations
    --------------------------------
    cr_co20 : float
        Initial CO2 concentration in reservoir (g/m³)
    cr_TOTC0 : float
        Initial TOTC concentration in reservoir (mol/m³)
    
    Mass transfer parameters
    ------------------------
    Kga : float
        Overall mass transfer coefficient for gas to liquid in gas phase units (s⁻¹)
    henry : array_like
        Henry's law constant coefficients [k_H at 25°C, d(ln(kH))/d(1/T)] as in NIST Chemistry WebBook
    
    Physical parameters
    -------------------
    temp : float
        Temperature (°C)
    dens_l : float
        Liquid density (kg/m³)
    pres : float, optional
        Total pressure (bar). Default is 1.0
    ssa : float, optional
        Particle specific surface area (m² surface / m³ bulk volume). Default is 1100
    
    Numerical parameters
    --------------------
    times : array_like
        Time points for which output is desired (s)
    
    Model options
    -------------
    kg : str, optional
        Method for gas phase mass transfer coefficient estimation. 
        Options: 'onda' (Onda's correlation). Default is 'onda'
    kl : str, optional
        Method for liquid phase mass transfer coefficient estimation.
        Options: 'onda' (Onda's correlation). Default is 'onda'
    ae : str, optional
        Method for effective interfacial area estimation.
        Options: 'onda' (Onda's correlation). Default is 'onda'
    v_res : float, optional
        Volume of reservoir for liquid phase (m³ per m² cross-sectional area). 
        Ignored if recirc = False. Default is 0
    typ : str, optional
        Type of filter material for Kim & Deshusses mass transfer estimation.
        Options: 'LR' (lava rock), 'PUF' (polyurethane foam), 'PR' (pall ring), 
        'PCB' (porous ceramic bead), 'PCR' (porous ceramic ring). Default is 'TBD'
    counter : bool, optional
        If True, simulate countercurrent flow. If False, simulate cocurrent flow.
        Default is True
    recirc : bool, optional
        If True, liquid is recirculated through a reservoir. If False, single-pass.
        Default is False
    enh_method : str, optional
        Method for calculating the enhancement factor due to chemical reaction.
        Options:
        - 'PFO' : Pseudo first-order reaction
        - 'RSO' : Reversible second-order reaction
        Default is 'PFO'
    
    Returns
    -------
    dict
        Dictionary containing simulation results with the following keys:
        
        gas_conc : ndarray
            Gas phase CO2 concentration over time and position (g/m³ gas)
            Shape: (nc, len(times))
        co2_liq_conc : ndarray
            Liquid phase CO2 concentration over time and position (g/m³ liquid)
            Shape: (nc, len(times))
        pH_profile : ndarray
            pH over time and position
            Shape: (nc, len(times))
        TOTC_liq_conc : ndarray
            Total inorganic carbon concentration in liquid phase (mol/m³ liquid)
            Shape: (nc, len(times))
        cell_pos : ndarray
            Center positions of each cell (m)
            Shape: (nc,)
        time : ndarray
            Time points corresponding to simulation output (s)
            Shape: (len(times),)
        inputs : dict
            Copy of all input parameters for reference
        pars : dict
            Calculated parameters including:
            - gas_rt : Gas phase retention time (s)
            - liq_rt : Liquid phase retention time (s)
            - Kga : Mass transfer coefficient (s⁻¹)
            - Kaw : Dimensionless air-water partition coefficient
            - ae : Effective interfacial area
            - kg : Gas phase mass transfer coefficient
            - kl : Liquid phase mass transfer coefficient
    
    Notes
    -----
    - All units are normalized per 1 m² filter cross-sectional area
    - The model accounts for temperature-dependent equilibrium constants
    - Chemical speciation (pH, carbonate species) is calculated using the 
      spec2_matrix function from an imported speciation model (sm)
    - The system of ODEs is solved using the 'Radau' method from scipy.integrate.solve_ivp
    
    References
    ----------
    - NIST Chemistry WebBook for Henry's law constants
    - Onda et al. (1968) for mass transfer correlations
    - Kim & Deshusses (2003) for filter material specific correlations
   """

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
   temp     = float(temp)
   dens_l   = float(dens_l)
   times    = np.array(times).astype(float)
   henry    = np.array(henry).astype(float)

   # Temperature and Henry's law constant
   TK   = temp + 273.15
   kh   = henry[0] * np.exp(henry[1] * (1/TK - 1/298.15))       # mol/kg-bar as liq:gas  (vant hoff)
   kh   = kh * dens_l / 1000                                    # mol/L-bar
   Kaw  = 1 / (kh * R * TK)                                     # dimensionless, gas:liq, e.g., g/L / g/L or g/m3 per g/m3
   
   # Temperature dependent constants
   KW        = 10**(-4.2195 - 2915.16/TK)                                                             # Water dissociation constant 
   rho_water = 1000 * (1 - ((temp + 288.9414) / (508929.2 * (temp + 68.12963))) * (temp - 3.9863)**2)   # density of water  [kg/m3]
   K2        = np.exp(-12092.1/TK -36.786 * np.log(TK) + 235.482) * rho_water                           # equlibrium constant for CO2 + H2O -> HCO3^- + H^+ [mol/m3]
   K4        = K2/ (KW * 1e6)                                                                                  # equlibrium constant for CO2 + OH^- --> HCO3^-
   # equilibrium constants from https://github.com/sashahafner/NH3-RTM kinSpec()
   K_hco3    = 10**(-353.5305 - 0.06092*TK + 21834.37/TK + 126.8339 * np.log10(TK) - 1684915/TK**2)
   K_co3     = 10 **(-461.4176 - 0.093448*TK + 26986.16/TK + 165.7595*np.log10(TK) - 2248629/TK**2)


   # Reaction rate constants
   # refrence:  Rate-Based Modeling of Reactive Absorption of CO2 and H2S into Aqueous Methyldiethanolamine
   #            Comprehensive Study of the Hydration and Dehydration Reactions of Carbon Dioxide in Aqueous Solution
   k3 = 4.315 * 10**13 * np.exp(- 6666/TK) / 1000               # Second-order rate constant for CO2 + OH⁻ → HCO3⁻  [m3/mol-s]
   k1 = 6.672 * 10**12 * np.exp(- 9.724 * 10**3 / TK)           # First-order rate constant for CO2 + H2O → H2CO3   [s^-1]
   k2 = 9.000 * 10**13 * np.exp(- 8.612 * 10**3 / TK)           # First-order rate constant for H2CO3  → CO2 + H2O  [s^-1]
   
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
   cg0      = gpm3_to_molpm3(cg0      , M_co2)
   cgin     = gpm3_to_molpm3(cgin     , M_co2)
   cl_co20  = gpm3_to_molpm3(cl_co20  , M_co2)
   clin_co2 = gpm3_to_molpm3(clin_co2 , M_co2)


   #============================== mol balance coming in ===============================
   Q_g = v_g * 1     # m/s * 1m^2 = m3/s
   m_gin = Q_g * cgin  # m3/s * mol/m3  = mol/s





   cr = np.array([
        gpm3_to_molpm3(cr_co20,  M_co2),
        cr_TOTC0,
    ]) * v_res  # mol

   # Uniform initial concentrations throughout reactor
   ccg    = np.full((nc), cg0)
   c_co2  = np.full((nc), cl_co20)   
   c_TOTC = np.full((nc), cl_TOTC0) 

   # convert concentration to moles
   ncg      = ccg    * vol_gas            # mol/m3 * m3 = mol
   ncl_co2  = c_co2  * vol_liq
   ncl_TOTC = c_TOTC * vol_liq


   # Initial state variable array
   y0 = np.concatenate([ncg, ncl_co2, ncl_TOTC, cr])
   #  y0 = np.append (y0, ncr)


   
   # Solve/integrate
   out = solve_ivp(rates, [0, max(times)], y0 = y0, 
                   t_eval = times, 
                   args = (v_g, v_l, cgin, clin_co2, clin_TOTC,
                           vol_gas, vol_liq, vol_tot,
                           k1, k2, k3, K4, K_hco3, K_co3, KW, ex_oh, Kga, v_res,
                           temp, henry, pres, ssa, dens_l, por_l, por_g, counter, recirc,
                           enh_method),
                   method = 'Radau')
   
   # Extract moles of compound [position, time]
   ncgt      = out.y[0        : nc]       
   ncl_co2t  = out.y[nc       : (2 * nc)]
   ncl_TOTCt = out.y[(2 * nc) : (3 * nc)]
   ncr       = out.y[(3 * nc) : (3 * nc + 2)]
   nctot     = np.sum(ncgt,0) + np.sum(ncl_co2t,0) + np.sum(ncl_TOTCt,0)  #total mass of compuond in the entire column as a function of time. why? 

   # Get concentrations vs. time
   # Gas
   ccgt      = ncgt / np.transpose(np.tile(vol_gas, (ncgt.shape[1], 1)))
   # co2
   ccl_co2t  = ncl_co2t / np.transpose(np.tile(vol_liq, (ncl_co2t.shape[1], 1)))
   # TOTC
   ccl_TOTCt = ncl_TOTCt / np.transpose(np.tile(vol_liq, (ncl_TOTCt.shape[1], 1)))
   
   # ====== equilibrium conc.=======
   ccl_co2t_eq = ccgt / Kaw
   # ================== mol balance out ================
   m_gout = Q_g * ccgt[-1,:]
   
   Q_l    = v_l * 1    # m3/s
   m_lout = Q_l * (ccl_co2t[-1,:] + ccl_TOTCt[-1,:])

   m_tout = m_gout + m_lout 


   # Total
  #  cctt = nctot / np.transpose(np.tile(vol_tot, (nctot.shape[1], 1)))

   nct = np.concatenate([ncgt, ncl_co2t,ncl_TOTCt])

   # integrating the speciation model
   pH_profile = np.zeros_like(ccl_TOTCt)
   for c in range(ccl_TOTCt.shape[0]):
       for t in range(ccl_TOTCt.shape[1]):
           res_out = sm.spec2_matrix(ccl_TOTCt[c,t], K_hco3, K_co3, KW, ex_oh)
           pH_profile[c,t] = res_out['pH']
    
   # ===== TOTC equilibrium ==== 
   c_h_profile  = 10**(-pH_profile) * 1000
   c_h2co3_eq   = k1/k2 * ccl_co2t_eq
   c_hco3_eq    = K_hco3 * 1000 * c_h2co3_eq /c_h_profile
   c_co3_eq     = K_co3  * 1e6 * c_h2co3_eq/c_h_profile**2
   TOTC_eq      = c_h2co3_eq + c_hco3_eq  + c_co3_eq

   # Return results as a dictionary
   return {'gas_conc'    : molpm3_to_gpm3(ccgt,M_co2),
          'co2_liq_conc' : molpm3_to_gpm3(ccl_co2t,M_co2),
          'pH_profile'   : pH_profile,
          'TOTC_liq_conc': ccl_TOTCt,
          'cell_pos'     : x,
          'time'         : times,
          'inputs'       : args_in, 
          'm_gin'        : m_gin,
          'm_gout'       : m_gout,
          'm_lout'       : m_lout,
          'm_tout'       : m_tout,
          'eq_conc'      : molpm3_to_gpm3(ccl_co2t_eq, M_co2),
          'TOTC_eq'      : TOTC_eq,  
          'pars'         : {'gas_rt': rt_gas, 'liq_rt': rt_liq, 'Kga': Kga, 'Kaw': Kaw, 'ae':ae,'kg':kg,'kl':kl}}