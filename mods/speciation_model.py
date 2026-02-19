import numpy as np
from scipy.optimize import root_scalar
# here is a function with same form as Sasha's but without the matrix...
def spec(TOTC,K_hco3,K_co3,KW):
    """
    Calculate carbonate system speciation and pH from total carbonate
    concentration

    Author: Osman Alsheghri

    Parameters/Input
    -------------
    TOTC  : float
            Total carbonate concentration (mol/m3) = [H2CO3] + [HCO3^-] + [CO3^2-]

    K_hco3: float
            Equilibrium constant for H2CO3 ⇌ HCO3^- + H^+

    K_co3 : float
            Equilibrium constant for HCO3^- ⇌ CO3 ^2-  + H^+

    KW    : float
            Water dissociation constant 

    Returns
    ----------
    Dictionary
        'c_h2co3' : float - H2CO3 concentration     (mol/m3)
        'c_hco3'  : float - HCO3^- concentration    (mol/m3)
        'c_co3'   : float - CO3^2- concentration    (mol/m3)
        'c_oh'    : float - OH^- concentration      (mol/m3)
        'c_h'     : float - H^+ concentration       (mol/m3)
        'pH'      : float - pH value
    """
    def residC(c_h2co3):
        def residH(c_h):
            c_hco3  = (K_hco3*c_h2co3)/c_h
            c_co3   = (K_hco3*K_co3*c_h2co3)/c_h**2
            c_oh    = KW/c_h
            return c_h - c_hco3 - 2*c_co3-c_oh
        c_h     = root_scalar(residH, bracket=[1E-14, 1], method='brentq').root
        c_hco3  = (K_hco3*c_h2co3)/c_h
        c_co3   = (K_hco3*K_co3*c_h2co3)/c_h**2
        c_oh    = KW/c_h
        return TOTC - c_h2co3 - c_hco3 - c_co3, c_h, c_hco3,c_co3,c_oh
    c_h2co3 = root_scalar(lambda x: residC(x)[0], bracket=[1E-10, TOTC], method='brentq').root
    residc, c_h, c_hco3,c_co3,c_oh = residC(c_h2co3)
    pH = - np.log10(c_h)
    out = {
        'c_h2co3': c_h2co3,
        'c_hco3' : c_hco3,
        'c_co3'  : c_co3,
        'c_oh'   : c_oh,
        'c_h'    : c_h,
        'pH'     : pH

    }
    return out








# Sasha said maybe we could find an algerbraic equation for c_h2co3 instead
# of solving it numerically.
# Maybe like this? where we use the TOTC equation to isolate c_h2co3, then 
# we plug it in the charge balance equation, this makes everything depend only of c_h??
# I actually did exactly what has been done in Sasha's R model to handle algebraic c_h2co3 :) 
def spec2(TOTC,K_hco3,K_co3,KW):
    """
    Calculate carbonate system speciation and pH from total carbonate
    concentration

    Author: Osman Alsheghri

    Parameters/Input
    -------------
    TOTC  : float
            Total carbonate concentration (mol/m3) = [H2CO3] + [HCO3^-] + [CO3^2-]

    K_hco3: float
            Equilibrium constant for H2CO3 ⇌ HCO3^- + H^+

    K_co3 : float
            Equilibrium constant for HCO3^- ⇌ CO3 ^2-  + H^+

    KW    : float
            Water dissociation constant 

    Returns
    ----------
    Dictionary
        'c_h2co3' : float - H2CO3 concentration     (mol/m3)
        'c_hco3'  : float - HCO3^- concentration    (mol/m3)
        'c_co3'   : float - CO3^2- concentration    (mol/m3)
        'c_oh'    : float - OH^- concentration      (mol/m3)
        'c_h'     : float - H^+ concentration       (mol/m3)
        'pH'      : float - pH value
    """
    def residH(c_h):
            c_h2co3 = TOTC/(1 + K_hco3/c_h + (K_co3*K_hco3)/c_h**2)
            c_hco3  = (K_hco3*c_h2co3)/c_h
            c_co3   = (K_hco3*K_co3*c_h2co3)/c_h**2
            c_oh    = KW/c_h
            return c_h -  c_hco3 - 2*c_co3-c_oh
    c_h     = root_scalar(residH, bracket=[1E-14, 1], method='brentq').root
    c_h2co3 = TOTC/(1 + K_hco3/c_h + (K_co3*K_hco3)/c_h**2)
    c_hco3  = (K_hco3*c_h2co3)/c_h
    c_co3   = (K_hco3*K_co3*c_h2co3)/c_h**2
    c_oh    = KW/c_h
    pH      = - np.log10(c_h) 
    out = {
        'c_h2co3': c_h2co3,
        'c_hco3' : c_hco3,
        'c_co3'  : c_co3,
        'c_oh'   : c_oh,
        'c_h'    : c_h,
        'pH'     : pH 
    }
    return out









# the same as spec2 but here with matrix apporach. When we are 100% sure it works we could just
# remove the two above or have a module with matrix approach only.  

def spec2_matrix(TOTC,K_hco3,K_co3,KW):
     """
    Calculate carbonate system speciation and pH from total carbonate
    concentration

    Author: Osman Alsheghri

    Parameters/Input
    -------------
    TOTC  : float
            Total carbonate concentration (mol/m3) = [H2CO3] + [HCO3^-] + [CO3^2-]

    K_hco3: float
            Equilibrium constant for H2CO3 ⇌ HCO3^- + H^+

    K_co3 : float
            Equilibrium constant for HCO3^- ⇌ CO3 ^2-  + H^+

    KW    : float
            Water dissociation constant 

    Returns
    ----------
    Dictionary
        'c_h2co3' : float - H2CO3 concentration     (mol/m3)
        'c_hco3'  : float - HCO3^- concentration    (mol/m3)
        'c_co3'   : float - CO3^2- concentration    (mol/m3)
        'c_oh'    : float - OH^- concentration      (mol/m3)
        'c_h'     : float - H^+ concentration       (mol/m3)
        'pH'      : float - pH value
     """
    # From the equilibriums we know that 
    # HCO3 = K1 * c_h2co3**1 * c_h**-1 
    # CO3  = K1*K2 * c_h2co3**1 * c_h **-2 
    # OH   = KW * c_h **-1 

    # define the stoichiometry matrix
    #        | H2CO3 |  H+
    #        --------|----
    # HCO3-  |  -1   |  1
    # CO3-2  |  -1   |  2
    #   OH-  |   0   |  1
     TOTC = TOTC/1000
     S = np.array([
          [-1 , 1], # HCO3^-
          [-1 , 2], # CO3^2-
          [0 , 1]  # OH^-
     ])
     K = np.array([K_hco3,K_co3,KW])
     def residH(c_h):
          # denominator for c_h2co3 equation --> denom = 1 + K1/c_h + K1*K2/c_h**2  
          denom = 1 + np.sum(K * (c_h**(-S[:,1])) * (- S[:,0])) # if we did not have the last part
          # then we would also say that denom depends on OH, it does not.
          c_h2co3 = TOTC / denom
          conc = K * (c_h2co3 ** (-S[:,0])) * (c_h ** (-S[:,1])) # conc is now an array with [c_hco3, c_co3, c_oh]
          return c_h - np.sum(conc * S[:,1])
     # find the concentration of H+
     c_h    = root_scalar(residH, bracket=[1E-14, 1], method='brentq').root
     # calculate the concentration of the equilibrium species
     denom   = 1 + np.sum(K * (c_h**(-S[:,1])) * (- S[:,0]))
     c_h2co3 = TOTC / denom
     conc    = K * (c_h2co3 ** (-S[:,0])) * (c_h ** (-S[:,1]))
     pH = - np.log10(c_h)
     out = {
        'c_h2co3': c_h2co3 * 1000,
        'c_hco3' : conc[0] * 1000,
        'c_co3'  : conc[1] * 1000,
        'c_oh'   : conc[2] * 1000,
        'c_h'    : c_h     * 1000,
        'pH'     : pH 
    }
     return out

