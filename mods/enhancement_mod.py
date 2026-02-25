import numpy as np

def enh_fac(c_oh, c_co2, k, K, kl, method ):
    """
    
    Calculating the CO2 mass-transfer enhancement
    factor for a reacting gas-liquid system,
    using either a pseudo first-order (PFO)
    or a fast second-order reversible reaction (RSO) 

    Author: Osman Alsheghri


    Parameters/Inputs
    -----------------
    c_oh   : float or ndarray
        OH^- concentration in the bulk liquid           [mol/m3]

    c_co2  : float or ndarray
        CO2 concentration at the interface              [mol/m3]

    k      : float
        Second order rate constant                      [m3/(mol*s)]
    
    K      : float  
        Equilibrium constant for CO2 + OH^- ⇌ HCO3^-
    
    kl     : float
        Liquid-side mass transfer coefficient            [m/s]
    
    method : 'PFO' or 'RSO'
        'PFO' for pseudo-first order reaction enhancement
        'RSO' for second-order reversible reaction enhancement
    
    Returns
    -------
    E      : float or ndarray
        Enhancement factor with the same shape as 
        c_oh and c_co2

    """
    D1 = 1.89 * 10 ** -9                    # Diffusivity of CO2    [m2/s]
    D2 = 5.40 * 10 ** -9                    # Diffusivity of OH^-   [m2/s]
    D3 = 1.20 * 10 ** -9                    # Diffusivity of HCO3^- [m2/s]
    if method == 'PFO':
        k_pfo = k * c_oh                    # pseudo first order rate constant  [1/s]
        Ha = np.sqrt(k_pfo * D1) / kl       # Hatta number
        E  = Ha / np.tanh(Ha)               # Enhancement factor for pfo reaction
    
    elif method == 'RSO':
        denom = 1.0 + K * (D3 / D2) * c_co2
        if np.any(denom <= 0.0):
            raise ValueError('Denominator in RSO expression is non-positive or 0')
        E  = 1.0 + (D3 / D1) * (K * c_oh / denom)
    else : 
        raise ValueError("method must be 'PFO' or 'RSO'.")
    return E