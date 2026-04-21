import numpy as np
import pandas as pd
import sys
from scipy.integrate import solve_ivp
import math

# Rates function
# Arg order: time, state variable, then arguments
# For CO2 absorption in aqueous NaOH solution at high pH (>10)
def rates(t, mc, v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, pH, Kga, v_res, k2, temp, henry, pres, pKa, ssa, dens_l, por_l, por_g, counter=True, recirc=False):
    # Interpolation function for time-variable inputs
    def interpolation(t, c):
        if isinstance(c, pd.DataFrame):
            c = np.interp(t, c.iloc[:, 0], c.iloc[:, 1])
        elif isinstance(c, np.ndarray):
            c = np.interp(t, c[0], c[1])
        elif not (isinstance(c, (int, float))):
            sys.exit(f'Error: input must be float, integer, numpy array, or pandas DataFrame, but is {type(c)}')
        return c

    # Interpolate time-variable concentrations/inputs at current time t
    cgin = interpolation(t, cgin)
    clin = interpolation(t, clin)
    pH   = interpolation(t, pH)

    # Reaction rate setup (k2 default = k for compatibility)
    if isinstance(k2, str) and k2.lower() == 'default':
        k2 = k

    # Calculate OH concentration from pH for CO2 + OH⁻ kinetics
    pOH   = 14 - pH
    M_OH  = 17.007  # g/mol, molar mass of OH⁻
    c_OH  = 10**(-pOH) * M_OH * 1000  # g/m³, OH⁻ concentration
    k_eff = k * c_OH  # 1/s, pseudo-first-order rate constant

    # Calculate Kga if 'onda' requested (Onda correlations)
    if isinstance(Kga, str) and Kga.lower() == 'onda':
        # Hard-coded constants
        g      = 9.81  # m/s²
        Dg     = 1.16e-5  # m²/s, gas diffusion coefficient
        Dliq   = 1.89e-9  # m²/s, liquid diffusion coefficient
        sigm_c = 0.072  # N/m, critical surface tension
        sigm_l = 0.072  # N/m, surface tension
        R      = 0.083144  # L·bar/(K·mol)
        mw_g   = 28.97  # g/mol, air molar mass

        dp = 6 * (1 - por_g) / ssa  # m, characteristic packing length
        # Empirical particle diameter
        dp_emp = 2.0 if dp < 15 else 5.23  # m

        TK     = temp + 273.15  # K
        dens_g = pres * mw_g / (R * TK)  # kg/m³, gas density
        visc_g = 9.1e-8 * TK - 1.16e-5  # Pa·s, gas viscosity
        visc_l = -2.55e-5 * TK + 8.51e-3  # Pa·s, liquid viscosity

        # Dimensionless numbers
        Re = dens_l * v_l / (ssa * visc_l)  # Reynolds
        Fr = v_l**2 * ssa / g  # Froude
        We = v_l**2 * dens_l / (sigm_l * ssa)  # Weber

        # Effective area
        ae = ssa * (1.0 - 2.71828**(-1.45 * (sigm_c / sigm_l)**0.75 *
                                    Re**0.1 * Fr**-0.05 * We**0.2))

        # Gas phase coefficient
        kg = (dp_emp * (v_g * dens_g / (ssa * visc_g))**0.7 *
              (visc_g / (dens_g * Dg))**(1/3) * (ssa * dp)**-2 * ssa * Dg)

        # Liquid phase coefficient
        kl = (0.0051 * (v_l * dens_l / (ae * visc_l))**(2/3) *
              (visc_l / (dens_l * Dliq))**(-0.5) * (ssa * dp)**0.4 *
              (dens_l / (visc_l * g))**(-1/3))

        # Henry's law (van't Hoff)
        kh  = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15))  # mol/(kg·bar)
        kh *= dens_l / 1000  # mol/(L·bar)
        Kaw = 1 / (kh * R * TK)  # dimensionless air-water partition (high pH: all CO2 reacts)
        Daw = Kaw

        # Hatta enhancement (fast reaction, Kasper & Feilberg 2019)
        Ha = np.sqrt(Dliq * k_eff) / kl
        E  = Ha / np.tanh(Ha)

        # Overall resistance (series, Seader et al. p181)
        Rtot = 1 / (kg * ae) + Daw / (kl * ae * E) + Daw / (k_eff * por_l)
        Kga = 1 / Rtot  # s⁻¹

    # Number of cells (integer division)
    nc = mc.shape[0] // 2

    # Separate state variables: gas CO2 (mcg), liquid CO2 (mcl), reservoir CO2 (mcr) [g/m²]
    mcg = mc[0:nc]
    mcl = mc[nc:2*nc]
    mcr = mc[2*nc:2*nc+1]

    # Concentrations [g/m³]
    ccg = mcg / vol_gas
    ccl = mcl / vol_liq

    dmcr = 0.0

    # Reservoir handling (if recirc=True)
    if recirc:
        oi = 0 if counter else nc - 1
        if v_res > 0:
            cr = mcr / v_res
            rxnr = k_eff * mcr  # g/(s·m²)
            dmcr = (ccl[oi] - cr) * abs(v_l) - rxnr  # g/(s·m²)
        elif v_res == 0:
            clin = ccl[oi]
        else:
            sys.exit('v_res must be positive float or 0')

    # Mass transfer rate [g/s]
    g2l = Kga * vol_tot * (ccg - ccl * Daw)

    # Gas phase derivatives [g/s] (advection only)
    cvec_gas  = np.insert(ccg, 0, cgin)
    advec_gas = -v_g * np.diff(cvec_gas)
    dmg       = advec_gas - g2l

    # Liquid phase derivatives [g/s] (advection + transfer - reaction)
    rxn = k_eff * mcl  # g/(s·m²)
    if counter:
        v_l = -v_l
        cvec_liq = np.insert(ccl, nc, clin)
    else:
        cvec_liq = np.insert(ccl, 0, clin)
    advec_liq = -v_l * np.diff(cvec_liq)
    dml = advec_liq + g2l - rxn

    # Combined derivatives
    dm = np.concatenate([dmg, dml])
    dm = np.append(dm, dmcr)
    return dm

# Main model function (TBR simulation)
def tfmod(L, por_g, por_l, v_g, v_l, nc, cg0, cl0, cgin, clin, k, Kga, henry, pKa, pH, temp, dens_l, times, kg='onda', kl='onda', ae='onda', v_res=0, k2='default', ccr=0, pres=1., ssa=1100, typ='TBD', counter=True, recirc=False):
    """
    Simulate CO2 absorption in trickle bed reactor with NaOH solution (high pH).

    Parameters
    ----------
    L : float
        Reactor height (m)
    por_g, por_l : float
        Gas/liquid porosity (m³ phase/m³ total)
    v_g, v_l : float
        Superficial gas/liquid velocities (m/s)
    nc : int
        Number of axial cells
    cg0, cl0 : float
        Initial gas/liquid CO2 (g/m³)
    cgin : float/array/DataFrame
        Gas inlet CO2 (g/m³)
    clin : float/array/DataFrame
        Liquid inlet CO2 (g/m³)
    k : float
        2nd-order rate constant CO2+OH⁻ (m³/(g·s))
    Kga : float/str
        Overall MTC (s⁻¹) or 'onda'
    henry : array
        [kH(25°C), dlnkH/d(1/T)] from NIST
    pKa : float
        Not used (high pH)
    pH : float/array/DataFrame
        Inlet/bulk pH (time-varying OK)
    temp : float
        Temperature (°C)
    dens_l : float
        Liquid density (kg/m³)
    times : array
        Output times (s)
    v_res : float
        Reservoir volume per m² (m³, 0=no recirc)
    k2 : str/float
        'default' → k2=k (compatibility)
    ccr : float
        Initial reservoir CO2 (g/m³)
    counter : bool
        True=countercurrent
    recirc : bool
        True=liquid recirculation

    Returns
    -------
    dict : Results with profiles (pos, time), masses, parameters
    """
    # Save inputs
    args_in = locals()
    R       = 0.083144  # L·bar/(K·mol)

    # Retention times [s]
    rt_gas = L * por_g / v_g
    rt_liq = L * por_l / v_l

    # Ensure numeric
    cg0, cl0, k, pKa, temp, dens_l = map(float, [cg0, cl0, k, pKa, temp, dens_l])
    times = np.array(times).astype(float)
    henry = np.array(henry).astype(float)

    # Henry's law [dimensionless gas:liq]
    TK  = temp + 273.15
    kh  = henry[0] * math.exp(henry[1] * (1/TK - 1/298.15))
    kh *= dens_l / 1000
    Kaw = 1 / (kh * R * TK)

    # Discretization
    x  = np.linspace(0, L, nc + 1)
    dx = np.diff(x)[0]
    x  = x[1:nc+1] - dx/2  # Cell centers [m]

    vol_tot = dx * 1  # m³/cell (A=1 m²)
    vol_gas = vol_tot * por_g
    vol_liq = vol_tot * por_l

    # Initial masses [g/m²]
    mcg0 = np.full(nc, cg0) * vol_gas
    mcl0 = np.full(nc, cl0) * vol_liq
    mcr0 = ccr * v_res
    y0   = np.concatenate([mcg0, mcl0, [mcr0]])

    # Integrate
    out = solve_ivp(rates, [0, max(times)], y0, t_eval=times,
                    args=(v_g, v_l, cgin, clin, vol_gas, vol_liq, vol_tot, k, pH, Kga, v_res, k2, temp, henry, pres, pKa, ssa, dens_l, por_l, por_g, counter, recirc),
                    method='Radau')

    # Extract/post-process [pos, time]
    nc    = out.y.shape[0] // 2  # Recompute for slicing
    mcgt  = out.y[0:nc]
    mclt  = out.y[nc:2*nc]
    mcrt  = out.y[2*nc:2*nc+1]
    mctot = np.sum(mcgt, 0) + np.sum(mclt, 0)  # Total column mass [g]

    ccgt = mcgt / np.transpose(np.tile(vol_gas, (mcgt.shape[1], 1)))  # g/m³ gas
    cclt = mclt / np.transpose(np.tile(vol_liq, (mclt.shape[1], 1)))  # g/m³ liq

    return {
        'gas_conc': ccgt, 'liq_conc': cclt, 'gas_mass': mcgt, 'liq_mass': mclt,
        'reservoir_mass': mcrt, 'tot_mass': mctot, 'cell_pos': x, 'time': times,
        'inputs': args_in,
        'pars': {'gas_rt': rt_gas, 'liq_rt': rt_liq, 'Kga': Kga, 'Kaw': Kaw}
    }