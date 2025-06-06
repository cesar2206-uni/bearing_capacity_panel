import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Nq_factor(phi):
    """
    Nq factor of the bearing capacity theory
    """
    phi_rad = np.radians(phi)
    return np.tan(np.pi/4 + phi_rad/2)**2 * np.exp(np.pi * np.tan(phi_rad))

def Nc_factor(phi):
    """
    Nc factor of the bearing capacity theory
    """
    N_q = Nq_factor(phi)
    phi_rad = np.radians(phi)
    return (N_q - 1) / np.tan(phi_rad)

def Ng_factor(phi, case = "rough"):
    """
    Ng factor of the bearing capacity theory
    In this case exist 2 types, smooth and rough surface. In general, the rough
    surface is used  
    """
    phi_rad = np.radians(phi)
    if case == "smooth":
        N_g = 0.0663 * np.exp(0.1623 * phi)
    elif case == "rough":
        N_g = 0.1054 * np.exp(0.1675 * phi)
    return N_g

def foundation_shape_Vesic1975(B_prime, L_prime, phi):
    """
    Modification factors for general bearing capacity equation (Vesic 1975)
    Equations for Foundation Shape
    """
    N_q = Nq_factor(phi)
    N_c = Nc_factor(phi)
    
    # Load Factor
    phi_rad = np.radians(phi)
    S_qs =  1 + (B_prime / L_prime) * np.tan(phi_rad)
    
    # Cohesion Factor 
    S_cs = 1 + (B_prime / L_prime) * (N_q / N_c)
        
    # Unit Weight Factor
    S_gs = 1 - 0.4 * (B_prime / L_prime)
    
    return S_cs, S_qs, S_gs

def inclined_loading_Vesic1975(B, L, B_prime, L_prime, horizontal_force, vertical_force, phi, cohesion, theta_direction = 90):
    """
    Modification factors for general bearing capacity equation (Vesic 1975)
    Equations for Inclined Loading
    """
    
    # Previous calculations 
    theta_rad = np.radians(theta_direction)
    phi_rad = np.radians(phi)
    m_B = (2 + B/L) / (1 + B/L)
    m_L = (2 + L/B) / (1 + L/B)
    N_c = Nc_factor(phi)
    
    if horizontal_force or vertical_force is False:
        return 1, 1, 1
    
    if theta_direction == 90:
        m = m_B
    elif theta_direction == 0:
        m = m_L
    else:
        m = m_L * (np.cos(theta_rad) ** 2) + m_B * (np.sin(theta_rad) ** 2)
    
    # Load Factor
    if phi == 0:
        S_qi = 1
    else:
        S_qi = (1 - horizontal_force / (vertical_force + B_prime * L_prime * cohesion / np.tan(phi_rad))) ** (m)
    
    # Cohesion Factor
    if phi == 0:
        S_ci = 1 - (m * horizontal_force) / (B_prime * L_prime * cohesion * N_c)
    else:
        S_ci = S_qi - (1 - S_qi) / (N_c * np.tan(phi_rad))
        
    # Unit Weight Factor
    S_gi = (1 - horizontal_force / (vertical_force + B_prime * L_prime * cohesion / np.tan(phi_rad))) ** (m + 1)
    
    return S_ci, S_qi, S_gi

def foundation_depth_Vesic1975(B, Df, phi):
    """
    Modification factors for general bearing capacity equation (Vesic 1975)
    Equations for Foundation Depth
    """
    
    # Previous calculations
    if Df / B <= 1:
        k = Df / B
    else:
        k = np.tan(Df / B)
    N_c = Nc_factor(phi)
    phi_rad = np.radians(phi)
        
    # Load Factor
    S_qd  = 1 + 2 * np.tan(phi_rad) * ((1 - np.sin(phi_rad)) ** 2) * k
    
    # Cohesion Factor
    if phi == 0:
        S_cd = 1 + 0.4 * k
    else:
        S_cd = S_qd - (1 - S_qd)/(N_c * np.tan(phi_rad))
        
    # Unit Weight Factor
    S_gd = 1
    
    return S_cd, S_qd, S_gd

def surface_slope_Vesic1975(beta, phi):
    """
    Modification factors for general bearing capacity equation (Vesic 1975)
    Equations for Surface Slope
    """
    
    # Previous calculations
    beta_rad = np.radians(beta)
    phi_rad = np.radians(phi)    
    N_c = Nc_factor(phi)
    
    # Load Factor
    S_qb = (1 - np.tan(beta_rad)) ** 2
    
    # Cohesion Factor
    if phi == 0:
        S_cb = 1 - 2 * beta_rad / (np.pi + 2)
    else:
        S_cb = S_qb - (1 - S_qb) / (N_c * np.tan(phi_rad))
    
    # Unit Weight Factor
    if phi == 0:
        S_gb = - 2 * np.sin(beta_rad)
    else: 
        S_gb = S_qb
        
    return S_cb, S_qb, S_gb

def base_inclination_Vesic1975(delta, phi):
    """
    Modification factors for general bearing capacity equation (Vesic 1975)
    Equations for base inclination
    """
    # Previous calculations
    delta_rad = np.radians(delta)
    phi_rad = np.radians(phi)
    N_c = Nc_factor(phi)
    
    # Load Factor
    S_qdl = (1 - delta_rad * np.tan(phi_rad)) ** 2
    
    # Cohesion Factor
    if phi == 0:
        S_cdl = 1 - 2 * delta_rad / (np.pi + 2)
    else:
        S_cdl = S_qdl - (1 - S_qdl) / (N_c * np.tan(phi_rad))
    
    # Unit Weight Factor
    S_gdl = S_qdl
    
    return S_cdl, S_qdl, S_gdl

def water_table_correction(gamma, B, Df, water_table = False):
    """
    Corrects unit weight (gamma) and surcharge (q_s) based on water table level.
    """
    # Constants (assumed)
    gamma_w = 9.81  # Unit weight of water [kN/m³]
    gamma_sub = gamma - gamma_w  # Submerged unit weight

    # Case 1: Water table below D + B (no effect)
    if water_table is False or water_table > Df + B:
        corrected_qs = gamma * Df
        corrected_gamma = gamma
    
    # Case 2: Water table at surface (z = 0)
    elif water_table == 0:
        corrected_qs = gamma_sub * Df
        corrected_gamma = gamma_sub

    # Case 3: Water table at foundation level (z = D)
    elif water_table == Df:
        corrected_qs = gamma * Df
        corrected_gamma = gamma_sub

    # Case 4: Water table between D and D + B
    elif Df < water_table < Df + B:
        z = water_table
        gamma_effective = gamma_sub + ((z - Df) / B) * (gamma - gamma_sub)
        corrected_qs = gamma * Df
        corrected_gamma = gamma_effective

    # Case 5: Water table above D (0 < z < D)
    elif 0 < water_table < Df:
        # Split the surcharge into two parts:
        # Submerged part from 0 to z, saturated from z to D
        corrected_qs = gamma_sub * water_table + gamma * (Df - water_table)
        corrected_gamma = gamma_sub

    else:
        raise ValueError("Invalid water_table depth provided.")

    return corrected_qs, corrected_gamma

def allowed_bearing_capacity_Canadian2024(B, L, B_prime, L_prime, Df, beta, delta, horizontal_force, vertical_force,
                                          cohesion, phi, gamma, FS, water_table = False, theta_direction = 90):
    
    # Bearing capacity factors
    N_q = Nq_factor(phi)
    N_c = Nc_factor(phi)
    N_g = Ng_factor(phi, case = "rough")
    
    # Correction Factors
    S_cs, S_qs, S_gs = foundation_shape_Vesic1975(B_prime, L_prime, phi)
    S_ci, S_qi, S_gi = inclined_loading_Vesic1975(B, L, B_prime, L_prime,
                               horizontal_force, vertical_force, 
                               phi, cohesion, theta_direction)
    S_cd, S_qd, S_gd = foundation_depth_Vesic1975(B, Df, phi)
    S_cb, S_qb, S_gb = surface_slope_Vesic1975(beta, phi)
    S_cdl, S_qdl, S_gdl = base_inclination_Vesic1975(delta, phi)
    
    S_c = S_cs * S_ci * S_cd * S_cb * S_cdl
    S_q = S_qs * S_qi * S_qd * S_qb * S_qdl
    S_g = S_gs * S_gi * S_gd * S_gb * S_gdl
    
    # Previous calculation
    corrected_qs, corrected_gamma = water_table_correction(gamma, B, Df, water_table)
    
    # Ultimate bearing capacity
    q_u = cohesion * N_c * S_c + corrected_qs * N_q * S_q + 0.5 * corrected_gamma * B_prime * N_g * S_g

    # Allowed bearing capacity
    q_a = q_u / FS
    
    # Print factors
    # correction_factors = pd.DataFrame({
    # "Correction": ["Shape", "Inclination", "Depth", "Slope", "Base Inclination"],
    # "S_c": [S_cs, S_ci, S_cd, S_cb, S_cdl],
    # "S_q": [S_qs, S_qi, S_qd, S_qb, S_qdl],
    # "S_γ": [S_gs, S_gi, S_gd, S_gb, S_gdl]
    # })
    # print(correction_factors)
    
    # bearing_factors = pd.DataFrame({
    # "Factor": ["N_c", "N_q", "N_γ"],
    # "Value": [N_c, N_q, N_g]
    # })
    # print(bearing_factors)
    
    return q_a

def barplot_bearing_capacity_Canadian2024(L_over_B_values, B_values, q_applied, Df,  beta, delta, horizontal_force, vertical_force,
                                          cohesion, phi, gamma, FS, water_table, theta_direction):
    """
    Obtain a Barplot of Bearing capacity
    """
    
    # Declaring the maximum size of the vector
    hatch_patterns = ['//', '\\\\', '||', '--', '++', 'xx', 'oo', '..']
    
    if len(L_over_B_values) > 8 or len(B_values) > 8:
        raise ValueError("Invalid water_table depth provided.")
    
    # Creation of barplots
    
    fig, ax = plt.subplots(figsize = (12, 6))
    
    x = np.arange(len(B_values))
    bar_width = 0.1
        
    for i, L_over_B in enumerate(L_over_B_values):
        q_allows = [allowed_bearing_capacity_Canadian2024(B_value, L_over_B * B_value, B_value, 
                                                          L_over_B * B_value, Df, 
                                                          beta, delta, horizontal_force, vertical_force,
                                                          cohesion, phi, gamma, FS, water_table, theta_direction) for B_value in B_values]
        bars = ax.bar(x + i * bar_width, q_allows, width = bar_width, label=f'L/B = {L_over_B}',
                      color='white', edgecolor='black', hatch=hatch_patterns[i % len(hatch_patterns)])    
    
    
    
    # Formatting
    x_labels = [str(B) for B in B_values]
    ax.set_xlabel("B Values")
    ax.set_ylabel("Allowable Bearing Capacity (kPa)")
    ax.set_title("Allowable Bearing Capacity vs B for Various L/B ratios (Canadian Manual, 2024)")
    ax.set_xticks(x + bar_width * (len(L_over_B_values) - 1) / 2)
    ax.set_xticklabels(x_labels)
    ax.axhline(y=q_applied, color='red', linestyle='--', linewidth=2, label='Applied Pressure')
    ax.legend(title="L/B Ratio", bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    plt.show()