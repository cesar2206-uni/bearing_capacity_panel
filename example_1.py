from functions import *

B = 2  # Width
L = 2 # Length
B_prime = B  # Width corrected by eccentricity
L_prime = L  # Width corrected by eccentricity
Df = 0  # Depth
beta = 0  # Ground slope
delta = 0  # Delta tilt
horizontal_force = False  # No horizontal force
vertical_force = False  # No vertical force
cohesion = 0  # Apparent cohesion
phi = 40  # Friction angle
gamma = 23  # Unit Weight (KN/m3)
FS = 3  # Factor of Safety (Normally is 3)
water_table = False  # No water table
theta_direction = 90  # Force inclined in B direction

q_a = allowed_bearing_capacity_Canadian2024(B, L, B_prime, L_prime, Df, beta, delta, horizontal_force, vertical_force,
                                          cohesion, phi, gamma, FS, water_table, theta_direction)

print(q_a * 0.0101972)


L_over_B_values=[1, 2, 3, 4, 8, 12, 16]
B_values=[2, 3, 4, 5, 6]
q_applied = 200
barplot_bearing_capacity_Canadian2024(L_over_B_values, B_values, q_applied, Df,  beta, delta, horizontal_force, vertical_force,
                                          cohesion, phi, gamma, FS, water_table, theta_direction)