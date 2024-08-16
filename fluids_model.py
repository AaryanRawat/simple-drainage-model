import numpy as np
from scipy.optimize import fsolve
import cmath

def reynolds_number(v, D, rho=984, mu=0.001):
    return (rho * v * D) / mu

def colebrook_white(Re, D, epsilon=0.0001):
    def equation(f):
        return 1/np.sqrt(f) + 2 * np.log10((epsilon / (3.7 * D)) + (2.51 / (Re * np.sqrt(f))))
    
    f_initial = 0.02  # Initial guess for the friction factor
    friction_factor, = fsolve(equation, f_initial)
    return friction_factor

def velocity_with_head_loss(h, D, L, g=9.81, rho=984, mu=0.001, epsilon=0.0001):
    # Initial velocity approximation using Torricelli's Law
    v_initial = np.sqrt(2 * g * h)
    
    # Calculate Reynolds number with initial velocity
    Re_initial = reynolds_number(v_initial, D, rho, mu)
    
    # Calculate friction factor using the initial Reynolds number
    f = colebrook_white(Re_initial, D, epsilon)
    
    # Calculate velocity considering head loss
    v = np.sqrt((2 * g * h) / (1 + f * (L / D)))
    
    return v

def flow_rate(v, D):
    A = np.pi * (D / 2) ** 2  # Cross-sectional area of the discharge tube
    Q = A * v  # Flow rate
    return Q

def drainage_time(A_tank, D, L, h_initial, h_final, dt=0.001):
    h = h_initial
    t1 = 0
    
    while h > h_final:
        v = velocity_with_head_loss(h, D, L)
        Q = flow_rate(v, D)
        dh = Q * dt / A_tank
        h -= dh
        t1 += dt
    
    return t1

def t_joint_flow_rate(v_inlet, D_inlet, D_outlet, L_outlet, K_t_joint, rho=984, mu=0.001, epsilon=0.0001, g=9.81):
    # Calculate flow rate in the inlet pipe
    Q_inlet = flow_rate(v_inlet, D_inlet)
    
    # Assuming the flow splits evenly between the two branches
    Q_outlet = Q_inlet/2
    
    # Cross-sectional area of the outlet pipes
    A_outlet = np.pi * (D_outlet / 2) ** 2
    
    # Calculate velocity in the T-joint branches using the continuity equation
    v_outlet = Q_outlet / A_outlet
    
    # Calculate Reynolds number for the outlet
    Re_outlet = reynolds_number(v_outlet, D_outlet, rho, mu)
    
    # Calculate friction factor for the outlet pipe
    f_outlet = colebrook_white(Re_outlet, D_outlet, epsilon)
    
    # Calculate head loss in each branch due to friction and the T-joint
    h_loss_friction = f_outlet * (L_outlet / D_outlet) * (v_outlet ** 2) / (2 * g)
    h_loss_t_joint = K_t_joint * (v_outlet ** 2) / (2 * g)
    
    # Total head loss
    h_loss_total = h_loss_friction + h_loss_t_joint

    # Effective head
    h_eff = h_initial - h_loss_total
    
    # Effective velocity considering total head loss via Toricelli's Law
    v_outlet_eff = cmath.sqrt(v_outlet * g * h_eff)

    # Effective Flow Rate considering total head loss
    Q_outlet_eff = v_outlet_eff.real * A_outlet
    
    return(Q_outlet_eff * 2)  # Return the velocity and the effective flow rate for both branches combined

def drainage_time_t_joint(A_tank, D_inlet, D_outlet, L_inlet, L_outlet, K_t_joint, h_initial, h_final, dt=0.01):
    h = h_initial
    t2 = 0
    
    while h > h_final:
        v_inlet = velocity_with_head_loss(h, D_inlet, L_inlet)
        Q_eff = t_joint_flow_rate(v_inlet, D_inlet, D_outlet, L_outlet, K_t_joint)
        dh = Q_eff * dt / A_tank
        h -= dh
        t2 += dt
    
    return t2

A_tank = 0.0832  # Cross-sectional area of the tank (m^2)
D = 0.00794  # Diameter of the discharge tube (m)
L = 0.2  # Length of the discharge tube (m)
L_outlet = 0.02 # T joint branch length
h_initial = 0.1  # Initial height of the water (m)
h_final = 0.02  # Final height of the water (m)
D_outlet = 0.01125  # Diameter of the T-joint outlet pipes (11.25 mm in meters) not used right now
K_t_joint = 0.5 # Assumed loss co-efficient

# Calculate time to drain the tank
t_drain = drainage_time(A_tank, D, L, h_initial, h_final)
t_drain_t_joint = drainage_time_t_joint(A_tank, D, D_outlet, L, L_outlet, K_t_joint, h_initial, h_final)

print(f"Time to drain the tank: {t_drain} seconds")
print(f"Time to drain the tank with a T-joint outlet: {t_drain_t_joint} seconds")
