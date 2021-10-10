import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N
N = 1
# Number infected
I_0 = 1.27 * 10**(-6) 
# Number recovered 
R_0 = 0
# Initial number susceptible to infection
S_0 = (N - I_0 - R_0)

# b is contacts per day needed to spread disease
# b = 1/2
b = 2
# k is the fixed fraction of recovery per day 
k = 1/3

alpha = 0 # 0 to 1 for 0 to 100% intensity

# time points (in days)
t = np.linspace(0, 100, 140) 

# Initial conditions 
y0 = S_0, I_0, R_0

# The SIR model differential equations.
def model(y, t, N, b, k, alpha):
    S,I,R = y
    dS = -b * (S * I)/N
    dI = (b * S * I / N - k * I) * (1 - alpha)
    dR = k * I
    return dS, dI, dR

# Integrate the SIR equations over the time grid, t.
y = odeint(model, y0, t, args=(N, b, k, alpha))

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.plot(t,y)
plt.xlabel('time (days)')
plt.ylabel('Number of cases')
plt.legend(['S = Suseptible', 'I = Infected', 'R = Recovered'], loc='right')
plt.show()