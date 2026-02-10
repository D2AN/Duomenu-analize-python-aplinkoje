import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
print("Pradiniai duomenys")
v0=float(input("Greitis: (m/s)"))
r0=float(input("Atstumas: (m)"))
# alfa=float(input("Skridimo kampas: (deg)"))
# sukimasis=float(input("Sukimasis: (deg)"))
m_met=float(input("Mase: (g)"))
r_met=float(input("Radijusas: (m)"))

inclination =math.radians(45)

Cd=0.47 #Drag koefecientas sferos

G=6.6743*(10**(-11))
R=6.371*10**6
M=5.9722*(10**24)
rho0=1.225 # sea level density
H=8500
omegaa=2*math.pi / 86164
dt = 1
theta = 0
x_lst = []
y_lst = []
z_lst = []

# phi=math.radians(sukimasis) 
# alfa = math.radians(alfa)

r=r0
vr=0  #cia sittu atveju graziai sonu skrenda
vt=v0
# L=r*vt #kampinis momentas
t=0
t_lim = 200000


A=math.pi*(r_met**2)

h_lst=[]
v_lst=[]
t_lst=[]

orbital_poz=0

while r>R and t < t_lim:
    alt=r-R
    g=(G*M) / (r**2)
    
    a_drag=0
    v_total = math.sqrt(vr**2 + vt**2)
    
    if alt < 120000:  # Kai kunas atmosferoj
        rho = rho0 * math.exp(-alt / H)
        # Drag formula: a = (0.5 * rho * v^2 * Cd * A) / mass
        a_drag= (0.5 * rho * v_total**2 * Cd * A) / (m_met/1000)
    
    ar = -g + (vt**2 / r) - (a_drag * (vr / v_total) if v_total > 0 else 0) #radianinis pagreitis
    at = -(a_drag * (vt / v_total) if v_total > 0 else 0) # tangentinis pagreitis
    vr += ar*dt
    vt += at*dt
    r += vr*dt
    
    d_angle = (vt/r) * dt
    orbital_poz+= d_angle
    # theta += d_angle * math.cos(alfa)
    # phi += d_angle * math.sin(alfa)
    
    v_total = math.sqrt(vr**2 + vt**2)
    
    x = r * math.cos(orbital_poz)
    y = r * math.sin(orbital_poz) * math.cos(inclination)
    z = r * math.sin(orbital_poz) * math.sin(inclination) 
    x_lst.append(x)
    y_lst.append(y)
    z_lst.append(z)


    # ve=omegaa*R*math.cos(fi)
    # vimpact=math.sqrt(vr**2+(vt-ve)**2)
    
    h_lst.append(r-R)
    v_lst.append(v_total)
    t_lst.append(t)
    
    t += dt

print("tasku skaicius:", len(t_lst))


print(t_lst)
print(h_lst)

plt.figure()
plt.plot(t_lst, h_lst)
plt.xlabel("Laikas (s)")
plt.ylabel("Aukštis (m)")
plt.title("Meteorito atstumas laiko atžvilgiu")
plt.grid()

plt.figure()
plt.plot(t_lst, v_lst)
plt.xlabel("Laikas (s)")
plt.ylabel("Greitis (m/s)")
plt.title("Meteorito greitis laiko atžvilgiu")
plt.grid()


fig = plt.figure(figsize=(8, 8)) # Use a square figure size
ax = fig.add_subplot(111, projection='3d')

# 1. Draw the Earth (stays the same)
phi_grid, theta_grid = np.mgrid[0:np.pi:20j, 0:2*np.pi:20j]
xe = R * np.sin(phi_grid) * np.cos(theta_grid)
ye = R * np.sin(phi_grid) * np.sin(theta_grid)
ze = R * np.cos(phi_grid)
ax.plot_wireframe(xe, ye, ze, color="blue", alpha=0.2)

# 2. Plot the Meteorite Path
ax.plot(x_lst, y_lst, z_lst, color='red', label='Meteorito kelias')

# --- THE PROPORTIONALITY FIX ---
# Combine all data to find the absolute limits
all_x = np.array(x_lst + [R, -R])
all_y = np.array(y_lst + [R, -R])
all_z = np.array(z_lst + [R, -R])

# Find the center of the path
mid_x = (all_x.max() + all_x.min()) * 0.5
mid_y = (all_y.max() + all_y.min()) * 0.5
mid_z = (all_z.max() + all_z.min()) * 0.5

# Find the largest distance traveled in any one direction
max_range = max(all_x.max()-all_x.min(), 
                all_y.max()-all_y.min(), 
                all_z.max()-all_z.min()) * 0.5

# Set all axes to use that SAME maximum distance
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Force the box to be a cube
ax.set_box_aspect((1, 1, 1)) 
# -------------------------------

plt.legend()
plt.show()

# TO DO, reikia padaryt jog atkreiptu emesi i trinti
# Padaryti jei sudega, per kiek laiko sudegs, jei nesudega, su kokia jega atsitrenks i zeme, ir ar isnaikins zmonyja
