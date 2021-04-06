# Testing a whispering Gallery Mode

import meep as mp
import numpy as np
from scipy import special
import matplotlib.pyplot as plt

# import matplotlib.pyplot as plt

# Characteristic lengthscale
a = 1  # micrometers

resolution = 20  # pixels/um

# Creating perfectly matched layers
dpml = 1.0
pml_layers = [mp.PML(dpml)]

# Vacuum Wavelength of source
wvln = 1.538  # micrometers

# Gaussian Pulse Source parameters
freq = 1 / wvln
fcen = freq
print("Frequency = " + str(fcen))
# cen = 0.6496971259430409
df = 0.002

# Setting real and complex permitivities, and material loss
complexPerm = 0
realPerm = 12
cond = 2 * np.pi * freq * complexPerm / realPerm

# Refractive index and ratio with index outside resonator
n = np.sqrt(realPerm)
q = n / 1

# Angular mode number
# Radial mode number
m = 30
nrad = 1  # Set to 1 to get fundamental whispering gallery modes

# Which zero of the airy function we need
azs = special.ai_zeros(nrad)[0]
az = azs[nrad-1]

# Size of resonator and waveguide, using schillers approx
r = 1 / (2 * np.pi * freq) * (m + 1 / 2 + az * ((m + 1) / 3) ** (1 / 3)
                              - q / np.sqrt(q ** 2 - 1) + 3 * az / (2 ** (2 / 3) * 10 * (m + 1 / 2) ** (1 / 3))
                              + q ** 3 * az / (3 * 2**(1/3) * (q**2-1)**(3/2) * (m+1/2) ** (3/2)))
print("Radius = " + str(r))
sep = 0.1
w = 0.1
pad = 2

sx = 2 * (r + pad + dpml)  # size of cell in X direction
sy = sx  # size of cell in Y direction
cell = mp.Vector3(sx, sy, 0)

geometry = [mp.Cylinder(material=mp.Medium(epsilon=12, D_conductivity=cond),
                        center=mp.Vector3(),
                        radius=r,
                        height=mp.inf),  # Whispering gallery resonator
            #mp.Block(mp.Vector3(w, mp.inf, mp.inf),
            #         material=mp.Medium(epsilon=12),
            #         center=mp.Vector3(-(r + sep + w / 2), 0, 0)
            #         )
            ]

# Waveguide Source
#sources = [mp.Source(mp.GaussianSource(frequency=fcen, width=df),
#                     component=mp.Ez,
#                     center=mp.Vector3(-(r + sep + w / 2), -sy / 2 + pad),
#                     size=mp.Vector3(w, 0))
#           ]

# Scattering Source in resonator
sources = [mp.Source(mp.GaussianSource(frequency=fcen, width=df),
                     component=mp.Ez,
                     center=mp.Vector3(-r + 0.1, 0))
           ]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=False)

sim.use_output_directory()
h = mp.Harminv(mp.Ez, mp.Vector3(- r + 0.1), fcen, df)

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(h),
        until_after_sources=10000)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()

QAtThisR = [m.Q for m in h.modes]

# for i in np.arange(0,len(h.modes)-1,1):
#    mode = h.modes[i]
#    print(mode.Q)
#    QAtThisR = np.append(QAtThisR,mode.Q)

print(QAtThisR)
print(max(QAtThisR))
