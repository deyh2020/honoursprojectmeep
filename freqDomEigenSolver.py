# Testing a whispering Gallery Mode

import meep as mp
import numpy as np

# import matplotlib.pyplot as plt

# Characteristic lengthscale
a = 1  # micrometers

resolution = 10  # pixels/um

# Creating perfectly matched layers
dpml = 1.0
pml_layers = [mp.PML(dpml)]

# Source Wavelength
wvln = 1.538  # micrometers

fcen = 1 / wvln

n = 3.4483 + 1j * 1.0901e-13

# Setting real and complex permitivities, and material loss
# complexPerm = 0
# realPerm = 12
perm = np.power(n, 2)
complexPerm = np.imag(perm)
realPerm = np.real(perm)
cond = 2 * np.pi * fcen * complexPerm / realPerm

# Some integer
m = 10

# Size of resonator and waveguide
r = m * wvln / (2 * np.pi * np.sqrt(12))
sep = 0.1
w = 0.1
pad = 2

sx = 2 * (r + pad + dpml)  # size of cell in X direction
sy = sx  # size of cell in Y direction
cell = mp.Vector3(sx, sy, 0)

geometry = [mp.Cylinder(material=mp.Medium(epsilon=realPerm, D_conductivity=cond),
                        center=mp.Vector3(),
                        radius=r,
                        height=mp.inf),  # Whispering gallery resonator
            # mp.Block(mp.Vector3(w, mp.inf, mp.inf),
            #         material=mp.Medium(epsilon=12),
            #         center=mp.Vector3(-(r + sep + w / 2), 0, 0)
            #         )
            ]

# Waveguide Source
# sources = [mp.Source(mp.ContinuousSource(frequency =f cen),
#                     component = mp.Ez,
#                     center = mp.Vector3(-(r + sep + w / 2), -sy / 2 + pad),
#                     size = mp.Vector3(w, 0))
#           ]

# Scattering Source in resonator
sources = [mp.Source(mp.ContinuousSource(frequency=fcen),
                     component=mp.Ez,
                     center=mp.Vector3(-r + 0.1, 0))
           ]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=True)

# sim.use_output_directory()

tol = 1e-7
maxiters = 100
cwtol = tol * 1e-3
cwmaxiters = 10000
L = 20

sim.init_sim()
eigfreq = sim.solve_eigfreq(tol, maxiters, fcen, cwtol, cwmaxiters, L)
Q = eigfreq.real / (-2 * eigfreq.imag)

print(eigfreq)
print(Q)
