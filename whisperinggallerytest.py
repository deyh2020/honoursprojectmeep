# Testing a whispering Gallery Mode

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# Characteristic lengthscale
a = 1  # micrometers

resolution = 20  # pixels/um

sx = 16  # size of cell in X direction
sy = 16  # size of cell in Y direction
cell = mp.Vector3(sx, sy, 0)

# Creating perfectly matched layers
dpml = 1.0
pml_layers = [mp.PML(dpml)]

# Size of resonator and waveguide
r = 6

# Whispering gallery resonator
geometry = [mp.Cylinder(material=mp.Medium(epsilon=12),
                        center=mp.Vector3(),
                        radius=r,
                        height=mp.inf)]

# Source Wavelength
wvln = 1.538  # micrometers

# Scattering source
sources = [mp.Source(mp.ContinuousSource(frequency=a / wvln),
                     component=mp.Ez,
                     center=mp.Vector3(-6, 0))]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

sim.use_output_directory()

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
        until=200)

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
