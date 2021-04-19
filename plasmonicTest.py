# Plasmonic Nanoantennae?

import meep as mp
import matplotlib.pyplot as plt
from meep.materials import Ag

w = 0.1
length = 0.4
gap = 0.05
pad = 0.5
dpml = 0.5

# um_scale = 0.1 * um_scale

cell = mp.Vector3(2*length + gap + 2*pad + 2*dpml, w + pad + dpml, 0)

geometry = [mp.Block(mp.Vector3(length, w, mp.inf),
                     center=mp.Vector3((length + gap)/2, 0),
                     material=Ag),
            mp.Block(mp.Vector3(length, w, mp.inf),
                     center=mp.Vector3(-(length + gap)/2, 0),
                     material=Ag)
            ]

sources = [mp.Source(mp.ContinuousSource(frequency=1/730),
                     component=mp.Ez,
                     center=mp.Vector3(-gap/2, 0))]

pml_layers = [mp.PML(dpml)]

resolution = 20

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution,
                    force_complex_fields=False)

sim.run(until=200)

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
