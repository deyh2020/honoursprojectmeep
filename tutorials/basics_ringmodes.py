import meep as mp
import numpy as np
import matplotlib.pyplot as plt

a = 1
wvln = 1.530
fcen = 1 / wvln  # pulse center frequency

complexPerm = 0.101
realPerm = 12
cond = 2 * np.pi * fcen * complexPerm/realPerm

w = 1  # width of waveguide
r = 10  # inner radius of ring
pad = 4  # padding between waveguide and edge of PML
dpml = 2  # thickness of PML
sxy = 2 * (r + w + pad + dpml)  # cell size

cell = mp.Vector3(sxy, sxy)

c1 = mp.Cylinder(radius=r + w, material=mp.Medium(epsilon=12, D_conductivity=cond))
c2 = mp.Cylinder(radius=r)

ring = [c1, c2]
cylinder = [c1]

df = 0.01  # pulse frequency width
src = mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r + 0.1))

sim = mp.Simulation(cell_size=cell,
                    geometry=cylinder,
                    sources=[src],
                    resolution=20,
                    boundary_layers=[mp.PML(dpml)])

sim.use_output_directory()

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(r + 0.1), fcen, df)),
        until_after_sources=300)

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
