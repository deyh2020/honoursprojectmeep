import numpy as np
from scipy import special
import matplotlib.pyplot as plt

wvln = np.arange(1480,1551,10)
neff1 = np.array([0.10385,0.10329,0.36788,0.31480,0.24872,0.15356,0.13057,0.27017])
neff2 = np.array([0.28062,0.40801,0.36254,0.27022,0.24156,0.14183,0.12626,0.16779])

fig1, ax1 = plt.subplots()

ax1.plot(wvln, neff1, label='Excited Face')
ax1.plot(wvln, neff2, label='Non-excited Face')

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("Effective index")

ax1.set_title("Dispersion Diagram of Photonic Crystal Waveguide")
ax1.legend()
plt.show()