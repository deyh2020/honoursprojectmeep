import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator
import re

# Get rid of the header and rename the columns of COMSOL's csv before entering it in
data = pd.read_csv('Sparrow_3_Mode.csv')
#print(data)


#print(data.wvln)
#print(data.neff)
#print(np.real(data.neff))

fig1, ax1 = plt.subplots()

# data.neff = np.real(data.neff)

#print(data.neff)

wvln = data.wvln.to_numpy()
neff = data.neff.to_numpy()
conf = data.confinement.to_numpy() # This was an integral that got rid of antisymmetric modes

# print(float(neff))
reg_neff = [re.sub("i","j", neff[i]) for i in range(len(neff))] 
#Replaces i with j because python can't $#%*ing read i as an imaginary number

complex_neff = [complex(reg_neff[i]) for i in range(len(neff))]
real_neff = np.real(complex_neff)

real_wvln = [float(wvln[i]) for i in range(len(wvln))]
real_conf = [abs(float(conf[i])) for i in range(len(conf))]
# print(real_neff)

#cond = real_conf > 7e-18*np.ones(len(conf))
cond = True

arr_wvln = np.array(real_wvln)
arr_neff = np.array(real_neff)

print(type(arr_neff))

ax1.plot(arr_wvln[cond], arr_neff[cond], 'b.')

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("Effective index")
ax1.yaxis.set_major_locator(MaxNLocator(5))

ax1.set_title("Dispersion Diagram of SPARROW")

plt.show()