# Animating 2D HDF5 data output of code into .gif format

import os
from tables import *

# Name of the python script, don't add .py on the end
filename = "whisperinggallerytest"  # input("Script Name: ")

# Uses pytables to determine the number of time slots
h5file = open_file(filename + "-out/" + filename + "-ez.h5", mode="r")
shape = h5file.root.ez.shape
h5file.close()

length = str(shape[2] - 1)
print(length)

# os.system("cd " + filename + "-out")
os.chdir(filename + "-out")

# Using h5l "to convert to png
# os.system("ls")
os.system(
    "h5topng -t 0:" + length + " -R -Zc dkbluered -a yarg -A " + filename + "-eps-000000.00.h5 " + filename + "-ez.h5")
# Using imagemagick to convert to gif
os.system("convert " + filename + "-ez.t*.png ez.gif")
