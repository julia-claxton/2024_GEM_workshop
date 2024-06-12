import os
import sys
# Activate python virtual environment for this project to sandbox its libraries
os.system("source " + sys.path[0] + "/python_environment/bin/activate")
# Library includes
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import numpy as np

"""
convert_gei_to_lla.py
This script reads a datafile of GEI x/y/z coordinates from <project top level>/data/TEMP_locations_gei.npy and
converts them to lat/long, saving the results to <project top level>/data/locations_latlong.npy.
The GEI location file is deleted by this script.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see ../LICENSE.txt for full license)
"""

# Load EMIC locations
filepath = os.path.realpath(__file__)
projectdir = os.path.dirname(os.path.dirname(filepath))

locations = np.load(f"{projectdir}/data/TEMP_locations_gei.npy")
latitude = np.zeros(len(locations))
longitude = np.zeros(len(locations))

for i in range(len(locations)):
    print(f"\r\t{round((i/len(locations)) * 100, 1)}%", end = "")
    coordinates = coord.Coords(locations[i], 'GEI', 'car', use_irbem = False)
    coordinates.ticks = Ticktock(['2000-01-01T12:00:00'], 'UTC')

    newcoords = coordinates.convert('GDZ', 'sph') # Convert from GEI to ALL (altitude/lat/long)

    latitude[i] = newcoords.lati[0]
    longitude[i] = newcoords.long[0]

# Save converted coords for Julia script
print("")
np.savez(f"{projectdir}/data/locations_latlong.npz", lat = latitude, lon = longitude)
os.remove(f"{projectdir}/data/TEMP_locations_gei.npy")