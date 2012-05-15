"""Compare accuracy of the moon_phase function with the ephem library"""

import ephem
import datetime
import moon_phase
import numpy as np

start_time = datetime.datetime(2011, 12, 1, hour=12)
delta_t = datetime.timedelta(hours=1)

norrk = ephem.Observer()
norrk.lat, norrk.lon = '58.5875', '16.1875'
norrk.pressure = 0
time_t = start_time

tlist = []
ephem_pha = []
pha = []
for idx in range(1000):
    norrk.date = time_t.strftime('%Y/%m/%d %H:%M')
    m = ephem.Moon(norrk)
    ephem_phase = m.phase
    phase = moon_phase.moon_phase(time_t)


    #diff = 100*phase - ephem_phase
    tlist.append(time_t)
    #print time_t.strftime("%Y-%m-%d %H:%M"), diff
    ephem_pha.append(ephem_phase)
    pha.append(phase)
    time_t = time_t + delta_t


pha = np.array(pha) * 100
ephem_pha = np.array(ephem_pha)

# Plot the moon phases to compare:
import matplotlib.pyplot as plt

fig = plt.figure()

ax = fig.add_subplot(111)
ax.scatter(ephem_pha, ephem_pha - pha)

ax.set_xlabel('Ephem moon phase (%)', fontsize=20)
ax.set_ylabel('Difference - Ephem-myprogam (%)', fontsize=20)
ax.set_title('Comparing Moon phase derivation')
ax.grid(True)

plt.savefig('./moonphase_compare.png')
