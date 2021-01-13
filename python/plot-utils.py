import sys
import numpy as np
import itertools as IT
import matplotlib
from matplotlib import rc, font_manager
matplotlib.rcParams['text.usetex'] = True
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import pathlib

name = []
P1 = 0
P2 = 1

area = 50
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rc('grid', linestyle=":", color='grey')
plt.title(r'\textbf{Design/Parameter space}')
plt.grid(True)

snap = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=7, label='Snapshots')
targ = mlines.Line2D([], [], color='crimson', marker='X', linestyle='None', markersize=7, label='Targets')
plt.legend(handles=[snap, targ])


file = pathlib.Path('list-snapshots.csv')
if file.exists ():
	Sdataset = pd.read_csv('list-snapshots.csv')
	for col in Sdataset.columns:
		name.append(col)

	params_value = pd.DataFrame(Sdataset)
	Sparam1 = params_value.iloc[:,P1]
	Sparam2 = params_value.iloc[:,P2]

	# Plot snapshots
	plt.scatter(Sparam1, Sparam2, marker='o', s=area, facecolors='k', edgecolors='k')
    
else:
    print ('')
    print (' MESSAGE: File list-snapshots.csv does not exist ...')
    snE = True


file = pathlib.Path('list-targets.csv')
if file.exists ():
	Tdataset = pd.read_csv('list-targets.csv')
	for col in Tdataset.columns:
		name.append(col)

	params_value = pd.DataFrame(Tdataset)
	Tparam1 = params_value.iloc[:,P1]
	Tparam2 = params_value.iloc[:,P2]

	# Plot targets
	plt.scatter(Tparam1, Tparam2, marker='X', s=area, color='crimson')
    
else:
    print ('')
    print (' MESSAGE: File list-targets.csv does not exist ...')
    tgE = True


plt.xlabel(name[P1])
plt.ylabel(name[P2])

plt.show()
print('')