from glob import glob
import numpy as np
import h5py as hf

fnames = np.sort(glob('/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v01/*.hf'))

sc_marker = [ 'd', '^', 'v', 'o', 'P', 'X', 'p', '<', '>', '*', 's', 'D'  ]
sc_color = [ 'tab:brown', 'tab:orange', 'tab:olive', 'm', 'k', 'tab:gray',
             'tab:green', 'tab:purple', 'tab:pink', 'tab:red', 'tab:blue', 'tab:cyan']

plt.clf()

for i, f in enumerate(fnames) :

	dat = hf.File(f)

	print('file read for %s' %(f[51:-51]) )

	plt.plot( dat['sc_r'], dat['heliographicLatitude'], marker=sc_marker[i], color=sc_color[i], ms=0.1, label=f[51:-51] )

plt.xlim( 0.1, 100 )
plt.xscale( 'log' )

plt.xlabel('Distance from Sun (AU)', fontsize=20 )
plt.ylabel('Heliographic Latitude', fontsize=20 )

plt.legend( fontsize=16, loc=2, ncol=3 )
plt.show()

