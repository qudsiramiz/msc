from glob import glob
import numpy as np
import h5py as hf
import pandas as pd

fnames = np.sort(glob('/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v05/*.p'))

sc_marker = [ 'd', '^', 'v', 'o', 'P', 'X', 'p', '<', '>', '*', 's', 'D'  ]
sc_color = [ 'tab:brown', 'tab:orange', 'tab:olive', 'm', 'k', 'tab:gray',
             'tab:green', 'tab:purple', 'tab:pink', 'tab:red', 'tab:blue', 'tab:cyan']
#plt.figure()
plt.clf()

for i, f in enumerate(fnames[:]) :

	if ( fnames[0][-1] == 'f' ) :

		df = hf.File(f)

		print('file read for %s' %(f[51:-51]) )

	elif (fnames[0][-1] == 'p' ):

		df = pd.read_pickle(f)

		print('file read for %s' %(f[51:-50]) )

	plt.plot( df['sc_r'], df['heliographicLatitude'], marker=sc_marker[i], 
	          color=sc_color[i], ms=0.1, label=f[51:-50])

#df = pd.read_pickle(fnames[-1])

#print('file read for %s' %(fnames[11][51:-50]) )
#plt.plot( df['sc_r'], df['heliographicLatitude'], marker=sc_marker[11], 
#          color=sc_color[11], ms=0.1, label=f[51:-50])

plt.xlim( 0.1, 100 )
plt.xscale( 'log' )

plt.xlabel('Distance from Sun (AU)', fontsize=20 )
plt.ylabel('Heliographic Latitude', fontsize=20 )

plt.legend( fontsize=16, loc=2, ncol=3 )
plt.show()

'''
## For pickle files

for i, f in enumerate(fnames) :

	df = pd.read_pickle(f)

	print('file read for %s' %(f[51:-51]) )

	plt.plot( df['sc_r'], df['heliographicLatitude'], marker=sc_marker[i], color=sc_color[i], ms=0.1, label=f[51:-50] )

plt.xlim( 0.1, 100 )
plt.xscale( 'log' )

plt.xlabel('Distance from Sun (AU)', fontsize=20 )
plt.ylabel('Heliographic Latitude', fontsize=20 )

plt.legend( fontsize=16, loc=2, ncol=3 )
plt.show()
'''
