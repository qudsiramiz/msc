data_dir = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/'
save_dir = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/'

fnames = sort( glob(data_dir + '*.hf') )

fnames = list( fnames )

for f in fnames[-4:] :

	print(f)
	d = {}

	dat = hf.File( f )

	sc_n.append( f[36:-47] )

	for key in dat.keys() :

		d[key] = dat[key]

	df = pd.DataFrame( d )

	df.index = pd.to_datetime( df.datetime, unit='s' )


