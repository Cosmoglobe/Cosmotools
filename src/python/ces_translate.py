import quiet
import os

_default_l2filename="/data4/quiet/runlist_l2_new.txt"
_current_l2filename=_default_l2filename
_l2data_c2r=None
_l2data_r2c=None

_default_chicagofilename="/data4/quiet/calib_data/chicago/ces_defs.txt"
_current_chicagofilename=_default_chicagofilename
_chicagodata_c2r=None
_chicagodata_r2c=None

def ces_translate(pairs,
	reverse=False,
	l2filename=_default_l2filename,
	chicagofilename=_default_chicagofilename):
	"""
	Translate segment numbers between the Chicago system and Oslo system.
	This is modelled after the ces_translate cpp utility by Sigurd but is accessible from python more easily.
	It generates the same results except, bizarrely for CES (378,0), which is rejected by 
	
	Takes a LIST OF PAIRS, as the argument i.e.:
	pairs=[(2000,1),(2000,2),(2000,3)]
	ces_translate(pairs)
	
	This means it is slightly ugly to translate a single pair:
	ces_translate( [(2000,1)] )
	
	The return will be a list of the translated pairs.  If there is no translation the run number 
	returned will be -1, so do check for this.
	
	Default is to translate chicago->l2.
	You can also specify reverse=True to translate l2-> chicago.

	By default it uses the files
	l2filename = "/data4/quiet/runlist_l2_new.txt"
	chicagofilename = "/data4/quiet/calib_data/chicago/ces_defs.txt"
	for the translations, but you can specify these as keyword arguments.
	
	The data from these files is saved as module-level variables, so the function will
	be slightly slower the first time you run it.  If you specify new filenames it will reload the data.
	
	"""
	global _l2data_c2r,_l2data_r2c,_current_l2filename
	global _chicagodata_c2r,_chicagodata_r2c,_current_chicagofilename


	#load the level 2 data
	if (_l2data_c2r is None) or (_l2data_r2c is None) or (_current_l2filename!=l2filename):
		_current_l2filename=l2filename
		segments=quiet.L2Runlist(_current_l2filename).segments()
		_l2data_c2r={}
		_l2data_r2c=[]
		for ID,segment in segments.iteritems():
			_l2data_c2r[ID]=(segment.startDate,segment.endDate)
			_l2data_r2c.append((segment.startDate,segment.endDate,ID[0],ID[1]))

	#Load the chicago data
	if (_chicagodata_r2c is None ) or (_chicagodata_c2r is None) or (_current_chicagofilename!=chicagofilename):
		_current_chicagofilename=chicagofilename
		_chicagodata_c2r={}
		_chicagodata_r2c=[]
		for line in open(_current_chicagofilename):
			run,seg,start,end=line.split()
			run=int(run)
			seg=int(seg)
			start=float(start)
			end=float(end)
			_chicagodata_c2r[(run,seg)]=(start,end)
			_chicagodata_r2c.append((start,end,run,seg))


	if reverse:
		from_data = _l2data_c2r
		to_data = _chicagodata_r2c
	else:
		from_data =_chicagodata_c2r
		to_data = _l2data_r2c
	output=[]
	for pair in pairs:
		try:
			n=len(pair)
			assert n==2
		except (TypeError,ValueError):
			raise ValueError("ces_translate takes a list of (run,seg) pairs")
		pair=tuple(pair)
		if pair not in from_data:
			if reverse:
				output.append((-1,1))
			else:
				output.append((-1,0))
			continue
		mjd_range=from_data[pair]
		mjd_mid=(mjd_range[0]+mjd_range[1])/2.0
		for start,end,run,seg in to_data:
			if start<mjd_mid<end:
				output.append((run,seg))
				break
		else:  #funky little python for..else construct.  This gets called if the loop does not break
			if reverse:
				output.append((-1,0))
			else:
				output.append((-1,1))
	return output	

