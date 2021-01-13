import quiet

#set the name of the level 2 dir (so we can work out where the files are)
level2 = "/data4/quiet/level2/"

#load the runlist into a structure
runlist = quiet.L2Runlist("/data4/quiet/runlist_l2.txt")

#get the full table of segments from the runlist
segments = runlist.segments()

#loop through the segment IDs
for (run_id,segment_id) in segments:

	#from the table, get the segment info (including the filename)
	segment_info = segments[(run_id,segment_id)]

	#get the full path to the file with this segment data in
	filename = segment_info.filename(level2)

	# load the data from the file
	data = quiet.l2_read(filename)

	#loop through the modules in the data we just loaded.
	for module in data:
		
		#compute some number x from the module data (the TOD is stored in module.tod.Q1, etc. e.g.
		x = module.tod.Q1.std()  #compute the standard deviation of the module data.

		#output.  we could also just save the number and do the plotting in python.
		print run_id, segment_id, module.number[0], x
		
