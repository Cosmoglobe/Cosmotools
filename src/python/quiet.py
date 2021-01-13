"""
$Id$
"""
import numpy as np
import struct
import os
import glob
import sys
import StringIO
import bisect
import pyfits
import copy
#import secateur
euler_mascheroni = 0.57721566490153286060
N_MODULE_MAX=19
N_MODULE_MAX_Q=19
N_MODULE_MAX_W=90
NDIODE_MODULE=4
ALL_MODULES=np.arange(19)
POL_MODULES=np.arange(17)
TEMP_MODULES=np.arange(17,19)
SAMPLING_FREQUENCY=50.0

patch_coords = {
"patch_2":( 185, -42),
"patch_2a": (181, -39),
"patch_4": (75, -40),
"patch_4a": (78, -39),
"patch_6": (10, -45),
"patch_6a": (12, -48),
"patch_7": (330, -30),
"patch_7a": (332, -36),
"patch_a": (148, 6),
"patch_b": (21, -48),
"patch_gc":(-92, -30)
}


class FortranFile(file):
	"""File with methods for dealing with fortran unformatted data files"""
	def __init__(self,fname, mode='r', buf=0,endian="=",field_dtype=np.int32):
		"""Open the file for writing, defaults to big endian."""
		file.__init__(self, fname, mode, buf)
		self.endian=endian
		self.field_dtype=field_dtype
		try:
			self.swap = {
			"=":      {True:False, False:False},  #native byte order; never swap
			"!":      {True:True,  False:True},   #non-native byte order; always swap
			">":      {True:True,  False:False},  #big endian; swap if system little-endian
			"<":      {True:False, False:True},   #little endian; swap if system big-endian
			}[endian] [np.little_endian]
		except KeyError:
			raise ValueError("Endianness must be in [=,<,>,!], not %s" % endian)
		

	def readString(self):
		"""Read in a string with error checking"""
		l = struct.unpack(self.endian+'i',self.read(4))[0]
		str = self.read(l)
		if  struct.unpack(self.endian+'i',self.read(4))[0] != l:
			raise IOError('Error reading string from data file')
		return str

	def writeString(self,s):
		"""Write a string
        :Parameters:
          - `s`: the string to write
        """
		self.write(struct.pack(self.endian+'i',len(s)))
		self.write(s)
		self.write(struct.pack(self.endian+'i',len(s)))

	def readArrays(self,dtypes):
		return [self.readArray(dtype) for dtype in dtypes]

	def writeArrays(self,arrays):
		for data in arrays:
			self.writeArray(data)

	def advance(self,n=1):
		for i in xrange(n):
			nb=np.fromfile(self,self.field_dtype,1)
			if not len(nb)==1:
				raise IOError("Fortran array format not correct")			
			if self.swap:
				nb.byteswap(True)
			self.seek(nb[0],1) #seek from the current position
			nb_check = np.fromfile(self,self.field_dtype,1)
			if not len(nb)==1:
				raise IOError("Fortran array format not correct")
			if self.swap:
				data.byteswap(True)
				nb_check.byteswap(True)
			if not np.array_equal(nb,nb_check):
				raise IOError("Fortran array format not correct")


	def readArray(self,dtype):
		nb=np.fromfile(self,self.field_dtype,1)
		if not len(nb)==1:
			raise IOError("Fortran array format not correct")			
		if self.swap:
			nb.byteswap(True)
		n = nb/np.nbytes[dtype]
		data=np.fromfile(self,dtype,n)
		if not len(data)==n:
			raise IOError("Fortran array format not correct")
		nb_check = np.fromfile(self,self.field_dtype,1)
		if not len(nb)==1:
			raise IOError("Fortran array format not correct")
		if self.swap:
			data.byteswap(True)
			nb_check.byteswap(True)
		if not np.array_equal(nb,nb_check):
			raise IOError("Fortran array format not correct")
		return data

	def writeArray(self,data):
		data=np.array(data)
		nb=np.array(data.nbytes,dtype=self.field_dtype)
		if self.swap:
			nb = nb.byteswap()
			data = data.byteswap()
		nb.tofile(self)
		data.tofile(self)
		nb.tofile(self)

	
#JZ Should used namedtuple here but not always available
class Pointing(object):
	"""
	This class contains the pointing from a level 2 file, for a single scan, for a single module.
	It is a simple structure with theta and phi locations and a psi sky rotation angle.
	"""
	def __init__(self,theta,phi,psi):
		self.theta=theta
		self.phi=phi
		self.psi=psi
		
class TOD(object):
	"""
	This class contains the TOD from a level 2 file, for a single scan, for a single module.
	It is a simple structure with four elements, Q1,U1,U2,Q2, corresponding to the diodes.
	"""
	def __init__(self,Q1,U1,U2,Q2):
		self.Q1=Q1
		self.U1=U1
		self.U2=U2
		self.Q2=Q2
	def __getitem__(self,n):
		return (self.Q1,self.U1,self.U2,self.Q2)[n]
	def __iter__(self):
		return iter( (self.Q1,self.U1,self.U2,self.Q2))
	def __setitem__(self,n,x):
		if n==0:
			self.Q1=x
		elif n==1:
			self.U1=x
		elif n==2:
			self.U2=x
		elif n==3:
			self.Q2=x
		else:
			raise IndexError("Indices into TOD objects run from 0..3 : (Q1,U1,U2,Q2)")


class L2Map(object):
	"""Maps of format stored in L2 files"""
	def __init__(self, nhits,mask2map,nside):
		super(L2Map, self).__init__()
		self.nhits = nhits
		self.mask2map = mask2map
		self.nside=nside

class Module(object):
	"""
	module_number,nsamp,reduction_factor,num_reduced,pointing,tod,time,tp,tp_var
	This class contains the data from a single module in a level 2 file.
	It is a simple structure with a module number (number), number of samples (nsamp), pointing information (pointing), tim-ordered data (tod), and MJD (time)
	"""
	def __init__(self,number,nsamp,reduction_factor,num_reduced,pointing,tod,time,tp,tp_var,map):
		self.number=number
		self.reduction_factor=reduction_factor
		self.num_reduced=num_reduced
		self.nsamp=nsamp
		self.pointing=pointing
		self.tod=tod
		self.time=time
		self.tp=tp
		self.tp_var=tp_var
		self.map=map

	@classmethod
	def empty_module(cls,n,nreduced=0):
		nr=nreduced
		if nr==0: nr=n/100
		reduction_factor=n/nr
		tod=TOD(np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n))
		tp=TOD(np.zeros(nr),np.zeros(nr),np.zeros(nr),np.zeros(nr))
		tp_var=TOD(np.zeros(nr),np.zeros(nr),np.zeros(nr),np.zeros(nr))
		time=np.zeros(n)
		pointing=Pointing(np.zeros(n),np.zeros(n),np.zeros(n))
		return cls(0,n,reduction_factor,nr,pointing,tod,time,tp,tp_var)



class Scan(object):
	def __init__(self,id,targetName,targetType,nFile,status,startDate,endDate):
		self.targetName=targetName
		self.id=id
		self.targetType=targetType
		self.files=[]
		self.fileNums=[]
		self.fileStartDates=[]
		self.fileEndDates=[]
		self.startDate=startDate
		self.endDate=endDate
		self.nFile=nFile
		self.newFormat="Unknown"
	def rectify(self):
		self.nFile = len(self.files)
	def addFile(self,line):
		words=line.split()
		if len(words)==1:
			self.files.append(words[0])
			self.newFormat=False
		else:
			assert len(words)==4, "Either wrong format for l1 runlist or Joe misunderstood format."
			self.newFormat=True
			self.fileNums.append(int(words[0]))
			self.fileStartDates.append(float(words[1]))
			self.fileEndDates.append(float(words[2]))
			self.files.append(words[3])
	def finish(self):
		if not self.newFormat:  #not (True or unknown)  might be unknown if scan empty
			del self.fileNums,self.fileStartDates,self.fileEndDates
		assert len(self.files)==self.nFile

def numerify(iterable):
	"""
	A utility function to convert strings (or similar) from an iterable to the simplest 
	number type it can be - it tries integer first, then float, then gives up and yields the item 
	unchanged.
	"""
	for it in iterable:
		try:
			yield int(it)
		except:
			try:
				yield float(it)
			except:
				yield it

class L2Segment(object):
	def __init__(self, runID,targetName,targetType,id,startDate,endDate,az,el,dk,time):
		"""
		Info on a single segments of a CES, as read from a level 2 runlist.
		Attributes are:
		runID  -  scan ID, eg 431
		targetName - target name, eg patch_2a
		targetType - target type, eg cmb
		id - segment ID, e.g. 2
		startDate - starting MJD
		endDate - ending MJD
		az - scan azimuth
		el - scan elevation
		dk - deck angle
		time - time of day, I think.
		"""
		super(L2Segment, self).__init__()
		self.runID = runID
		self.targetName = targetName
		self.targetType = targetType
		self.id=id
		self.startDate = startDate
		self.endDate = endDate
		self.az = az
		self.el = el
		self.dk = dk
		self.time = time
	def filename(self,l2dir=""):
		"""
		Generate the file name of the level 2 file in the standard directory structure.
		Optional argument gives the level2 directory to use, otherwise assumed to be
		relative to current dir.
		e.g. /data4/quiet/level2/cmb/patch_2a/cmb_patch_2a_431_seg002.unf
		"""
		s=os.path.sep
		return l2dir+ "%s/%s/%s_%s_%d_seg%.3d.unf" % (
					   self.targetType,
					      self.targetName,
					         self.targetType,
					            self.targetName,
					               self.runID,
					                     self.id
										     )
				
		
class L2Scan(object):
	"""
	A run, consisting of several segments and associated information (usually read from a runlist)
	targetName - target name, eg patch_2a
	targetType - target type, eg cmb
	id - run ID, e.g. 431
	startDate - starting MJD
	endDate - ending MJD
	nmodule - number of modules in this scan
	nSegment - the number of segment scans in this run
	
	"""
	def __init__(self,id,targetName,targetType,nSegment,nmodule,startDate,endDate):
		self.targetName=targetName
		self.id=id
		self.targetType=targetType
		self.segments=[]
		self.startDate=startDate
		self.endDate=endDate
		self.nSegment=nSegment
	def addSegment(self,id,startDate,endDate,az,el,dk,time):
		"""
		Function called by the code that reads the runlist from a file, that adds segments on.
		If you call this manually to adjust the runlist, be sure to call "rectify" afterwards
		"""
		seg = L2Segment(self.id,self.targetName,self.targetType,id,startDate,endDate,az,el,dk,time)
		self.segments.append(seg)
	def finish(self):
		"""
		Having read in this scan from a runlist, make sure that the full number of segments the 
		file said would be present have been read in. 
		"""
		assert len(self.segments)==self.nSegment
	def rectify(self):
		"""
		Set the number of segments attribute to equal the number of contained segments.
		Use this after you modify the runlist.
		"""
		self.nSegment = len(self.segments)


class Target(object):
	"""
	Target information from a runlist (L1 or L2).
	name - eg patch_2a
	targetType - eg cmb
	nScan - number of scans that hit the patch
	scans - a dictionary of the runs that hit this target.
			keys are run IDs, values are scan objects.
	
	"""
	def __init__(self,name,targetType,nScan):
		self.name=name
		self.nScan=nScan
		self.scans={}
		self.targetType=targetType
	def addScan(self,scan):
		"""
		Called by the runlist read-in code to add a scan
		"""
		self.scans[scan.id]=scan
	def finish(self):
		"""
		Having read in this target from a runlist, make sure that the full number of runs the 
		file said would be present have been read in. 
		"""
		assert len(self.scans)==self.nScan
	def rectify(self):
		"""
		Set the number of scans attribute to equal the number of contained scans.
		Use this after you modify the runlist.
		"""		
		self.nScan = len(self.scans)
	def segments(self):
		"""
		Return a dictionary of all the segments that hit this scan.
		This is derived from the scan objects in self.scans.
		
		You cannot modify this directly - it will not affect the runlist.
		Instead you need to access the scans in self.scans directly.
		"""
		segments={}
		for scanID,scan in self.scans.items():
			for segment in scan.segments:
				segmentID=segment.id
				segments[(scanID,segmentID)]=segment
		return segments
	

class TargetType(object):
	"""
	A target type from a runlist, eg 'cmb' or 'calib'
	name - the name of the type
	nTarget - number of targets in the type
	targets - a dictionary of targets.
        keys are string, like 'patch_2a', vales are of class 'Target'
	"""
	def __init__(self,name,ntarget):
		self.name=name
		self.nTarget=ntarget
		self.targets={}
	def addTarget(self,target):
		self.targets[target.name]=target
	def removeTarget(self,targetName):
		del self.targets[targetName]
		self.rectify()
	def finish(self):
		assert len(self.targets)==self.nTarget
	def rectify(self):
		self.nTarget = len(self.targets)


		
class InOutFiles(object):
	"""
	Helpful super-class (should probably be a mix-in but I don't understand those) for classes
	that read and write from files that are stdin/stdout by default.
	Use the properties "infile" and "outfile", which can be either strings (to open a file of that 
	name) or objects with ".read" or ".write" methods.  The latter can be useful to, for example, 
	read from URLs or write to StringIO buffers. 
	"""
	def getInfile(self):
		if hasattr(self,"infile_"):
			return self.infile_
		else:
			return sys.stdin
	def setInfile(self,infile):
		if hasattr(infile,"read"):
			self.infile_=infile
		else:
			self.infile_=open(infile,"r")
	def delInfile(self):
		if hasattr(self,"infile_"):
			del self.infile_
			
	infile=property(getInfile,setInfile,delInfile,"Input file property")

	def getOutfile(self):
		if hasattr(self,"outfile_"):
			return self.outfile_
		else:
			return sys.stdout
	def setOutfile(self,outfile):
		if hasattr(outfile,"write"):
			self.outfile_=outfile
		else:
			self.outfile_=open(outfile,"w")
	def delOutfile(self):
		if hasattr(self,"outfile_"):
			del self.outfile_
	outfile=property(getOutfile,setOutfile,delOutfile,"Output file property")
	
		
class L1Runlist(InOutFiles):
	"""
	A class correspnding to a level 1 runlist file.  It stores the contents of the file as a hierarchial collection
	of objects: TargetType -> Target -> Scan.
	
	These objects are stored in dictionaries.  This means that iterating through them with, for example,
	for ttype in runlist.targetTypes:
		print ttype
		#do work
	will yield the keys (the names) rather than the objects.
	
	To get both, use 
	for name,ttype in runlist.targetTypes.items():
		print name,ttype
		#do work
	
	Initialize with the name of the input file or a file-like object.
	"""
	def __init__(self,infile):
		if hasattr(infile,'name'):
			self.title=infile.name
		else:
			self.title=str(infile)
		self.infile=infile
		line=self.infile.readline()
		self.ntype=int(line)
		self.targetTypes={}
		for i in range(self.ntype):
			type=self.readType()
			self.targetTypes[type.name]=type
		del self.infile
	
	def removeTargetType(self,item):
		"""
		Remve a target type completely from the runlist, either by name or by object.
		e.g.
		runlist.removeTargetType("calib")
		"""
		if isinstance(item,TargetType):
			del self.targetTypes[item.name]
		else:
			del self.targetTypes[item]
		self.ntype=len(self.targetTypes)

	def removeTarget(self,target,targetType=None):
		"""
		Remove a given patch (target) from the runlist.
		Again, specify by name or by object.
		
		If you pass in the target name as a string then you also need to specify
		targetType; otherwise this is optional.
		"""
		if targetType is None:
			targetType=self.targetTypes[target.targetType]
		elif isinstance(targetType,TargetType):
			pass
		else:
			#Assume we have been passed the name of the target type
			targetType=self.targetTypes[targetType]
		if isinstance(target,Target):
			target=target.name
		targetType.removeTarget(target)
		
	def removeScan(self,scan):
		"""
		Remove a scan by object from the runlist.
		i.e., you do not pass an ID number to this function, you pass the scan object.
		(for now - I may upgrade this).
		"""
		targetType=self.targetTypes[scan.targetType]
		target=targetType.targets[scan.targetName]
		del target.scans[scan.id]
		target.rectify()
	
	def targets(self):
		"""
		Return a dictionary of the name:target_object pairs for all the targets
		in all the types in the runlist. 
		
		Editing the returned dict will not change the runlist, BUT changing the things in it will.
		"""
		targets = {}
		for targetType in self.targetTypes.values():
			for targetName,target in targetType.targets.items():
				targets[targetName] = target
		return targets
	
	def scans(self):
		"""
		Return a dictionary of all the run_number:scan_object pairs for all the targets
		for all the types.
		
		Editing the returned dict will not change the runlist, BUT changing the things in it will.
		"""
		scans={}
		for target in self.targets().values():
			for id,scan in target.scans.items():
				scans[id] = scan
		return scans
		
			
	def readType(self):
		"""
		Internal function used to read a type from the runlist file.
		You probably do not want to use this manually.
		"""
		line=self.infile.readline()
		typeName,ntarg=line.split()
		ntarg=int(ntarg)
		targetType=TargetType(typeName,ntarg)
		for i in range(ntarg):
			target=self.readTarget(typeName)
			targetType.addTarget(target)
		targetType.finish()
		return targetType

	def readTarget(self,typeName):
		"""
		Internal function used to read a target from the runlist file.
		You probably do not want to use this manually.
		"""
		line=self.infile.readline()
		words=line.split()
		targetName=words[0]
		nscan=int(words[1])
		target=Target(targetName,typeName,nscan)
#		self.targets[targetName]=target
		for i in xrange(nscan):
			scan = self.readScan(targetName,typeName)
			target.addScan(scan)
		target.finish()
		return target

	def toFile(self,filename):
		"""
		Output this runlist to a file.
		You can pass in either a filename (which will open and write to the file),
		or a file-like object (like sys.stdout, which will print to stdout), or anything
		else with a "write" method.  In particular, you might want to use the
		StringIO object like this:
		import StringIO
		buffer = StringIO.StringIO
		runlist.toFile(buffer)
		#The buffer now contains the printed runlist, but you can add more 
		#to the end of it.  To print it:
		buffer.seek(0)
		print buffer.read()
		
		"""
		self.outfile=filename
		self.outfile.write("%d\n"%self.ntype)
		targetTypes = self.targetTypes.values()
		targetTypes.sort(key=lambda targetType: targetType.name) #sort the types by name
		for targetType in targetTypes:
			self.writeType(targetType)
		del self.outfile

	def writeType(self,targetType):
		"""
		Internal function used to write a type to a runlist file.
		You probably do not want to use this manually.		
		"""
		self.outfile.write("%s %d\n" % (targetType.name,targetType.nTarget))
		targets = targetType.targets.values()
		targets.sort(key=lambda target:target.name)  #Sort the targets by name
		for target in targets:
			self.writeTarget(target)

	def writeTarget(self,target):
		"""
		Internal function used to write a target to a runlist file.
		You probably do not want to use this manually.		
		"""
		self.outfile.write("  %s %d\n" % (target.name,target.nScan))
		scans = target.scans.values()
		key = lambda scan: scan.id
		scans.sort(key=key)
		for scan in scans:
			self.writeScan(scan)

	def readScan(self,targetName,typeName):
		"""
		Internal function used to read a scan from a runlist file.
		You probably do not want to use this manually.		
		"""
		line=self.infile.readline()
		words=line.split()
		scanId=int(words[0])
		startDate=float(words[1])
		endDate=float(words[2])
		nFile=int(words[3])
		status=words[4]
		scan=Scan(scanId,targetName,typeName,nFile,status,startDate,endDate)
#		self.scans[scanId]=scan
		for i in xrange(nFile):
			filename=self.infile.readline().strip()
			scan.addFile(filename)
		scan.finish()
		return scan
	
	def writeScan(self,scan):
		"""
		I have not implemented writing out level 1 runlists, only reading them.
		If you need this, let me know.
		"""
		raise NotImplementedError("Writing L1 runlists not implemented")
	
	def __str__(self):
		"""
		The full text form of the runlist.
		In interactive mode you will this if you type:
		print runlist
		
		but not just:
		runlist
		"""
		output=StringIO.StringIO()
		self.toFile(output)
		output.seek(0)
		return output.read()
		

class L2Runlist(L1Runlist):
	"""
	A class correspnding to a level 2 runlist file.  Since there is a lost of similiar functionality,
	this class inherits from L1Runlist.
	
	It stores the contents of the file as a hierarchial collection
	of objects: TargetType -> Target -> Scan -> Segment.
	These objects are stored in dictionaries.  This means that iterating through them with, for example,
	for ttype in runlist.targetTypes:
		print ttype
		#do work
	will yield the keys (the names) rather than the objects.
	
	To get both, use 
	for name,ttype in runlist.targetTypes.items():
		print name,ttype
		#do work
	
	Initialize with the name of the input file or a file-like object.
	"""
	def __init__(self,infile):
		if isinstance(infile,str) or isinstance(infile,unicode):
			if os.path.isdir(infile):
				initDir=True
			else:
				initDir=False
		else:
			initDir=False
		if initDir:
			self.initDir(infile)
		else:	
			L1Runlist.__init__(self,infile)

	def initDir(self,infile):
		"""
		If needed I could code something up to generate a runlist from a directory structure.
		Right now this is not implemented.
		"""
		raise NotImplementedError

	def segments(self):
		"""
		Return a dictionary of the segments in all the scans in all the targets in all the types.
		
		Editing the returned dict will not change the runlist, BUT changing the things in it will.
		"""
		scans=self.scans()
		segments={}
		for scan in scans.values():
			for segment in scan.segments:
				segments[(segment.runID,segment.id)]=segment
		return segments


			
	def readScan(self,targetName,typeName):
		"""
		Internal function to read a scan from the runlist file.
		
		You probably do not want to use this manually.
		"""
		line=self.infile.readline()
		words=line.split()
		scanId=int(words[0])
		startDate=float(words[1])
		endDate=float(words[2])
		nSeg=int(words[3])
		nModule=int(words[4])
		scan=L2Scan(scanId,targetName,typeName,nSeg,nModule,startDate,endDate)
#		self.scans[scanId]=scan
		for i in xrange(nSeg):
			data=self.infile.readline().split()
			scan.addSegment(*numerify(data))
		scan.finish()
		return scan
		
	def writeScan(self,scan):
		"""
		Internal function to write a scan to a runlist file.
		
		You probably do not want to use this manually.
		"""
		line = "    %d %f %f %d %d\n"  % (scan.id,scan.startDate,scan.endDate,scan.nSegment,0)
		self.outfile.write(line)
		for seg in scan.segments:
			line = "      %d %f %f %f %f %f %f\n" % (seg.id,seg.startDate,seg.endDate,seg.az,seg.el,seg.dk,seg.time)
			self.outfile.write(line)

	def deleteGappyScans(self,l2dir):
		"""
		Delete scan segments after a missing segments.
		For example, if a scan has segments with IDs 1,2,3,5,6
		then segments 5 and 6 are removed.
		Scans missing segment 1 are removed completely
		"""
		self.deleteMissingSegments(l2dir)
		badScans=[]
		for scan in self.scans().values():
			segNumbers = [seg.id for seg in scan.segments]
			badSegStart=None
			if segNumbers[0]!=1:
				badSegStart=0
			else:
				for i in xrange(1,len(segNumbers)):
					if segNumbers[i] - segNumbers[i-1]!=1:
						badSegStart=i
						break
			if badSegStart is not None:
				badSegments=scan.segments[badSegStart:]
			else:
				badSegments=[]
			for seg in badSegments:
				scan.segments.remove(seg)
			scan.rectify()
			if scan.nSegment==0:
				badScans.append(scan)
		self.deleteScans(badScans)

	def filenames(self,l2dir):
		"""
		Return a list of filenames for all the segments in the scan, for
		the specified level 2 directory.
		"""
		output = []
		for scan in self.scans().values():
			for segment in scan.segments:
				output.append(l2dir + os.path.sep + segment.filename() )
		return output

	def deleteMissingScans(self,l2dir):
		"""
		Delete any scans that have any missing segment files from the runlist
		"""
		badScans = []
		for scan in self.scans().values():
			for segment in scan.segments:
				filename = l2dir + os.path.sep + segment.filename()
				if not os.path.exists(filename):
					badScans.append(scan)
					break
		self.deleteScans(badScans)
	def deleteMissingSegments(self,l2dir):
		"""
		Delete any segments from the runlist that do not have a corresponding file in the 
		segments file.
		"""
		badScans=[]
		for scan in self.scans().values():
			badSegments=[]
			for segment in scan.segments:
				filename = l2dir + os.path.sep + segment.filename()
				if not os.path.exists(filename):
					badSegments.append(segment)
			for segment in badSegments:
				scan.segments.remove(segment)
			scan.rectify()
			if scan.nSegment==0:
				badScans.append(scan)
		self.deleteScans(badScans)
			
	def deleteScans(self,scanList):
		"""
		Delete a list of scan objects from the runlist. (not by scan ID)
		"""
		targets=self.targets()
		for scan in scanList:
			target = targets[scan.targetName]
			del target.scans[scan.id]
			target.rectify()
		self.removeEmptyTargets()
		
	def removeEmptyTargets(self):
		"""
		Remove any targets from the runlist that have had all their scans removed.
		You should call this after manually fiddling with the scans.
		It is automatically called by the deleteScans method.
		"""
		emptyTargets=[]
		for targetName,target in self.targets().items():
			if not target.scans:
				emptyTargets.append((targetName,target.targetType))
		emptyTypes=set()
		for targetName,targetTypeName in emptyTargets:
			targetType = self.targetTypes[targetTypeName]
			del targetType.targets[targetName]
			targetType.rectify()
			if targetType.nTarget==0:
				emptyTypes.add(targetType.name)
		for typeName in emptyTypes:
			del self.targetTypes[typeName]
		self.ntype=len(self.targetTypes)
				

def acceptedListFromRunList(runlist):
	"""
	A function that builds and returns an AcceptedScanList
	from all the scans in the runlist.
	Nothing is excluded from any of the scans in the runlist.
	"""
	def make_scanlist():
		for scan in runlist.scans().values():
			runID=scan.id
			for seg in scan.segments:
				segID=seg.id
				yield (runID,segID)
	scanlist=make_scanlist()
	return AcceptedScanList(scanlist=scanlist)
			

class AcceptedScanList(object):
	"""
	A class representing an accept list - a list of scans that
	are accepted, and within 
	each, a list of diodes that are excluded.

	You can either load one of these from a file, build one
	from scratch with a list of scans,
	or use the utility function acceptedListFromRunList.
	"""
	def __init__(self,infile=None,scanlist=None):
		skip=False
		self.title=''
		if hasattr(infile,"read"):
			infile=infile
		elif infile is None:
			skip=True
		else:
			self.title=os.path.abspath(infile)
			infile=open(infile)
		self.scans={}
		if skip:
			for scan,seg in scanlist:
				shape = (N_MODULE_MAX,NDIODE_MODULE)
				included=np.ones(shape,dtype=bool)
				self.scans[(scan,seg)]=included
		else:
			for line in infile:
				if not line.strip():
					continue
				if line.startswith("#"):
					continue
				words=line.split()
				scan=int(words[0])
				seg=int(words[1])
				nex=int(words[2])
				shape = (N_MODULE_MAX,NDIODE_MODULE)
				included=np.ones(shape,dtype=bool)
				for i in xrange(nex):
					modnum=words[3+2*i]
					diode=words[4+2*i]
					included[modnum,diode]=False
				self.scans[(scan,seg)]=included

	def __eq__(self,other):
		if self.title != other.title: return False
		if len(self.scans)!= len(other.scans): return False
		for key in self.scans.iterkeys():
			if key not in other.scans: return False
		for key in self.scans.iterkeys():
			if not (self.scans[key]==other.scans[key]).all(): return False
		return True
	def explain_inequality(self,other):
		if self.title != other.title: return "Titles of accept lists are different: '%s'  vs  '%s' " % (self.title,other.title)
		if len(self.scans)!= len(other.scans): return "Lengths of accept lists are different: %d  vs  %d" %(len(self.scans),len(other.scans))
		for key in self.scans.iterkeys():
			if key not in other.scans: return "Scan (run,seg)=(%d,%d) is missing in second accept list"%key
		for key in self.scans.iterkeys():
			if not (self.scans[key]==other.scans[key]).all(): return "Scan (run,seg)=(%d,%d) has different exclusions in lists"%key
		return ""

	def __str__(self):
		output=StringIO.StringIO()
		self.toFile(output)
		output.seek(0)
		return output.read()


	def toFile(self,filename,header=None):
		if hasattr(filename,'read'):
			f=filename
		else:
			f=open(filename,"w")
		if header is not None: self.writeHeader(f,header)
		scans = self.scans.keys()
		scans.sort()
		for scan,seg in scans:
			included=self.scans[(scan,seg)]
			if not np.any(included): continue
			nex = included.size - included.sum()
			f.write("%d %d %d " % (scan,seg,nex))
			for modnum in range(N_MODULE_MAX):
				for diode in range(NDIODE_MODULE):
					if not included[modnum,diode]:
						f.write("%d %d " % (modnum,diode))
			f.write("\n")

        def writeHeader(self,f,header):
		f.write(header)

	def isIncluded(self,scan,seg,module,diode):
		return self.scans[(scan,seg)][module,diode]
	def isExcluded(self,*args):
		return not self.isIncluded(*args)
	def include(self,*args):
		"""
		Include some segment,module or diode from the list.
		number of args			meaning
			2 scan,segment: include all modules and diodes
			3 scan,segment,module: include all diodes
			4 scan,segment,module,diode
		"""
		if len(args)==2:
			scan,seg=args
			self.scans[(scan,seg)][:]=True
		elif len(args)==3:
			scan,seg,module=args
			self.scans[(scan,seg)][module,:]=True
		elif len(args)==4:
			scan,seg,module,diode=args
			self.scans[(scan,seg)][module,diode]=True
		else:
			 raise ValueError("AcceptedScanList.include should be called with 2-4 args - see docstring")
	def exclude(self,*args):
		"""
		Exclude some segment,module or diode from the list.
		number of args			meaning
			2 scan,segment: exclude all modules and diodes
			3 scan,segment,module: exclude all diodes
			4 scan,segment,module,diode
		"""
		if len(args)==2:
			scan,seg=args
			self.scans[(scan,seg)][:]=False
		elif len(args)==3:
			scan,seg,module=args
			self.scans[(scan,seg)][module,:]=False
		elif len(args)==4:
			scan,seg,module,diode=args
			self.scans[(scan,seg)][module,diode]=False
		else:
			 raise ValueError("AcceptedScanList.exclude should be called with 2-4 args - see docstring")
	def excludeSegmentDiode(self,scan,seg,module,diode):
		self.scans[(scan,seg)][module,diode] = False
	def includeSegmentDiode(self,scan,seg,module,diode):
		self.scans[(scan,seg)][module,diode] = True
	def excludeSegmentModules(self,scan,seg,module):
		for diode in range(NDIODE_MODULE):
			self.scans[(scan,seg)][module,diode] = False
	def includeSegmentModules(self,scan,seg,module):
		for diode in range(NDIODE_MODULE):
			self.scans[(scan,seg)][module,diode] = True
	def excludeModule(self,module):
		for ces in self.scans.itervalues():
			ces[module,:]=False
	def excludeDiode(self,module,diode):
		for ces in self.scans.itervalues():
			ces[module,diode]=False
	def includeModule(self,module):
		for ces in self.scans.itervalues():
			ces[module,:]=True
	def includeDiode(self,module,diode):
		for ces in self.scans.itervalues():
			ces[module,diode]=True
	def removeScan(self,scan,seg):
		del self.scans[scan,seg]
	def removeScansIn(self,scanList):
		bad=set()
		for scan in self.scans:
			if scan in scanList:
				bad.add(scan)
		for scan in bad:
			del self.scans[scan]
	def removeScansNotIn(self,scanList):
		bad=set()
		for scan in self.scans:
			if scan not in scanList:
				bad.add(scan)
		for scan in bad:
			del self.scans[scan]
	def removeEmptyScans(self):
		IDs = self.scans.keys()
		for ID in IDs:
			included=self.scans[ID]
			if not np.any(included):
				del self.scans[ID]

	def countIncludedDiodes(self,scan,seg,modules=ALL_MODULES):
		"For a CES, count number of included diodes"
		if (scan,seg) in self.scans.keys():
			return self.scans[(scan,seg)][modules].sum()
		else:
			return None
	def countExcludedDiodes(self,scan,seg,modules=ALL_MODULES):
		"For a CES, count number of excluded diodes"
		if (scan,seg) in self.scans.keys():
			inc=self.scans[(scan,seg)][modules]
			return inc.size-inc.sum()
		else:
			return None
	def calculateDiodeStats(self):
		"Count diodes in/out of each CES and return as dictionaries"
		indict={}
		outdict={}
		self.removeEmptyScans()
		for scan,seg in self.scans.keys():
			#if not np.any((scan,seg)): continue # Needed?
			ninc=self.countIncludedDiodes(scan,seg)
			nexcl=self.countExcludedDiodes(scan,seg)
			if ninc not in indict.keys():
				indict[ninc]=[]
			indict[ninc].append((scan,seg))
			if nexcl not in outdict.keys():
				outdict[nexcl]=[]
			outdict[nexcl].append((scan,seg))

		return indict,outdict

	def printSummary(self,verbose=False,silent=False):
		self.removeEmptyScans()
		included = self.countIncludedScanDiodes()
		excluded = self.countExcludedScanDiodes()
		if not silent:
			if verbose:
				fractionin=float(included)/(float(included+excluded))
				print ' acceptlist has %i CESs with diodes %i in %i out (%2.1f per cent)' \
				      % (self.acount(),self.countIncludedScanDiodes(),\
					 self.countExcludedScanDiodes(),100.0*fractionin)
			else:
				print ' acceptlist has %i CESs with diodes %i in %i out' \
				      % (self.acount(),self.countIncludedScanDiodes(),\
					 self.countExcludedScanDiodes())
		return included

	def countSummary(self):
		"Summarize numbers of included/excluded diodes" # (ugly)
		self.removeEmptyScans()
		incl,excl=self.calculateDiodeStats()
		for key in incl.keys():
			if key==0: continue # Ignore any with none included
			inx=incl[key]
			for key2 in excl.keys():
				if excl[key2]==inx:
				# assume incl+excl is ok
					print ' %i CESs have %i in %i out' % \
					      (len(incl[key]),key,key2)
		return

	def countIncludedScanDiodes(self):
		"count number of included CES-diodes"
		ncesdiode=0
		for scan,seg in self.scans.keys():
			ncesdiode+=self.countIncludedDiodes(scan,seg)
		return ncesdiode

	def countExcludedScanDiodes(self):
		"count number of excluded CES-diodes"
		ncesdiode=0
		for scan,seg in self.scans.keys():
			ncesdiode+=self.countExcludedDiodes(scan,seg)
		return ncesdiode

	def acount(self):
		"""Count number of uncut CESs in this acceptlist
			To avoid changing the accept list in any way this does not call 
			removeEmptyScans first, it just counts the included scans only.
			
			This may be unecessary but this makes the count non-destructive ie.e
			you could re-include a scan again later easily.
		"""
		used_scans = [ID for ID in self.scans if np.any(self.scans[ID]) ]
		return len(used_scans)
	#def CutTime(self,patch,flavour):
	#	"""Jackknife acceptlist in time"""
	#	import secateur
	#	return secateur.CutTime(self,flavour)

def buildL1Filename(l1dir,ttype,tname,scan,seg):
	return "%s/%s/%s/%s_%s_%d_seg%.3d.unf" % (l1dir,ttype,tname,ttype,tname,scan,seg)

def buildL2Filename(l2dir,ttype,tname,scan,seg):
	return "%s/%s/%s/%s_%s_%d_seg%.3d.unf" % (l2dir,ttype,tname,ttype,tname,scan,seg)

paramType_COMMON=1
paramType_TODPROC=2
paramType_TOD2MAP=3
paramType_TODSIM=4

def fortranfloat(x):
	try:
		return float(x)
	except:
		pass
	return float(x.replace('d','e'))


class ParameterDataClassic(object):
	"""
	
	JAZ This is for the old style parameter files which used parameters to indicate that other parameter files should be included.
	The new style files use the @INCLUDE syntax.  These should be read with the ParameterData class.
		
	A class that contains all the parameters from a set of parameter files - a common file, a tod_proc file, a tod2map file and a todsim file.
	It keeps track of the values of all the parameters and also which type of parameter file they came from, so you can write them 
	to the full set of files again.
	
	All the additional files are optional, and will only be read if the name is set in the common file.
	
	Initialize this with a common file and it will read the names of the other files from that and so read all the parameters.
	If you want to write out to a new set of files you need to set the values of the parameters:
	PARFILE_TOD_PROC, PARFILE_TOD2MAP, PARFILE_TODSIM
	using:
	parameters.setParam("PARFILE_TOD2MAP",paramType_COMMON,"new_filename.txt")
	etc.
	
	"""
	def __init__(self,commonFile):
		self.paramVals={}
		self.paramType={}
		self.getParams(commonFile,paramType_COMMON)
		if "PARFILE_TOD_PROC" in self.paramVals and os.path.exists(self.paramVals['PARFILE_TOD_PROC']):
			self.getParams(self.paramVals["PARFILE_TOD_PROC"],paramType_TODPROC)
		if "PARFILE_TOD2MAP" in self.paramVals and os.path.exists(self.paramVals['PARFILE_TOD2MAP']):
			self.getParams(self.paramVals["PARFILE_TOD2MAP"],paramType_TOD2MAP)
		if "PARFILE_TODSIM" in self.paramVals and os.path.exists(self.paramVals['PARFILE_TODSIM']):
			self.getParams(self.paramVals["PARFILE_TODSIM"],paramType_TODSIM)
			
	def __getitem__(self,p):
		return self.paramVals[p]
	def __setitem__(self,p,v):
		self.paramVals[p]=v
	def __delitem__(self,p):
		del self.paramVals[p]
	
	def setParam(self,name,ptype,value):
		self.paramType[name]=ptype
		self.paramVals[name]=value
		
	def toFiles(self,outfilename):
		f=open(outfilename,"w")
		for param,val in sorted(self.paramVals.iteritems()):
			if self.paramType[param]==paramType_COMMON:
				f.write("%s = %s\n" % (param,self.param2txt(val)))
		f.close()
		if "PARFILE_TOD_PROC" in self.paramVals:
			f = open(self.paramVals["PARFILE_TOD_PROC"],"w")
			for param,val in sorted(self.paramVals.iteritems()):
				if self.paramType[param]==paramType_TODPROC:
					f.write("%s = %s\n" % (param,self.param2txt(val)))
			f.close()
		if "PARFILE_TOD2MAP" in self.paramVals:
			f = open(self.paramVals["PARFILE_TOD2MAP"],"w")
			for param,val in sorted(self.paramVals.iteritems()):
				if self.paramType[param]==paramType_TOD2MAP:
					f.write("%s = %s\n" % (param,self.param2txt(val)))
			f.close()
			

	def getParams(self,filename,ptype):
		for line in open(filename):
			line=line.strip()
			if line=="" or line.startswith("#"):
				continue
			if line.startswith('@INCLUDE'):
				raise ValueError("You have tried to use a parameter file with @INCLUDE directives in the ParameterDataClassic class.  Use the ParameterData class for this.")
			data=line.split("=")
			param_name=data[0].strip()
			param_val_raw = ('='.join(data[1:])).strip()
			param_val = self.txt2param(param_val_raw)
			self.setParam(param_name,ptype,param_val)

	@staticmethod
	def param2txt(p):
		if p is True:
			return ".true."
		elif p is False:
			return ".false."
		elif isinstance(p,str):
			return "'%s'" % p
		else:
			return str(p)
	@staticmethod
	def txt2param(txt):
		txt=txt.strip()
		if txt==".true.":
			return True
		elif txt==".false.":
			return False
		else:
			try:
				p=int(txt)
				return p
			except:
				pass
			try:
				p=fortranfloat(txt)
				return p
			except:
				pass
			return txt.strip().strip("'")



import warnings
class ParameterData(object):
	"""
	JAZ This is the new style ParameterData object which works with @INCLUDE parameters.
	It will not follow the old style PARFILE_TOD2MAP type parameters.  If you need both of these at once let me know.
	
	A class that contains all the parameters from a set of parameter files.
	It keeps track of the values of all the parameters and also which parameter file they came from, so you can write them 
	to the full set of files again.
	
	If you try anything stupid like files included from multiple places or infinite inclusion loops then you'll get what you deserve.
	"""
	def __init__(self,filename):
		self.baseFilename=filename
		self.paramVals={}
		self.filenames=[]
		self.paramFiles={}
		self.inclusions={}
		self.getParams(filename)
		for param in ["PARFILE_TOD2MAP","PARFILE_TODSIM","PARFILE_TOD_PROC"]:
			if param in self.paramVals:
				warnings.warn("WARNING:  You have loaded a parameter file with the parameter %s set.  This parameter will not be followed by the new ParameterData class you used.  If you need the old behavior use the class ParameterDataClassic." % param)
	
	def filenameFor(self,param):
		"""Return the filename that the given parameter will be written to (the one it came from, unless you changed it)"""
		return self.filenames[self.paramFiles[param]]
	
	def renameParameterFileFor(self,param,newname):
		"""Rename the filename assoicated with a particular parameter, so that the parameter and all the ones from the same file will be saved to the new filename"""
		oldname = self.filenameFor(param)
		self.renameParameterFile(oldname,newname)
		
	def renameParameterFile(self,oldname,newname):
		"""Rename a parameter file so that parameters read from it will be saved to a new file."""
		filenumber = self.filenumberForFilename(oldname)
		self.filenames[filenumber]=newname
	
	def __getitem__(self,p):
		return self.paramVals[p]
	def __setitem__(self,p,v):
		self.setParam(p,v)
	def __delitem__(self,p):
		del self.paramVals[p]
		del self.paramFiles[p]

	def addInclude(self,parent,child):
		"""
		Create a new included file child from the parent file.
		The parent can be a number or name.
		"""
		parent_number = self.filenumberForFilename(parent)
		self.filenames.append(child)
		child_number = self.filenames.index(child)
		if parent_number not in self.inclusions:
			self.inclusions[parent_number] = []
		self.inclusions[parent_number].append(child_number)
		

	def filenumberForFilename(self,filename):
		if isinstance(filename,int):
			if filenumber<0 or filenumber >= len(self.filenames):
				raise ValueError("The file  number you specified is not in the current list.  To add a new one use addInclude")		
			return filename
		try:
			filenumber = self.filenames.index(filename)
		except ValueError:
			raise ValueError("The filename you specified is not in the current filename list.  To add a new one use addInclude")		
		return filenumber
		
	def setParam(self,name,value=None,file=None):
		"""
		Set a new value or file (or both) for a parameter, or create a new parameter with a new value and file.
		You can use either a filenumber or a filename for the file argument.
		
		If the parameter name does not already exist you must specify a value and file.
		"""
		if file is not None:
			filenumber = self.filenumberForFilename(file)
			
		if name in self.paramVals:
			if file is not None: self.paramFiles[name] = filenumber
			if value is not None: self.paramVals[name]=value
		else:
			if file is None:
				raise ValueError("When creating a new parameter you must specify which existing file it is to go in with the filename_or_number argument.  To add a new file use addInclude")
			if value is None:
				raise ValueError("You must give a parameter a value when you create it.")
			self.paramVals[name]=value
			self.paramFiles[name]=filenumber
			
			
			
	def _setParam(self,name,value,filenumber):
		self.paramVals[name]=value
		self.paramFiles[name]=filenumber

	def save(self):
		"""
		Save the parameter set to the currently filenames.
		"""
		self.saveTo(self.filenames)
	
	def saveAs(self,filenames):
		"""
		Save the parameter set to the given filenames and use them as the default filenames from now on.
		Pass in a list of filenames.
		"""
		self.filenames=list(filenames)
		self.save()

	def saveTo(self,filenames):
		"""
		Save the parameter set to the given filenames, but do not set them as the default filenames for use from now on.
		Pass in a list of filenames.
		"""
		if len(filenames) != len(self.filenames):
			error_message = "When using the saveAs command you need to supply the same number of filenames as there are parameter filenames, in this case %d" % len(self.filenames)
			raise ValueError()
		output_files = [open(filename,"w") for filename in filenames]
		for includer,included_list in self.inclusions.items():
			for included in included_list:
				output_files[includer].write("@INCLUDE '%s'\n" % filenames[included])
		for paramName in sorted(self.paramVals.keys()):
			filenumber = self.paramFiles[paramName]
			paramVal = self.paramVals[paramName]
			output_files[filenumber].write('%s = %s\n' % (paramName,self.param2txt(paramVal))  )

	def getParams(self,filename,includer=None):
		if filename not in self.filenames:
			self.filenames.append(filename)
		filenumber = self.filenames.index(filename)
		if includer is not None:
			if includer not in self.inclusions:
				self.inclusions[includer]=[]
			self.inclusions[includer].append(filenumber)
		for line in open(filename):
			line=line.strip()
			if line=="" or line.startswith("#"):
				continue
			if line.startswith('@INCLUDE'):
				if '#' in line:
					raise ValueError("No support for # comments after @INCLUDE directives in quiet.py.  Ask Joe if you need this.")
				filename=' '.join(line.split()[1:])
				filename = filename.strip("'")
				self.getParams(filename,includer=filenumber)
				continue
			data=line.split("=")
			param_name=data[0].strip()
			param_val_raw = ('='.join(data[1:])).strip()
			param_val = self.txt2param(param_val_raw)
			self._setParam(param_name,param_val,filenumber)

	@staticmethod
	def param2txt(p):
		if p is True:
			return ".true."
		elif p is False:
			return ".false."
		elif isinstance(p,str):
			return "'%s'" % p
		else:
			return str(p)
	@staticmethod
	def txt2param(txt):
		txt=txt.strip()
		if txt==".true.":
			return True
		elif txt==".false.":
			return False
		else:
			try:
				p=int(txt)
				return p
			except:
				pass
			try:
				p=fortranfloat(txt)
				return p
			except:
				pass
			return txt.strip().strip("'")





class Any_:
	"""
	An instance of Any simply responds "True" to "x in Any()" for any x.
	when print it yields the word 'all'
	Don't try to do anything else with it.
	"""
	def __contains__(self,x):
		return True
	def __str__(self):
		return 'all'

Any=Any_()


class Cut(object):
	"""A single cut, corresponding to a line in a cutfile"""
	def __init__(self, start,end,runs,segments,modules,diodes,comments):
		self.start = start
		self.end = end
		self.runs = runs
		self.segments = segments
		self.modules = modules
		self.diodes = diodes
		self.comments = comments
	def __str__(self):
		return "cut %f %f %s %s %s %s %s" % (
			self.start,self.end,
			printRange(self.runs),printRange(self.segments),printRange(self.modules),printRange(self.diodes),self.comments)
		
def parseRange(word):
	"""
	Parse a range from the format specified in the cut lists.
	Examples of the format:
	1
	1-3
	1,2,5
	1,3-4,9
	all
	none
	
	Returns a list of the numbers specified, or the special type Any which responds to "x in Any" with True for all x.
	Jon seems to have extended this to non-integer numbers.
	"""
	import re
	#if ',' in word and '-' not in word:
	#	return [int(i) for i in word.split(',')]
	if ',' in word or '-' in word:
		wordlist=[]
		blocks=word.split(',')
		for block in blocks:
			if re.search('-',block) == None:
				wordlist.append(int(block))
			else:
				ends=block.split('-')
				#for end in ends:
				#	end=int(round(float(end)))
				for end in range(int(round(float(ends[0]))),int(round(float(ends[-1])))+1):
	# *** Round CES.seg to nearest run number - NB BIT DODGY ***
					wordlist.append(int(round(end)))
		return wordlist
	elif word=='all':
		return Any
	elif word=='none':
		return []
	else:
		return [int(word)]

def printRange(R):
	if R is Any:
		return 'all'
	elif R is []:
		return 'none'
	else:
		return ','.join([str(r) for r in R])

class CutList(InOutFiles):
	def __init__(self,*infiles):
		self.cuts=[]
		for infile in infiles:
			self.infile=infile
			for line in self.infile:
				if (not line.startswith('#')) and line.strip():
					self.addCutLine(line,insert=False)
			del self.infile
		self.cuts.sort() #Sort the cuts by t1 then t2
	def addCutLine(self,line,insert=True):
		words=line.split()
		if (not words) or (not words[0])=='cut':return
		t1=float(words[1])
		t2=float(words[2])
		runs=parseRange(words[3])
		segs=parseRange(words[4])
		mods=parseRange(words[5])
		diodes=parseRange(words[6])
		comments = ' '.join(words[7:])
		cut = Cut(t1,t2,runs,segs,mods,diodes,comments)
		if insert:
			bisect.insort(self.cuts,cut)
		else:
			self.cuts.append(cut)

	def shouldCut(self,segment,module,diode,verbose=False):
		startDate=segment.startDate
		endDate=segment.endDate
		shouldCut=False
		for cut in self.cuts:
#JAZ get these lines to work some time, if the code is too slow.  Sorting is wrong, I think
#			if cut.end < startDate: continue  #The cuts are sorted by start date then end date
#			if cut.start > endDate: break     #so these two lines skip all the irrelevant ones
			if (cut.start<startDate<cut.end) or (cut.start<endDate<cut.end) or (cut.start<startDate and cut.end>endDate) :
				if (module in cut.modules) and (diode in cut.diodes):
					if verbose:
						print "Cutting segment (%d,%d), diode (%d,%d), reason: %s" % (segment.runID,segment.id,module,diode,cut.comments)
					shouldCut=True
					break
		return shouldCut
	
	def pruneAcceptList(self,runlist,accept,verbose=False):
		segments = runlist.segments()
		for (scan_id,segment_id) in accept.scans:
			segment = segments[(scan_id,segment_id)]
			for module in xrange(N_MODULE_MAX):
				for diode in xrange(NDIODE_MODULE):
					if self.shouldCut(segment,module,diode,verbose):
						accept.excludeSegmentDiode(scan_id,segment_id,module,diode)


class L2Info(object):
	def __init__(self,filename,nmod,nsamp,time,modules):
		self.filename=filename
		self.nModules=nmod
		self.nSamples=nsamp
		self.time=time
		self.modules=np.array(modules)
	def __repr__(self):
		return "<L2 File %s: %d modules, %d samples, time:%f - %f>" % (self.filename,self.nModules,self.nSamples,self.time[0],self.time[-1])

def l2_info(filename):
	f=FortranFile(filename)
	nmod=f.readArray(np.int32)
	nsamp=f.readArray(np.int32)
	time=f.readArray(np.float64)
	module_numbers=[]
	for i in xrange(nmod):
		module_number=f.readArray(np.int32)
		module_numbers.append(module_number)
		f.advance(2)
	f.close()
	return L2Info(filename,nmod,nsamp,time,module_numbers)


	

class ModuleList(list):
	def __init__(self,nsamp,reduction_factor,num_reduced,time=None,az=None,el=None,dk=None):
		list.__init__(self)
		self.nsamp=nsamp
		self.reduction_factor=reduction_factor
		self.num_reduced=num_reduced
		self.time=time
		self.az=az
		self.el=el
		self.dk=dk
	def __repr__(self):
		if self.time is not None:
			return "<A list of %d module data length %d starting at MJD %f>" % (len(self),self.nsamp,self.time[0])
		else:
			return "<A list of %d module data length %d>" % (len(self),self.nsamp)

def l2_read(filename,modules=None,demod=True,pointing=True,TP=False,time=True,horizon=False,maps=True,full=False):
	"""
	Load from a level 2 quiet file.

	Various parameters control which data is read from the file.  The default is just to read the demod data and module pointing, 
	but I change that if it's annoying.
	
	This function now returns a subclass of the standard python list, with additional attributes containing the file-level
	data. Always: nsamp,num_reduced,reduction_factor
		  If you ask for them: time,az,el,dk (otherwise set to None) 

	Required parameter:
		filename - the full path of the level 2 file
	Optional parameters:
		modules - a list (or other sequence) of modules to read (others are ignored). Default: read all.
		
		demod [default: True] - set to False to ignore demod info and load only total power.
		pointing [default: True] - set to False to ignore pointing data and set it to None.
		maps [default: True] - set to False to ignore the maps that now come with L2
		TP [default: False] - set to True to load the total power and variance (one array per module).
		time [default: True] - set to True to read the time array (one array for whole file).
		horizon [default: False] - set to True to read the original horizon co-ordinate (az,el,dk).  One set for whole file.
		
		full [defaut: False] - set all these flags to True (ignore the passed flags)
	"""
	if full:
		demod=True
		pointing=True
		TP=True
		time=True
		horizon=True
	f=FortranFile(filename)
	nmod=f.readArray(np.int32)
	nsamp=f.readArray(np.int32)
	reduction_factor=f.readArray(np.int32)
	num_reduced=f.readArray(np.int32)
	nside=f.readArray(np.int32)
	module_numbers=f.readArray(np.int32)
	if modules is None:
		modules=Any
	map_data=[]
	npix=f.readArray(np.int32)  #Used only to get sizes for fortran - do not need to save.
	for i in xrange(nmod):
		if maps and (module_numbers[i] in modules):
			mask2map = f.readArray(np.int32)
			nhit = f.readArray(np.int32)
			l2map = L2Map(nhit,mask2map,nside)
			map_data.append(l2map)
		else:
			f.advance(3)			
			map_data.append(None)

	if time:
		#		time=f.readArray(np.float32)
		time=f.readArray(np.float64)
	else:
		time=None
		f.advance()
	if horizon:
		az,el,dk=f.readArray(np.float32).reshape((nsamp,3)).transpose()
	else:
		f.advance()
		az,el,dk=None,None,None
	output=ModuleList(nsamp,reduction_factor,num_reduced,time,az,el,dk)
	#read pointing
	#import pdb;pdb.set_trace()
	pointings = []
	for i in xrange(nmod):
		if pointing and (module_numbers[i] in modules):
			phi,theta,psi=f.readArray(np.float32).reshape((nsamp,3)).transpose()
			pointing_data=Pointing(theta,phi,psi)
			pointings.append(pointing_data)
		else:
			f.advance()
			pointings.append(None)

	#Read TOD
	tods=[]
	for i in xrange(nmod):
		if demod and (module_numbers[i] in modules):
			q1,u1,u2,q2=f.readArray(np.float32).reshape((nsamp,4)).transpose()
			tod=TOD(q1,u1,u2,q2)
			tods.append(tod)
		else:
			f.advance()
			tods.append(None)
	
	#Read TP
	tp_data=[]	
	for i in xrange(nmod):
		if TP and (module_numbers[i] in modules):
			tp_q1,tp_u1,tp_u2,tp_q2=f.readArray(np.float32).reshape((num_reduced,4)).transpose()
			tp=TOD(tp_q1,tp_u1,tp_u2,tp_q2)
			tpv_q1,tpv_u1,tpv_u2,tpv_q2=f.readArray(np.float32).reshape((num_reduced,4)).transpose()
			tp_var=TOD(tpv_q1,tpv_u1,tpv_u2,tpv_q2)
			tp_data.append((tp,tp_var))
		else:
			f.advance(2)
			tp_data.append((None,None))

	#Collect together
	for (module_number,pointing_data,tod,(tp,tp_var),map) in zip(module_numbers,pointings,tods,tp_data,map_data):
		if module_number in modules:
			module=Module(module_number,nsamp,reduction_factor,num_reduced,pointing_data,tod,time,tp,tp_var,map)
			output.append(module)
	f.close()
	return output

def l2_read_old(filename,modules=None,demod=True,pointing=True,TP=False,time=True,horizon=False,full=False):
	"""
	Load from a level 2 quiet file.

	Various parameters control which data is read from the file.  The default is just to read the demod data and module pointing, 
	but I change that if it's annoying.

	This function now returns a subclass of the standard python list, with additional attributes containing the file-level
	data. Always: nsamp,num_reduced,reduction_factor
		  If you ask for them: time,az,el,dk (otherwise set to None) 

	Required parameter:
		filename - the full path of the level 2 file
	Optional parameters:
		modules - a list (or other sequence) of modules to read (others are ignored). Default: read all.

		demod [default: True] - set to False to ignore demod info and load only total power.
		pointing [default: True] - set to False to ignore pointing data and set it to None.
		TP [default: False] - set to True to load the total power and variance (one array per module).
		time [default: True] - set to True to read the time array (one array for whole file).
		horizon [default: False] - set to True to read the original horizon co-ordinate (az,el,dk).  One set for whole file.

		full [defaut: False] - set all these flags to True (ignore the passed flags)
	"""
	if full:
		demod=True
		pointing=True
		TP=True
		time=True
		horizon=True
	f=FortranFile(filename)
	nmod=f.readArray(np.int32)
	nsamp=f.readArray(np.int32)
	reduction_factor=f.readArray(np.int32)
	num_reduced=f.readArray(np.int32)
	if time:
		time=f.readArray(np.float64)
	else:
		time=None
		f.advance()
	if horizon:
		az,el,dk=f.readArray(np.float64).reshape((nsamp,3)).transpose()
	else:
		f.advance()
		az,el,dk=None,None,None
	output=ModuleList(nsamp,reduction_factor,num_reduced,time,az,el,dk)
	for i in xrange(nmod):
		module_number=f.readArray(np.int32)
		if demod and ((modules is None) or (module_number in modules)):
			if pointing:
				phi,theta,psi= f.readArray(np.float64).reshape((nsamp,3)).transpose()
				pointing_data=Pointing(theta,phi,psi)
			else:
				f.advance()
				pointing_data=None
			if demod:
				q1,u1,u2,q2=f.readArray(np.float64).reshape((nsamp,4)).transpose()
				tod=TOD(q1,u1,u2,q2)
			else:
				tod=None
				f.advance(1)
			if TP:
				tp_q1,tp_u1,tp_u2,tp_q2=f.readArray(np.float64).reshape((num_reduced,4)).transpose()
				tp=TOD(tp_q1,tp_u1,tp_u2,tp_q2)
				tpv_q1,tpv_u1,tpv_u2,tpv_q2=f.readArray(np.float64).reshape((num_reduced,4)).transpose()
				tp_var=TOD(tpv_q1,tpv_u1,tpv_u2,tpv_q2)
			else:
				tp=None
				tp_var=None
				f.advance(2)
			module=Module(module_number,nsamp,reduction_factor,num_reduced,pointing_data,tod,time,tp,tp_var)
			output.append(module)
		else:
			f.advance(4)
	f.close()
	return output
	

def l2_read_time(filename):
	"""
	Reads the time data array contained in the level 2 file called filename.  This is faster than reading the whole file if you do not need 
	the rest of the data.
	Returns an array.
	"""
	nmod=f.readArray(np.int32)
	nsamp=f.readArray(np.int32)
	reduction_factor=f.readArray(np.int32)
	num_reduced=f.readArray(np.int32)
	time=f.readArray(np.float64)
	return time
	

def check_for_none(value,name):
	if value is None:
		raise ValueError("No data '%s' in data to be used in l2_write - must have full data to write L2 file" % name)

def l2_write(modules,filename):
	raise NotImplemented("Have not yet written code to l2_write new file format.  Please ask if you require.")

def l2_write_old(modules,filename):
	"""
	Write a level 2 file to "filename" from a list of modules.
	"""
	try:
		nmod = np.array([len(modules)],dtype=np.int32)
		nsamp = modules.nsamp
		reduction_factor=modules.reduction_factor
		num_reduced=modules.num_reduced
		time=modules.time
		az=modules.az
		el=modules.el
		dk=modules.dk
	except AttributeError, (A):
		raise AttributeError(str(A) + "\nYou should probably be calling l2_write using the new ModuleList class")
	check_for_none(time,"time")
	check_for_none(az,"az")
	check_for_none(el,"el")
	check_for_none(dk,"dk")
	for m in modules:
		check_for_none(m.tod,"tod")
		check_for_none(m.pointing,"pointing")
		check_for_none(m.tp,"tp")
		check_for_none(m.tp_var,"tp_var")
		
	f=FortranFile(filename,"w")
	horizon = np.array([az,el,dk]).transpose()
	f.writeArrays([nmod,nsamp,reduction_factor,num_reduced,time,horizon])
	for m in modules:
		f.writeArray(m.number)
		pointing=np.array([m.pointing.phi,m.pointing.theta,m.pointing.psi]).transpose()
		f.writeArray(pointing)
		tod=np.array([m.tod.Q1,m.tod.U1,m.tod.U2,m.tod.Q2]).transpose()
		f.writeArray(tod)
		tp=np.array([m.tp.Q1,m.tp.U1,m.tp.U2,m.tp.Q2]).transpose()
		f.writeArray(tp)
		tp_var=np.array([m.tp_var.Q1,m.tp_var.U1,m.tp_var.U2,m.tp_var.Q2]).transpose()
		f.writeArray(tp_var)
	f.close()
		

class Calibration(InOutFiles):
	def __init__(self,infile):
		self.infile=infile
		self.params={}
		for line in self.infile:
			self.process_line(line)
		self.sort()
		del self.infile

	def sort(self):
		for param in self.params.itervalues():
			for val in param.itervalues():
				if type(val)==list:
					val.sort()
					
		
	def __getattr__(self,name):
		if hasattr(self,"params"):
			if name in self.params:
				return self.params[name]
		raise AttributeError, name
		
	def process_line(self,line):
		words=tuple(numerify(line.split()))
		item,startDate,endDate,scan,seg,mod,diode,value=words[:8]
		if item=='corr':
			return #For now.
		if item not in self.params:
			self.params[item]={}
			def get_item(self):
				return self.params[item]
		data=self.params[item]
		if scan==-1:
			#This item is indexed by date
			if (mod,diode) not in data:
				data[(mod,diode)]=[]
			data[(mod,diode)].append((startDate,endDate,value))
#			data[(mod,diode)].sort()
		else:
			#This item is indexed by scan number
			if (scan,seg) not in data:
				data[(scan,seg)]=np.zeros((N_MODULE_MAX,NDIODE_MODULE)) + np.nan
			scan_data=data[(scan,seg)]
			scan_data[mod,diode]=value
	def lookup_date(self,item,module,diode,date,default=np.nan):
		item_data=getattr(self,item)
		try:
			diode_data=item_data[(module,diode)]
		except:
			return default
		for (start,end,value) in diode_data:
			if start<date<end:
				return value
		#otherwise linearly interpolate
		for i in xrange(len(diode_data)-1):
			first_end=diode_data[i][1]
			next_start=diode_data[i+1][0]
			if first_end<date<next_start:
				r=(date-first_end)/(next_start-first_end)
				return r*diode_data[i][2] + (1-r)*diode_data[i+1][2]
		return default
	def lookup_scan(self,item,module,diode,scan,seg,default=np.nan):
		item_data=getattr(self,item)
		try:
			result = item_data[(scan,seg)][module,diode]
		except:
			result = default
		if result==np.nan:
			return default
		return result


def read_pixfile(filename):
	f=open(filename)
	lines=[line.strip() for line in f.xreadlines()]
	nside=int(lines[0].split("=")[1])
	npix=int(lines[1].split("=")[1])
	pix=np.zeros(npix,dtype=int)
	vals=np.zeros(npix,dtype=float)
	for i,line in enumerate(lines[2:]):
		p,x=line.split()
		p=int(p)
		x=float(x)
		pix[i]=p
		vals[i]=x
	f.close()
	return nside,pix,vals



def write_fits_partial(nside,pix,vals,filename,nest=False):
	col1=pyfits.Column(name="PIXEL",format="1J",array=pix)
	col2=pyfits.Column(name="RMSNOISE",format="1E",array=vals)
	coldefs=pyfits.ColDefs([col1,col2])
	tbhdu = pyfits.new_table(coldefs)
	tbhdu.header.update('PIXTYPE','HEALPIX','HEALPIX pixelisation')
	if nest:
		tbhdu.header.update('ORDERING','NESTED','Pixel ordering scheme, either RING or NESTED')
	else:
		tbhdu.header.update('ORDERING','RING','Pixel ordering scheme, either RING or NESTED')

	tbhdu.header.update('EXTNAME','DATA','name of this binary table extension')
	tbhdu.header.update('NSIDE',nside,'Resolution parameter of HEALPIX')
	tbhdu.header.update('GRAIN',1,'Indexing is explicit and is given by PIXEL col')
	tbhdu.header.update('OBJECT','PARTIAL','Sky coverage is partial')
	tbhdu.header.update('INDXSCHM','EXPLICIT',
						'Indexing: IMPLICIT or EXPLICIT')
	tbhdu.header.update('OBS_NPIX',len(pix),
						'Number of observed pixels')
	tbhdu.writeto(filename,clobber=True)

