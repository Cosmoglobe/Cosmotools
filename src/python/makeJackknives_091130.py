#!/usr/bin/python

"""
Jackknife acceptlist using .cut lists plus cutlist from (e.g.) weather cut\n
Call from within python thus:
    import makeJackknives
    makeJackknives.makeJackknife(cutces,flavours,jack_outdir='.',cut_indir='.',deadfiles=[],patch='patch_6a'):
"""

#-----------------------------------------------------------

def makeJackknife(cutces,jack_outdir='.',cut_indir='.',deadfiles=[],patch=None,*flavours):

      import quiet,secateur,re,copy,sys

      """
      Jackknife acceptlist using input .cut lists and weather-cut list
      """

      # Handle requested jackknives
      if len(flavours)==1 and type(flavours)==dict:
            pass
      #elif len(flavours)==1 and type(flavours)==str:
      else:
            flavia={}
            jtype,flavia[jtype]=re.split(' ',flavours[0])
            flavours=copy.deepcopy(flavia)
      #     print 'Choking in makeJackknife, because of flavour defn', flavours
      #     sys.exit(1)

      # Fetch master runlist and acceptlist
      runlist,fulllist=secateur.GenerateCMBRunlists(patch=patch)

      # Remove irrelevant patches, if none is supplied explicitly
      if patch is not None:
            runlist[patch]=secateur.CleanL2Runlist(runlist[patch],patch)

      # Read the 'weather' acceptlist cut
      cutlist={}
      cutlist[patch]=quiet.AcceptedScanList(cutces)

      print 'number of CESs in full list', fulllist[patch].acount()
      print 'number of CESs in cut list',  cutlist[patch].acount()

      # Apply the 'weather' acceptlist cut to the active acceptlist
      remaining={}
      remaining[patch]=secateur.multizap(fulllist[patch],cutlist[patch])

      # Apply the static and other cuts to the active acceptlist
      for deadfile in deadfiles:
            dead=quiet.CutList(deadfile)
            #dead.pruneAcceptList(runlist[patch],fulllist[patch]) # bug?
            # JTLZ TO CHECK THIS LINE!!!!!
            dead.pruneAcceptList(runlist[patch],remaining[patch])

      print 'number of CESs in runlists', remaining[patch].acount()

      # Loop over types of jackknife
      for jtype in flavours.keys(): 
           split=jtype+' '+flavours[jtype]
           aname,bname=secateur.ApplySingleJackknife(split,jack_outdir,remaining,runlist)

      return aname,bname

#-----------------------------------------------------------

##############################
#  Main program starts here  #
##############################

if __name__ == '__main__':

      import sys,os

      narg=4
      if len(sys.argv) <= narg:
            print >> sys.stderr, '%s <cutlist> <patch> <jack_outdir> <cut_indir> <.cut list of dead diodes.cut> ... <.cut list of bad times>' % sys.argv[0]
            sys.exit(2)

      patch='patch_6a'

      cutces,patch,jack_outdir,cut_indir = sys.argv[1:narg+1]

      deadfiles = []
      for arg in sys.argv[narg+1:]:
            deadfiles.append(os.path.join(cut_indir,arg))

      flavours={}
      flavours['time'] = ['season','battenberg','interleaved','tandem']
      flavours['module'] = ['alternate','radius','horizontal','vertical1','vertical2','stripes1','stripes2','stripes3','mab1','mab2','mab3','mab12','mab13','mab23','diodes','frequency']
#cut1, cut2 = makeJackknives_091130.makeJackknife(cutlist, jackknife, '.','.',['dead.cut','static.cut'])
      makeJackknife(cutces,flavours,jack_outdir,cut_indir,deadfiles,patch)

      sys.exit(0)

#-----------------------------------------------------------
