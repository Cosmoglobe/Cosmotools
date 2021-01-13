#!/usr/bin/python

"""
Jackknife acceptlist using .cut lists plus cutlist from (e.g.) weather cut\n
Call from within python thus:
    import makeJackknives
    makeJackknives.makeJackknife(cutces,flavours,jack_outdir='.',cut_indir='.',deadfiles=[]):
"""

#-----------------------------------------------------------

def makeJackknife(cutces,flavours,jack_outdir='.',cut_indir='.',deadfiles=[]):

      import quiet,secateur

      """
      Jackknife acceptlist using input .cut lists and weather-cut list
      """

      patch='patch_6a'

      runlist,fulllist=secateur.GenerateCMBRunlists()
      remaining={}
      cutlist={}
      cutlist[patch]=quiet.AcceptedScanList(cutces)

      print 'number of CESes in full list', fulllist[patch].acount()
      print 'number of CESes in cut list', cutlist[patch].acount()

      remaining[patch]=secateur.multizap(fulllist[patch],cutlist[patch])

      for deadfile in deadfiles:
            dead=quiet.CutList(deadfile)
      #quiet.CutList.pruneAcceptList(dead,runlist['patch_6a'],fulllist['patch_6a'])
            dead.pruneAcceptList(runlist[patch],fulllist[patch])

      print 'number of CESes in runlists', remaining[patch].acount()

      #cutfile1, cutfile2 = secateur.CutTime(remaining,flavour,jack_outdir)
      #cutfile1, cutfile2 = secateur.ApplySingleJackknife(remaining,flavour,jack_outdir)
      # Working here
      # Need to work in optional accept/runlists
      secateur.ApplyMultipleOperations(flavours,patch)
      #remaining.CutTime(flavour)

      #return cutfile1, cutfile2
      return

#-----------------------------------------------------------

##############################
#  Main program starts here  #
##############################

if __name__ == '__main__':

      import sys

      if len(sys.argv) <= 4:
            print >> sys.stderr, '%s <cutlist> <list of jackknifes> <jack_outdir> <cut_indir> <.cut list of dead diodes.cut> ... <.cut list of bad times>' % sys.argv[0]
            sys.exit(2)

      cutces,flavours,jack_outdir,cut_indir = sys.argv[1:4]
      #flavour = sys.argv[2]
      #jack_outdir = sys.argv[3]
      #cut_indir = sys.argv[4]

      deadfiles = []
      for arg in sys.argv[5:]:
            deadfiles.append(cut_indir+arg)

      makeJackknife(cutces,flavours,jack_outdir,cut_indir,deadfiles)

      sys.exit(0)

#-----------------------------------------------------------
