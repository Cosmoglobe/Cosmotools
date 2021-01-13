#!/usr/bin/python

"""
Set up map2cl directory and files
"""

import os
user = os.getenv('USER')
SVN='/usit/titan/u1/%s/quiet_svn/oslo/src/f90' % user
EXE_MAP2CL=os.path.join(SVN,'map2cl/map2cl')
DATA_DIR='/data4/quiet/%s/' % user
SETUP_DIR='/data4/quiet/jonathaz/map2cl_setup'
MAX_NCORES_PER_NODE=8

#-------------------------------------------------------------------

BoolPy2F={True:'.true.',False:'.false.'}

#-------------------------------------------------------------------

def generateSimulationList(dir,listf='sims.txt',expr='*lcut5_map.fits',\
                           noffset=1,ipatch=1):

    """
    For a dir containing tod2map lcut simulation realisations,
    listr these realisations to file and
    make the filename available for the parfile
    A chance to try out the pythonic glob
    """

    import os,glob

    cmbfiles=glob.glob(os.path.join(os.path.abspath(dir),expr))
    cmbfiles.sort()

    num_cmbfiles=len(cmbfiles)+noffset

    if os.path.exists(listf): os.rename(listf,'%s.old'%listf)
    list=open(listf,'w')

    for n,f in enumerate(cmbfiles):
        line="CMBFILE%02i_%04i='%s'" %(ipatch,n+1+noffset,f)
        #print line
        list.write('%s\n'%line)

    list.close()

    return os.path.join(dir,listf),num_cmbfiles

#-------------------------------------------------------------------

def writeSubstitutions(outf,substituted):
    """
    Write substituted template to file
    """

    if os.path.exists(outf): os.rename(outf,'%s.old'%outf)
    outfile=open(outf,'w')
    outfile.write(substituted)
    print 'Wrote %s' % outf
    outfile.close()

    return

#-------------------------------------------------------------------

def inject(patient,complaint,medicine):

    """
    Inject the substring medicine into the filename patient triggering
    on character complaint
    """

    import os

    path=os.path.dirname(patient)
    conv=list(os.path.basename(patient))
    conv.insert(conv.index(complaint),medicine)
    conv=os.path.join(path,''.join(conv))

    return conv

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,shutil
    from optparse import OptionParser

    print
    # Set up options
    usage = "./prepMap2Cl.py -h \n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
           help="specify QUIET CMB patch", type="string",\
           metavar='(2a|4a|6a|7b)')

    parser.add_option("-t", "--templatedir", dest="setupdir",\
           help="specify setup dir", type="string",\
           metavar=SETUP_DIR)

    parser.add_option("-f", "--parfile", dest="parfile",\
           help="specify template parfile", type="string",\
           metavar='parfile_template.txt')

    parser.add_option("-o", "--outdir", dest="outdir",\
           help="specify output dir for scripts", type="string",\
           metavar='output_dir')

    parser.add_option("-s", "--stage", dest="stage",\
           help="specify reduction stage", type="string",\
           metavar='[map2cl]')

    parser.add_option("-m", "--map", dest="mapf",\
           help="specify map", type="string",\
           metavar='map.fits')

    parser.add_option("-k", "--mask", dest="maskf",\
           help="specify maskf", type="string",\
           metavar='mask.fits')

    parser.add_option("-c", "--covar", dest="covarfs",\
           help="specify covariance matrix(es)", type="string",\
           metavar='covar.unf[,sqrtf,invsqrtf]')

    parser.add_option("-e", "--preexist", action="store_true", dest="exist",
                      help="If you have pre-existing sqrt, invsqrt covars")

    parser.add_option("-n", "--nside", dest="nside",\
           help="specify nside", type="int",metavar='128,256')

    parser.add_option("-u", "--nodes", dest="nodes",\
           help="specify number of nodes", type="int",\
           metavar='4')

    parser.add_option("-r", "--tod2map", type="string", dest="t2m",
                      help="dir for tod2map realisations",\
                      metavar="100000_sims")

    parser.add_option("-x", "--numcmbfiles", dest="num_cmbfiles_init",\
           help="specify num_cmbfiles", type="int",metavar='100')

    parser.add_option("-y", "--numsims", dest="num_sims_init",\
           help="specify num_sims", type="int",metavar='0')

    # Assert defaults
    parser.set_defaults(stage='map2cl',setupdir=SETUP_DIR,\
        patch=None,verbose=False,outdir='.',mapf=None,maskf=None,\
        covarfs=None,nodes=None,parfile=None,nside=128,exist=False,\
                        t2m=None,num_cmbfiles_init=None,num_sims_init=None)

    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    verbose=options.verbose
    setupdir=options.setupdir
    outdir=options.outdir
    stage=options.stage
    mapf=options.mapf
    maskf=options.maskf
    covarfs=options.covarfs
    nodes=options.nodes
    parfile=options.parfile
    nside=options.nside
    exist=options.exist
    t2m=options.t2m
    num_cmbfiles_init=options.num_cmbfiles_init
    num_sims_init=options.num_sims_init

    # Check we have all the inputs
    assert(patch is not None), 'patch required!'
    assert(mapf is not None), 'map required!'
    assert(maskf is not None), 'mask required!'
    assert(covarfs is not None), 'covars required!'
    #assert(npernode<=MAX_NCORES_PER_NODE), 'invalid npernode %f!' %npernode

    mapf=os.path.abspath(mapf)
    maskf=os.path.abspath(maskf)

    # Set up the covar files
    covarf=None; sqrtf=None; invsqrtf=None
    if ',' in covarfs:
        for c in covarfs.split(','):
            if 'sqrt' in c and 'inv' in c: invsqrtf=c
            if 'sqrt' in c and 'inv' not in c: sqrtf=c
            if 'sqrt' not in c and 'inv' not in c: covarf=c
            c=os.path.abspath(c)
    else:
        covarf=covarfs
        covarf=os.path.abspath(covarfs)

    if exist:
        sqrtf=inject(covarf,'N','sqrt_')
        invsqrtf=inject(covarf,'N','sqrt_inv_')

    # Check lcut has been applied!
    for c in [covarf,sqrtf,invsqrtf,mapf]:
        assert('lcut' in c), '%s not lcut!' % c

    print 'covars: ',covarf,sqrtf,invsqrtf

    # Set up the template parfile and batch file
    if parfile is not None:
        parfile=os.path.abspath(parfile)
    else:
        parfile='parfile_template.txt'
        parfile=os.path.join(setupdir,parfile)

    batchfile='batch.txt'
    batchfile=os.path.join(setupdir,batchfile)

    # Set up the output dir
    outdir=os.path.abspath(outdir)
    if not os.path.exists(outdir): os.mkdir(outdir)

    # Copy some files from the setup dir
    transfers=[t%nside for t in ['bins_d50_n%s.dat','pixel_window_n0%s.fits']]
    for f in transfers:
        fr=os.path.join(setupdir,f)
        to=os.path.join(outdir,f)
        if not os.path.exists(to): shutil.copy(fr,to)

    # If tod2map-supplied realisations
    # Working here
    if t2m is not None:
        t2m=os.path.abspath(t2m)
        cmbfiles,num_cmbfiles=generateSimulationList(t2m,noffset=1)
        include_cmbfiles=''
        if num_cmbfiles_init is not None:
            num_cmbfiles=num_cmbfiles_init
        numsims=0
    else:
        cmbfiles=''
        include_cmbfiles='#'
        num_cmbfiles=1
        numsims=100

    if num_sims_init is not None:
        numsims=num_sims_init

    # Set up the batch and par files, by substitution
    if stage=='map2cl':
        pass

    # Specify the cluster parameters
    exe=EXE_MAP2CL
    job='%s_%s' % (os.path.basename(outdir),patch)
    clock='96:00:00'
    if nodes is None: nodes = 8
    #ncores=104
    ncores=nodes*MAX_NCORES_PER_NODE
    npernode=-1 # not in use
    omp_num_threads=1
    #omp_num_threads=int(MAX_NCORES_PER_NODE/npernode)
    #assert(npernode*omp_num_threads==MAX_NCORES_PER_NODE),\
    #       'invalid threading combination'

    # Set up the output and working directories
    if outdir is not None:
        work_dir=os.path.abspath(outdir)
    else:
        work_dir=os.path.join(DATA_DIR,out_dir)

#       par: mapf,maskf,covarf,patch,setupdir
#       batch: outdir,exe,job,nodes,cores,omp_num_threads,users,workdir
    substitutions={'mapf':'map.fits',
                   'maskf':'mask.fits',
                   'covarf':'covar.unf',
                   'sqrtf':'sqrt.unf',
                   'invsqrtf':'invsqrt.unf',
                   'patch':patch,
                   'setupdir':setupdir,
                   'outdir':outdir,                       
                   'exe':exe,
                   'job':job,
                   'walllimit':clock,
                   'nodes':str(nodes),
                   'cores':str(ncores),
                   'omp_num_threads':str(omp_num_threads),
                   'users':user,
                   'workdir':work_dir,
                   'npernode':str(npernode),
                   'nside':str(nside),
                   'num_cmbfiles':str(num_cmbfiles),
                   'include_cmbfiles':include_cmbfiles,
                   'cmbfiles':cmbfiles,
                   'numsims':numsims
                   }

    # Make the substitutions
    files={'batch':'batch.txt','par':'parfile_n%i.txt'%nside}
    templatef={'batch':batchfile,'par':parfile}

    template={}; substituted={};
    for f in files.keys():
        templatef[f]=os.path.join(setupdir,templatef[f])
        template[f]=open(templatef[f],'r').read()
        substituted[f] = template[f] % substitutions
        if verbose:
            print "'ere we go:\n"
            print substituted[f]
        outf=os.path.join(outdir,files[f])
        writeSubstitutions(outf,substituted[f])

    # Link to the map, mask and covar files
    fdict={mapf:'map.fits',maskf:'mask.fits',covarf:'covar.unf',\
           sqrtf:'sqrt.unf',invsqrtf:'invsqrt.unf'}
    for f in fdict.keys():
        link=os.path.join(outdir,fdict[f])
        if os.path.exists(link):
            if os.path.islink(link) and (os.readlink(link) == f):
                pass
            else:
                print 'The symlink for %s already exists\
                    but points elsewhere (%s)' % (link,os.readlink(link))
                print 'Try deleting the existing link, and rerun'
                sys.exit(0)
        else:
            os.symlink(f,link)

    # Nearly finished...
    print 'Now do:'
    print 'cd %s' %outdir
    print 'sbatch batch.txt'

    sys.exit(0)
