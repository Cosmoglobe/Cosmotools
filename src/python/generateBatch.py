#!/usr/bin/python

"""
Take a tod2map parfile and generate batch scripts for finalmap and lcut

1. Identify jackknife pairs
2. Identify patch, run, dir
3. Fetch templates

"""
#-------------------------------------------------------------------

import os

user = os.getenv('USER')
SVN='/usit/titan/u1/%s/quiet_svn/oslo/src/f90' % user
EXE_FINALMAP=os.path.join(SVN,'scalapost/scalapost')
EXE_LCUT=os.path.join(SVN,'postmap/postmap')
DATA_DIR='/data4/quiet/%s/' % user

DEFAULT_ELL=5
#DEFAULT_ELL=25

if __name__ == '__main__':

    import sys,os
    import quiet,secateur
    from optparse import OptionParser
    import scanf

    print
    # Set up options
    usage = "./generateBatch.py -h \n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
           help="specify QUIET CMB patch", type="string",\
           metavar='(2a|4a|6a|7b)')

    parser.add_option("-f", "--parfile", dest="parfile",\
           help="specify input parfile", type="string",\
           metavar='parfile.txt')

    parser.add_option("-t", "--template", dest="templatef",\
           help="specify template batch", type="string",\
           metavar='template.txt')

    parser.add_option("-o", "--outdir", dest="outdir",\
           help="specify output dir for scripts", type="string",\
           metavar='output_dir')

    parser.add_option("-q", "--quick", action="store_true", dest="quick",
                      help="apply onlymap/nocov")

    parser.add_option("-s", "--stage", dest="stage",\
           help="specify reduction stage", type="string",\
           metavar='[finalmap|lcut|lcutsum|finalmapmask|\
            lcutmask|lcutsummask|lcutsims]')

    # Assert defaults
    parser.set_defaults(stage=None,quick=False,\
        patch=None,parfile=None,templatef=None,verbose=False,outdir='.')
    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    verbose=options.verbose
    parfile=options.parfile
    templatef=options.templatef
    outdir=options.outdir
    stage=options.stage
    quick=options.quick
    
    assert(patch is not None), 'patch required!'
    assert(parfile is not None), 'parfile required!'
    assert(stage is not None), 'stage required!'
    #assert(templatef is not None), 'template required!'

    #assert(stage is not 'lcutmask'), '%s not implemented' % stage
    assert(stage is not 'lcutsummask'), '%s not implemented' % stage

    # Hack here if filename is in non-standard form
    try:
        version,patch=scanf.sscanf(os.path.basename(parfile),'par-%7c-%2c.txt')
    except scanf.IncompleteCaptureError:
        version,patch,details=scanf.sscanf(os.path.basename(parfile),'par-%7c-%2c-%s')
        version='%s_%s' % (version,details.split('.')[0])

    #version='chic8'; patch='2a';details='multi';

    params=quiet.ParameterData(parfile)

    nsets=params['NUM_SETS']
    print nsets,version,patch

    accept={}
    pairs={}
    for n in range(1,1+nsets):
        accept[n]=params['ACCEPTED_SCANS%i'%n]
        flavour=accept[n].split('_')[-1].split('.')[0]
        if flavour not in pairs.keys(): pairs[flavour]=[]
        #pairs[flavour].append('%i'%n)
        pairs[flavour].append('%03i'%n)

    if outdir is not None:
        work_dir=os.path.abspath(outdir)
    else:
        out_dir=params['OUTPUT_DIR']
        work_dir=os.path.join(DATA_DIR,out_dir)

    if templatef is not None: template=open(templatef,'r').read()

    print 'splits are', pairs

    who='columbia'

    # Set up the point-source mask, if required
    mask='junk'
    if stage=='finalmapmask':
        nside=params['NSIDE_OUT']
        #maskdir='/data4/quiet/external_data'
        #maskf='wmap_temperature_source_mask_r9_5yr_v3_n%s_quiet.fits'%nside
        maskdir='/data4/quiet/auxilliary/masks'
        maskf='patch_%s_Qband_mask_3sigma_ptsrc_n%s_v1.fits' % (patch,nside)
        mask=os.path.join(maskdir,maskf)
    else:
        mask=''

    if stage=='lcutsims':
        work_dir_root=work_dir

    for split in pairs.keys():
        [seta,setb]=pairs[split]
        if stage=='finalmap':
            mpirunner='mpirun'
            exe=EXE_FINALMAP
            stagefunc='finalmap'
            run=version
            job='%s_pol_%s_%s_%s_finalmap' % (patch,split,who,run)
            clock='4:00:00'
            memory='2000M'
            ntaskspernode=8
            nodes=4
            do_omp='#'; omp_num_threads=8
            constraints='SBATCH --constraint=ib'
            prefix='%s_%s_%s' % (run,patch,split)
            files='cmb_patch_%(patch)s_set%(seta)s_sub001_P_eqn.unf cmb_patch_%(patch)s_set%(setb)s_sub001_P_eqn.unf' % {'patch':patch,'seta':seta,'setb':setb}
            mask=mask
            args=prefix

        if stage=='finalmapmask':
            mpirunner='mpirun'
            exe=EXE_FINALMAP
            stagefunc='finalmap'
            run=version
            job='%s_pol_%s_%s_%s_finalmap_mask' % (patch,split,who,run)
            clock='4:00:00'
            memory='2000M'
            ntaskspernode=8
            nodes=4
            do_omp='#'; omp_num_threads=8
            constraints='SBATCH --constraint=ib'
            prefix='%s_%s_%s_mask' % (run,patch,split)
            files='cmb_patch_%(patch)s_set%(seta)s_sub001_P_eqn.unf cmb_patch_%(patch)s_set%(setb)s_sub001_P_eqn.unf' % {'patch':patch,'seta':seta,'setb':setb}
            mask=mask
            args=prefix

        elif stage=='lcut':
            mpirunner='mpirun'
            exe=EXE_LCUT
            stagefunc='lcut'
            run=version
            job='%s_pol_%s_%s_%s_lcut' % (patch,split,who,run)
            clock='1:00:00'
            memory='24000M'
            ntaskspernode=1
            nodes=1
            do_omp=''; omp_num_threads=1
            constraints='SBATCH --constraint=hugemem --comment=resv:hugemem.0.0'
            prefix='%s_%s_%s' % (run,patch,split)
            files='%(prefix)s_map2mask.fits %(prefix)s_inv_N_diffmap.unf %(prefix)s_true_diffmap.fits' % {'prefix':prefix}
            outprefix='%s_%s' % (patch,split)
            ell=DEFAULT_ELL
            args='%s %i' % (outprefix,ell)

        elif stage=='lcutmask':
            mpirunner='mpirun'
            #onlymap=''
            if quick:
                onlymap='onlymap'
            else:
                onlymap=''
            exe=EXE_LCUT
            stagefunc='lcut'
            run=version
            job='%s_pol_%s_%s_%s_lcut_mask' % (patch,split,who,run)
            clock='1:00:00'
            memory='24000M'
            ntaskspernode=1
            nodes=1
            do_omp=''; omp_num_threads=1
            constraints='SBATCH --constraint=hugemem --comment=resv:hugemem.0.0'
            prefix='%s_%s_%s_mask' % (run,patch,split)
            files='%(prefix)s_map2mask.fits %(prefix)s_inv_N_diffmap.unf %(prefix)s_true_diffmap.fits' % {'prefix':prefix}
            outprefix='%s_%s_mask' % (patch,split)
            ell=DEFAULT_ELL
            args='%s %i %s' % (outprefix,ell,onlymap)

        elif stage=='lcutsum':
            mpirunner='mpirun'
            exe=EXE_LCUT
            stagefunc='lcut'
            run=version
            job='%s_pol_%s_%s_%s_lcutsum' % (patch,split,who,run)
            clock='1:00:00'
            memory='24000M'
            ntaskspernode=1
            nodes=1
            do_omp=''; omp_num_threads=1
            constraints='SBATCH --constraint=hugemem --comment=resv:hugemem.0.0'
            prefix='%s_%s_%s' % (run,patch,split)
            files='%(prefix)s_map2mask.fits %(prefix)s_inv_N.unf %(prefix)s_map.fits' % {'prefix':prefix}
            outprefix='%s_%s_sum' % (patch,split)
            ell=DEFAULT_ELL
            args='%s %i' % (outprefix,ell)

        elif stage=='lcutsims':
            # This uses a slight hijack of terminology..
            mpirunner='mpifor'
            onlymap='onlymap'
            exe=EXE_LCUT
            stagefunc="'%s lcut" % exe
            exe="''" # Now ditch this
            run=version
            job='%s_pol_%s_%s_%s_lcut_sims' % (patch,split,who,run)
            clock='1:00:00'
            memory='2000M'
            ntaskspernode=8
            nodes=1
            do_omp='#'; omp_num_threads=1
            constraints='SBATCH --constraint=ib'
            prefix='%s_%s_%s_mask' % (run,patch,split)
            files='../%(prefix)s_map2mask.fits ../%(prefix)s_inv_N_diffmap.unf $i ${i/_true_diffmap.fits/}' % {'prefix':prefix}
            ell=DEFAULT_ELL
            args="%i %s'" % (ell,onlymap)
            mask="'%(prefix)s_sim???_true_diffmap.fits'" % {'prefix':prefix}
            simsdir='%(prefix)s_sims' % {'prefix':prefix}
            work_dir=os.path.join(work_dir_root,simsdir)


        substitutions={'job':job,
                       'account':'quiet',
                       'walllimit':clock,
                       'memory':memory,
                       'ntaskspernode':str(ntaskspernode),
                       'nodes':str(nodes),
                       'constraints':constraints,
                       'users':user,
                       'tmp':'2G',
                       'do_omp':do_omp,
                       'omp_num_threads':str(omp_num_threads),
                       'stagefunc':stagefunc,
                       'work_dir':work_dir,
                       'exe':exe,
                       'files':files,
                       'mask':mask,
                       'args':args,
                       'mpirunner':mpirunner}

        substituted = template % substitutions
        if verbose:
            print "'ere we go:\n"
            print substituted

        outf='batch-%s-%s.txt' % (stage,split)
        #outf=os.path.abspath(outf)
        if not os.path.exists(outdir): os.mkdir(outdir)
        outf=os.path.join(outdir,outf)
        if os.path.exists(outf): os.rename(outf,'%s.old'%outf)
        outfile=open(outf,'w')
        outfile.write(substituted)
        print "batch for (%s's %s) written to %s" % (split,stage,outf)
        outfile.close()

    print 'Now do:'
    print 'cd %s' %outdir
    print 'for i in batch-%s-*.txt; do sbatch $i; echo $i; done' %(stage)

    sys.exit(0)
