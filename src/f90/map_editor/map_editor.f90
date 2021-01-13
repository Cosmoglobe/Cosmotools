program map_editor
  use healpix_types
  use map_editor_utils
  use map_editor_simple_ops_mod
  use map_editor_complex_ops_mod
  implicit none

  ! *******************************************************
  ! *                                                     *
  ! *          Utility for editing Healpix maps           *
  ! *                                                     *
  ! *  written by Hans Kristian Eriksen, November, 2004   *
  ! *                                                     *
  ! *         Copyright 2004. All rights reserved.        *
  ! *                                                     *
  ! *******************************************************

  !
  ! ChangeLog:
  ! 
  ! November 2, 2004   -- First version
  !
  ! August 27-31, 2018 -- Documented by Daniel Herman

  integer(i4b)       :: iargc, lmin, lmax, lcut, nside_out, unit, i, l, m, n, seed
  integer(i4b)       :: nside, ordering, nmaps, component, s_max, ncol, nsim, verbose
  character(len=3)   :: suffix
  character(len=256) :: mapname_in1, mapname_in2, infofile, mapname_out, maskfile, rmsfile
  character(len=256) :: beamfile_in, beamfile_out, option, bandpass_in, output_file, simfile_name
  character(len=256) :: string_real, string_int, operation, beaminfo, covartype, map2mask_file, outprefix
  real(dp)           :: value, sigma_0, rms, cl0, r_fill
  real(dp)           :: fwhm_in, fwhm_out, md(4), fact, f(3)

  real(dp),     pointer, dimension(:,:) :: map, map2, resmap
  logical(lgt), pointer, dimension(:)   :: mask
  character(len=80), dimension(180)  :: header

  if (iargc() == 0) then
     write(*,*) ' '
     write(*,*) 'Usage: map_editor [operation] [arg1] [arg2] ...'
     write(*,*) ' '
     write(*,*) '       ____ SIMPLE OPERATIONS  ____'
     write(*,*) '       ---- SINGLE MAP OPERATIONS ----'
     write(*,*) '       scale, add_offset, log, ln, exp, abs, inv, sqrt, '
     write(*,*) '       asinh, hitcount2rms, max_scalar, min_scalar, missing2mask'
     write(*,*) '       QU2P, rms2mask, amp2mask, hitcount2mask, invert_mask,'
     write(*,*) '       missing2val'
     write(*,*) ''
     write(*,*) '       ---- TWO MAP OPERATIONS  ----'
     write(*,*) '       add, subtract, multiply, divide, half_sum, half_diff, '
     write(*,*) '       max, min, splice'
     write(*,*) ''
     write(*,*) '       _____ COMPLEX OPERATIONS _____'
     write(*,*) '       ---- COMMON OPERATIONS ----'
     write(*,*) '       smooth, ud_grade, subtract_mono_dipole, fit_gain_offset,'
     write(*,*) '       fit_gain_offset_dipole, mask2misspix, fix_monopole,'
     write(*,*) '       smooth_rms, smooth_rms_quick'
     write(*,*) ''
     write(*,*) '       ---- PRINT OPERATIONS ----'
     write(*,*) '       print_stats (stats), print_stats_col, print_map_to_ascii,'
     write(*,*) '       print_two_maps_to_ascii, print_scaled_gal_avg,'
     write(*,*) '       print_isolatitude, print_isolatitude_var, badcount,'
     write(*,*) '       compute_mean_stddev, size_info, zero_count, value_count,'
     write(*,*) '       misspix_count'
     write(*,*) ''
     write(*,*) '       ---- MASK OPERATIONS ----'
     write(*,*) '       expand_mask, maskcount, ud_grade_mask, apply_mask, '
     write(*,*) '       create_source_mask_from_file'
     write(*,*) ''
     write(*,*) '       ---- ADDITIONAL OPERATIONS ----'
     write(*,*) '       fit_line_to_ASCII_data, fit_ideal_dust, compute_spectral_index_map,'
     write(*,*) '       convert_beam, output_pointsource, project2xy, crossproduct,'
     write(*,*) '       shift_columns, compute_weighted_sum, firas_calib'
     write(*,*) '       partrans, make_ptsrc_map, copy_T_to_QU,copy_pol_to_T'
     write(*,*) '       add_gaussian_noise, add_gaussian_noise_sqrt_N, '
     write(*,*) '       subtract_mono_dipole_highl, extract_multipole_range,'
     write(*,*) '       summarize_detector_angles, read_covmatrix, make_co_region_map,'
     write(*,*) '       qu_transport_map, merge_maps, median_filter_source_holes, '
     write(*,*) '       median_filter_specific_value, median_filter_misspix,'
     write(*,*) '       copy_T_to_QU, copy_pol_to_T, copy_pol_to_pol,'
     write(*,*) '       fill_zero_mask_udgrade, rescale_dust_amplitude_map,'
     write(*,*) '       generate_beam'
     write(*,*) ''
     write(*,*) '   For more information on any given operation, execute the '
     write(*,*) '   program with options "help + operation", e.g.:'
     write(*,*) ''
     write(*,*) '       map_editor help smooth'
     write(*,*) ''
     stop
  end if

  unit = 25
  call getarg(1,operation)
  n = len(trim(operation))
  suffix = operation(n-2:n)

  if (trim(operation) == 'scale' .or. trim(operation) == 'add_offset' .or. &
       & trim(operation) == 'log' .or. trim(operation) == 'ln' .or. &
       & trim(operation) == 'exp' .or. trim(operation) == 'sqrt' .or. &
       & trim(operation) == 'abs' .or. trim(operation) == 'inv' .or. &
       & trim(operation) == 'hitcount2rms' .or. trim(operation) == 'rms2mask' .or. &
       & trim(operation) == 'amp2mask' .or. trim(operation) == 'asinh' .or. &
       & trim(operation) == 't2a' .or. trim(operation) == 'a2t' .or. &
       & trim(operation) == 'max_scalar' .or. trim(operation) == 'min_scalar' .or. &
       & trim(operation) == 'hitcount2mask' .or. trim(operation) == 'QU2P' .or. &
       & trim(operation) == 'missing2mask' .or. trim(operation) == 'missing2val' .or. &
       & trim(operation) == 'invert_mask') then
     
     if (iargc() < 2) then
        write(*,*) ''
        write(*,*) '   The following operations act on a single map, and require '
        write(*,*) '   only an input filename and an output filename:'
        write(*,*) ''
        write(*,*) '        log, ln, exp, abs '
        write(*,*) ''
        write(*,*) '   The following operations require an additional floating'
        write(*,*) '   point argument:'
        write(*,*) ''
        write(*,*) '        scale, add_offset, missing2val, rms2mask ..'
        write(*,*) ''
        write(*,*) '   Example usage:'
        write(*,*) ''
        write(*,*) '        map_editor log inmap.fits outmap.fits'
        write(*,*) '        map_editor scale inmap.fits outmap.fits 2.0'
        write(*,*) ''
        stop
     end if
        
     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     if (iargc() == 4) then
        call getarg(4,string_real)
        read(string_real,*) value
     else
        value = 0.
     end if

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)

     call operate_on_single_map(map, operation, value)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'ring2nest' .or. trim(operation) == 'nest2ring') then

     if (iargc() /= 3) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor ring2nest/nest2ring [input map] [output map]'
        write(*,*) ''
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)

     if (trim(operation) == 'ring2nest') then

        if (ordering == 2) then
           write(*,*) 'Error: Map is already in NESTED format'
           stop
        end if

        do i = 1, nmaps
           call convert_ring2nest(nside, map(:,i))
        end do

        ordering = 2

     else

        if (ordering == 1) then
           write(*,*) 'Error: Map is already in RING format'
           stop
        end if

        do i = 1, nmaps
           call convert_nest2ring(nside, map(:,i))
        end do

        ordering = 1

     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')


  else if (trim(operation) == 'scale_TQU') then

     if (iargc() /= 6) then
        write(*,*) ''
        write(*,*) 'Usage: map_editor scale_TQU [infile] [outfile] [T scale] [Q scale] [U scale]'
        write(*,*) ''
        stop
     end if

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     call getarg(4,string_real)
     read(string_real,*) f(1)
     call getarg(5,string_real)
     read(string_real,*) f(2)
     call getarg(6,string_real)
     read(string_real,*) f(3)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     do i = 1, 3
        map(:,i) = map(:,i) * f(i)
     end do
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'apply_mask') then

     if (iargc() < 4) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor apply_mask [input map] [mask] [output map]'
        write(*,*) ''
        stop
     end if

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)
     if (iargc() == 5) then
        call getarg(5,string_real)
        read(string_real,*) fact
     end if
     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     if (iargc() == 4) then
        call apply_mask1(maskfile, nside, ordering, map)
     else if (iargc() == 5) then
        call apply_mask1(maskfile, nside, ordering, map, fact)
     else
        write(*,*) ''
        write(*,*) "Usage: map_editor apply_mask [input map] [mask] [output map]"
        write(*,*) ''
     end if
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'fix_monopole') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_in2)
     call getarg(4,maskfile)
     call getarg(5,mapname_out)
     if (iargc() == 5) then
        call fix_monopole(mapname_in1, mapname_in2, maskfile, mapname_out)
     else
        write(*,*) ''
        write(*,*) "Usage: map_editor fix_monopole [input map] [reference map] [mask] [output map]"
        write(*,*) ''
     end if

  else if (trim(operation) == 'compute_mean_stddev') then

     if (iargc() < 4) then
        write(*,*) ''
        write(*,*) '   This operation requires at least 3 arguments! See help for more details.'
        write(*,*) ''
        write(*,*) '   Usage: map_editor compute_mean_stddev [outprefix] [mask] [map1]'
        write(*,*) '                             [map2](optional) [map3](optional) ...'
        write(*,*) ''
        stop
     endif

     call getarg(2,outprefix)
     call getarg(3,maskfile)

     call output_mean_and_stddev(nside, ordering, nmaps, header, map, map2)
     call write_result_map(trim(outprefix)//'_mean.fits', nside, ordering, header, map, suffix=='_dp')
     call write_result_map(trim(outprefix)//'_stddev.fits', nside, ordering, header, map2, suffix=='_dp')

  else if (trim(operation) == 'add' .or. trim(operation) == 'subtract' .or. &
       & trim(operation) == 'multiply' .or. trim(operation) == 'divide' .or. &
       & trim(operation) == 'max' .or. trim(operation) == 'min' .or. &
       & trim(operation) == 'half_sum' .or. trim(operation) == 'half_diff' .or. &
       & trim(operation) == 'splice') then

     if (iargc() < 4) then
       write(*,*) ''
       write(*,*) '   The following operations act on two maps, and require '
       write(*,*) '   two input filenames and one output filename:'
       write(*,*) ''
       write(*,*) '        add, subtract, multiply, divide, half_sum, half_diff '
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor subtract inmap1.fits inmap2.fits outmap.fits'
       write(*,*) ''
       stop
     endif
       
     ! Simple two map operations
     call getarg(2,mapname_in1)
     call getarg(3,mapname_in2)
     call getarg(4,mapname_out)

     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     call initialize_second_map(mapname_in2, nside, ordering, nmaps, map2)

     call operate_on_two_maps(map, map2, operation)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'weighted_sum') then

     if (iargc() /= 3) then
         write(*,*) ''
         write(*,*) '   Usage: map_editor weighted_sum [infofile] [output map]'
         write(*,*) ''
         stop
     endif

     call getarg(2,infofile)
     call getarg(3,mapname_out)
     call compute_weighted_sum(infofile, nside, ordering, nmaps, resmap, header)
     call write_result_map(mapname_out, nside, ordering, header, resmap, suffix=='_dp')

  else if (trim(operation) == 'ptsrc') then

     call output_pointsource

  else if (trim(operation) == 'smooth') then

     if (iargc() < 9) then
        write(*,*) ''
        write(*,*) '   Usage:  map_editor smooth [beam id] [input filename] [lmin] [lmax]'
        write(*,*) '               [nside_out] [input beam] [output beam]'
        write(*,*) '               [output filename] [radius fill in pixels; optional]'
        write(*,*) ''
        stop
     end if

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside

     if (iargc() == 10) then
        call getarg(10,string_int)
        read(string_int,*) r_fill
     else
        r_fill = -1.d0
     end if

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'f2f_EB') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g_EB') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f_EB') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g_EB') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else 
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'smooth_zerospin') then

     if (iargc() < 9) then
        write(*,*) ''
        write(*,*) '   Usage:  map_editor smooth_zerospin [beam id] [input filename] [lmin] [lmax]'
        write(*,*) '               [nside_out] [input beam] [output beam]'
        write(*,*) '               [output filename] [radius fill in pixels; optional]'
        write(*,*) ''
        stop
     end if

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside

     if (iargc() == 10) then
        call getarg(10,string_int)
        read(string_int,*) r_fill
     else
        r_fill = -1.d0
     end if

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .false., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'f2f_EB') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g_EB') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f_EB') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g_EB') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_map_zerospin(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, .true., &
             & fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else 
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'smooth_rms_quick') then
     if (iargc() < 10) then
        write(*,*) ''
        write(*,*) '   Usage:  map_editor smooth_rms_quick [beam id] [input filename] [lmin] [lmax]'
        write(*,*) '               [nside_out] [input beam] [output beam]'
        write(*,*) '               [output filename] [seed] [radius fill in pixels; optional]'
        write(*,*) ''
        stop
     end if

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside
     call getarg(10,string_int)
     read(string_int,*) seed

     simfile_name='none'

     if (iargc() == 11) then
        call getarg(11,string_int)
        read(string_int,*) r_fill
     else
        r_fill = -1.d0
     end if

     nsim = 5 ! only use 5 sims for the fast estimation 
     verbose = 0

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'smooth_rms') then
     if (iargc() < 11) then
        write(*,*) ''
        write(*,*) '   Usage:  map_editor smooth_rms [beam id] [input filename] [lmin] [lmax]'
        write(*,*) '               [nside_out] [input beam] [output beam] [output filename]'
        write(*,*) '               [seed] [N_sims] [OPTIONS]'
        write(*,*) ''        
        write(*,*) '   Options: -verbose <integer>   : add more output to terminal (1: per sim, 2:debug) '        
        write(*,*) '            -rfill <integer>     : radius of pixels to fill in bad pixels'        
        write(*,*) '            -sim_rms <string>    : outputs the simulated RMS used for scaling to file given by input'        
        stop
     end if

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside
     call getarg(10,string_int)
     read(string_int,*) seed
     call getarg(11,string_int)
     read(string_int,*) nsim

     verbose = 0
     r_fill = -1.d0

     simfile_name='none'

     i = 12
     do while (i <= iargc())
        call getarg(i,string_int)
        if (trim(string_int) == '-verbose') then
           call getarg(i+1,string_int)
           read(string_int,*) verbose
           i = i+1
        else if (trim(string_int) == '-rfill') then
           call getarg(i+1,string_int)
           read(string_int,*) r_fill
           i = i+1
        else if (trim(string_int) == '-sim_rms') then
           call getarg(i+1,simfile_name)
           i = i+1
        else
           write(*,*) 'Unknown option', trim(string_int)
           stop
        end if
        i = i+1
     end do

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'smooth_rms_degrade') then
     if (iargc() < 11) then
        write(*,*) ''
        write(*,*) '   Usage:  map_editor smooth_rms_degrade [beam id] [input filename] [lmin]'
        write(*,*) '           [lmax] [nside_out] [input beam] [output beam] [output filename]'
        write(*,*) '           [seed] [N_sims] [OPTIONS]'
        write(*,*) ''        
        write(*,*) '   Options: -verbose <integer>   : add more output to terminal (1: per sim, 2:debug) '        
        write(*,*) '            -rfill <integer>     : radius of pixels to fill in bad pixels'        
        write(*,*) '            -sim_rms <string>    : outputs the simulated RMS used for scaling to file given by input'        
        stop
     end if

     call getarg(2,beaminfo)
     call getarg(3,mapname_in1)
     call getarg(4,string_int)
     read(string_int,*) lmin
     call getarg(5,string_int)
     read(string_int,*) lmax
     call getarg(6,string_int)
     read(string_int,*) nside
     call getarg(10,string_int)
     read(string_int,*) seed
     call getarg(11,string_int)
     read(string_int,*) nsim

     verbose = 0
     r_fill = -1.d0

     simfile_name='none'

     i = 12
     do while (i <= iargc())
        call getarg(i,string_int)
        if (trim(string_int) == '-verbose') then
           call getarg(i+1,string_int)
           read(string_int,*) verbose
           i = i+1
        else if (trim(string_int) == '-rfill') then
           call getarg(i+1,string_int)
           read(string_int,*) r_fill
           i = i+1
        else if (trim(string_int) == '-sim_rms') then
           call getarg(i+1,simfile_name)
           i = i+1
        else
           write(*,*) 'Unknown option', trim(string_int)
           stop
        end if
        i = i+1
     end do

     if (trim(beaminfo) == 'f2f') then
        ! Input beam from file, output beam from file
        call getarg(7,beamfile_in)
        call getarg(8,beamfile_out)
        call smooth_rms_true_degrade(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'f2g') then
        ! Input beam from file, output beam from gaussian
        call getarg(7,beamfile_in)
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true_degrade(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, beamfile_in=beamfile_in, fwhm_out=fwhm_out)
     else if (trim(beaminfo) == 'g2f') then
        ! Input beam from gaussian, output beam from file
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,beamfile_out)
        call smooth_rms_true_degrade(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, beamfile_out=beamfile_out)
     else if (trim(beaminfo) == 'g2g') then
        ! Input beam from gaussian, output beam from gaussian
        call getarg(7,string_real)
        read(string_real,*) fwhm_in
        call getarg(8,string_real)
        read(string_real,*) fwhm_out
        call smooth_rms_true_degrade(mapname_in1, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
             & seed, nsim, verbose, simfile_name, fwhm_in=fwhm_in, fwhm_out=fwhm_out)
     else
        write(*,*) 'Invalid beam option. Exiting.'
        stop
     end if

     call getarg(9,mapname_out)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')
     
  else if (trim(operation) == 'shift_columns') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)
     call getarg(4,string_int)
     read(string_int,*) ncol

     call shift_columns(mapname_in1, ncol, nside, ordering, header, map)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'add_gaussian_noise') then

     if (iargc() /= 6) then
        write(*,*) ''
        write(*,*) '   See map_editor help add_gaussian_noise for usage help'
        write(*,*) ''
        stop
     endif
     
     call getarg(2,mapname_in1)
     call getarg(3,covartype)
     call getarg(4,rmsfile)
     call getarg(5,string_int)
     read(string_int,*) seed
     call getarg(6,mapname_out)
     
     if (trim(covartype) == 'sqrt_N') then

        call getarg(7,map2mask_file)
        
        ! Assume that second file contains sqrt_N of Gaussian noise
        call add_gaussian_noise_sqrt_N(mapname_in1, rmsfile, map2mask_file, seed, nside, ordering, nmaps, map, header)

     else

        if (iargc() == 6) then

           ! Assume that second file contains RMS of Gaussian noise
           call add_gaussian_noise(mapname_in1, rmsfile, seed, nside, ordering, nmaps, map, header)

        else if (iargc() == 7) then
           ! Assume that second file contains Nobs of Gaussian noise
           call getarg(7,string_real)
           read(string_real,*) sigma_0
           call add_gaussian_noise(mapname_in1, rmsfile, seed, nside, ordering, nmaps, map, header, sigma_0)
        end if

     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'subtract_mono_dipole') then

     if (iargc() /= 4 .and. iargc() /= 8) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor subtract_mono_dipole [input map] [mask] [output map]'
        write(*,*) '                               (optional) [mono] [dipx] [dipy] [dipz]'
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)
     if (iargc() == 8) then
        call getarg(5,string_real)
        read(string_real,*) md(1)
        call getarg(6,string_real)
        read(string_real,*) md(2)
        call getarg(7,string_real)
        read(string_real,*) md(3)
        call getarg(8,string_real)
        read(string_real,*) md(4)
        call subtract_mono_dipole(mapname_in1, maskfile, nside, ordering, map, header, md)
     else
        call subtract_mono_dipole(mapname_in1, maskfile, nside, ordering, map, header)
     end if

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'subtract_mono_dipole_highl') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,string_int)
     read(string_int,*) lmax
     call getarg(5,string_int)
     read(string_int,*) lcut
     call getarg(6,mapname_out)
     call subtract_mono_dipole_highl(mapname_in1, maskfile, lmax, lcut, nside, ordering, map, header)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'partrans') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) nside_out
     call getarg(4,mapname_out)
     call initialize_single_map(mapname_in1, nside, ordering, nmaps, header, map)
     write(*,*) nside, ordering, nside_out
     call qu_transport_map(nside_out, map)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'make_co_region_map') then

     call getarg(2,mapname_in1)
     call getarg(3,mapname_out)

     call make_co_region_map(mapname_in1, nside, ordering, map)

     call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=size(map,2)==3)
     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'extract_multipole_range') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) lmin
     call getarg(4,string_int)
     read(string_int,*) lmax
     call getarg(5,mapname_out)

     call extract_multipole_range(mapname_in1, lmin, lmax, nside, ordering, map, header)

     call write_result_map(mapname_out, nside, ordering, header, map, suffix=='_dp')

  else if (trim(operation) == 'ud_grade') then

     if (iargc() /= 4) then
        write(*,*) 'Usage: map_editor ud_grade [input map] [nside_out] [output map]'
        stop
     end if

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) nside
     call getarg(4,mapname_out)
     call ud_grade_map_editor(mapname_in1, nside, mapname_out, suffix=='_dp')

  else if (trim(operation) == 'ud_grade_rms') then

     if (iargc() /= 4) then
        write(*,*) 'Usage: map_editor ud_grade_rms [input map] [nside_out] [output map]'
        stop
     end if

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) nside
     call getarg(4,mapname_out)
     call ud_grade_rms_map_editor(mapname_in1, nside, mapname_out, suffix=='_dp')

  else if (trim(operation) == 'median_filter_source_holes') then

     call median_filter_holes

  else if (trim(operation) == 'median_filter_specific_value') then

     call median_filter_specific_value

  else if (trim(operation) == 'median_filter_misspix' .or. trim(operation) == 'fillpix') then

     call median_filter_misspix

  else if (trim(operation) == 'median_filter_misspix_neighbor') then

     call median_filter_misspix_neighbor

  else if (trim(operation) == 'mask2misspix') then          
                         
     call mask2misspix 

  else if (trim(operation) == 'expand_mask') then

     call expand_mask

  else if (trim(operation) == 'print_map_to_ascii' .or. trim(operation) == 'map2ascii') then

     if (iargc() /= 4) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor print_map_to_ascii [mapfile in] [mask] [mapfile out]'
        write(*,*) ''
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_map_to_ascii(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'print_two_maps_to_ascii') then

     if (iargc() /= 5) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor print_two_maps_to_ascii [mapfile1] [mapfile2]'
        write(*,*) '                                              [mask] [mapfile out]'
        write(*,*) ''
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,mapname_in2)
     call getarg(4,maskfile)
     call getarg(5,mapname_out)

     call print_two_maps_to_ascii(mapname_in1, mapname_in2, maskfile, mapname_out)

  else if (trim(operation) == 'print_isolatitude') then

     if (iargc() /= 4) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor print_isolatitude [mapfile] [mask] [outfile]'
        write(*,*) ''
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_isolatitude(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'print_isolatitude_var') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)
     call getarg(4,mapname_out)

     call print_isolatitude_var(mapname_in1, maskfile, mapname_out)

  else if (trim(operation) == 'fit_ideal_dust') then

     call fit_ideal_dust

  else if (trim(operation) == 'help') then

     call getarg(2,option)        
     call print_help(option)

  else if (trim(operation) == 'print_maximum') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) component
     call getarg(4,string_int)
     read(string_int,*) fwhm_in

     call print_maximum(mapname_in1, component, fwhm_in)

  else if (trim(operation) == 'print_stats' .or. trim(operation) == 'stats') then
     
     if (iargc() == 2) then
        call getarg(2,mapname_in1)
        maskfile = 'dummy'
        call print_stats(mapname_in1,maskfile)
        stop
     endif
     
     if (iargc() > 3) then
        write(*,*) 'Usage: map_editor print_stats [input map] [mask]'
        stop
     endif

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)

     call print_stats(mapname_in1, maskfile)

  else if (trim(operation) == 'print_stats_col') then

     call getarg(2,mapname_in1)
     call getarg(3,maskfile)

     call print_stats_col(mapname_in1, maskfile)

  else if (trim(operation) == 'summarize_detector_angles') then

     call getarg(2,mapname_in1)
     call getarg(3,string_int)
     read(string_int,*) s_max
     call getarg(4,string_int)
     read(string_int,*) nside
     call getarg(5,mapname_out)

     call summarize_detector_angles(mapname_in1, s_max, nside, mapname_out)

  else if (trim(operation) == 'compute_spectral_index_map') then

     call compute_index_map

  else if (trim(operation) == 'maskcount') then
     
     if (iargc() /= 2) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor maskcount [mask]'
        write(*,*) ''
        stop
     endif

     call getarg(2,mapname_in1)
     call maskcount(mapname_in1)

  else if (trim(operation) == 'badcount') then

     if (iargc() /= 2) then
        write(*,*) ''
        write(*,*) 'Usage: map_editor badcount [input map]'
        write(*,*) ''
        stop
     end if
     call getarg(2,mapname_in1)
     call badcount(mapname_in1)

  else if (trim(operation) == 'merge_maps') then

     call merge_maps(suffix=='_dp')

  else if (trim(operation) == 'make_ptsrc_map') then

     call make_ptsrc_map

  else if (trim(operation) == 'make_procmask') then

     call make_procmask

  else if (trim(operation) == 'convolve_masks') then

     call convolve_masks

  else if (trim(operation) == 'expand_mask_neighbours') then

     call expand_mask_neighbours

  else if (trim(operation) == 'print_scaled_gal_avg') then

     call print_scaled_gal_avg
  
  else if (trim(operation) == 'fit_gain_offset') then

     call fit_gain_offset

  else if (trim(operation) == 'fit_gain_offset_dipole') then

     call fit_gain_offset_dipole

  else if (trim(operation) == 'fit_gain') then

     call fit_gain

  else if (trim(operation) == 'fit_line_to_ASCII_data') then
     
     call fit_line
  
  else if (trim(operation) == 'convert_beam') then

     if (iargc() < 5) then
        write(*,*) ''
        write(*,*) '  Usage: map_editor convert_beam [fwhm] [lmax] [nmaps] [outfile] [infile (opt)]'
        write(*,*) ''
        stop
     end if

     call getarg(2,string_real)
     read(string_real,*) fwhm_in
     call getarg(3,string_int)
     read(string_int,*) lmax
     call getarg(4,string_int)
     read(string_int,*) nmaps
     call getarg(5,beamfile_out)

     if (iargc() == 6) then
        call getarg(6,beamfile_in)
        call convert_beam(fwhm_in, lmax, nmaps, beamfile_out, beamfile_in)
     else
        call convert_beam(fwhm_in, lmax, nmaps, beamfile_out)
     end if

  else if (trim(operation) == 'generate_beam') then

     if (iargc() < 5) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor generate_beam [fwhm] [lmax] [nmaps] [outfile]'
        write(*,*) ''
        stop
     end if

     call getarg(2,string_real)
     read(string_real,*) fwhm_in
     call getarg(3,string_int)
     read(string_int,*) lmax
     call getarg(4,string_int)
     read(string_int,*) nmaps
     call getarg(5,beamfile_out)
     call convert_beam(fwhm_in, lmax, nmaps, beamfile_out)

  else if (trim(operation) == 'create_source_mask') then

     call create_source_mask_from_file

  else if (trim(operation) == 'ud_grade_mask') then

     call ud_grade_mask

  else if (trim(operation) == 'ud_grade_rms') then

     call ud_grade_rms

  else if (trim(operation) == 'fill_zero_mask_udgrade') then

     call fill_zero_mask_udgrade

  else if (trim(operation) == 'size_info') then

     call map_size_info

  else if (trim(operation) == 'read_header') then

     call map_read_header

  else if (trim(operation) == 'copy_T_to_QU') then

     call copy_T_to_QU

  else if (trim(operation) == 'copy_pol_to_T') then

     call copy_pol_to_T

  else if (trim(operation) == 'rescale_dust_amplitude_map') then

     call rescale_dust_amplitude_map

  else if (trim(operation) == 'zero_count') then

     call zero_count

  else if (trim(operation) == 'value_count') then

     call value_count

  else if (trim(operation) == 'misspix_count') then

     call misspix_count

  else if (trim(operation) == 'polarization_fraction') then

     call polarization_fraction

  else if (trim(operation) == 'rms_half_ring_diff') then
     call calculate_rms_half_ring_diff

  else if (trim(operation) == 'fit_CO_lineratio') then
     call fit_CO_lineratio

  else if (trim(operation) == 'firas_calib') then
     
     if (iargc() /= 4) then
        write(*,*) ''
        write(*,*) '   Usage: map_editor firas_calib [bandpass file] [map weights (output)]'
        write(*,*) '                                 [output map]' ![units (MJy or uK)] '
        write(*,*) ''
        stop
     end if

     call getarg(2,bandpass_in)
     call getarg(3,output_file)
     call getarg(4,mapname_out)
     ! call getarg(5,option)

     call firas_calib(bandpass_in,output_file,mapname_out) !,option)   

  else
     write(*,*) 'Unknown operation. Exiting.'
     stop 
  end if


  ! Clean up arrays and exit 
  ! Compaq's F90 compiler doesn't seem to allow allocated() calls on pointers
!  if (allocated(map))    deallocate(map)
!  if (allocated(map2))   deallocate(map2)
!  if (allocated(resmap)) deallocate(resmap)
!  if (allocated(mask))   deallocate(mask)


contains

  subroutine print_help(option)
    implicit none

    character(len=*), intent(in) :: option

    if (trim(option) == 'scale' .or. trim(option) == 'add_offet' .or. &
         & trim(option) == 'log' .or. trim(option) == 'ln' .or. &
         & trim(option) == 'exp' .or. &
         & trim(option) == 'abs') then
       
       write(*,*) ''
       write(*,*) '   The following operations act on a single map, and require '
       write(*,*) '   only an input filename and an output filename:'
       write(*,*) ''
       write(*,*) '        log, ln, exp, abs '
       write(*,*) ''
       write(*,*) '   The following operations require an additional floating'
       write(*,*) '   point argument:'
       write(*,*) ''
       write(*,*) '        scale, add_offset'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor log inmap.fits outmap.fits'
       write(*,*) '        map_editor scale inmap.fits outmap.fits 2.0'
       write(*,*) ''

  else if (trim(option) == 'add' .or. trim(option) == 'subtract' .or. &
       & trim(option) == 'multiply' .or. trim(option) == 'divide' .or. &
       & trim(option) == 'half_sum' .or. trim(option) == 'half_diff') then
    
       write(*,*) ''
       write(*,*) '   The following operations act on two maps, and require '
       write(*,*) '   two input filenames and one output filename:'
       write(*,*) ''
       write(*,*) '        add, subtract, multiply, divide, half_sum, half_diff '
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor subtract inmap1.fits inmap2.fits outmap.fits'
       write(*,*) ''

    else if (trim(option) == 'ring2nest' .or. trim(option) == 'nest2ring') then

       write(*,*) ''
       write(*,*) '   The ring2nest and nest2ring operations change the pixel '
       write(*,*) '   numbering scheme of the input map from RING or NESTED to the '
       write(*,*) '   NESTED or RING respectively.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor ring2nest [input map] [output map]'
       write(*,*) ''

    else if (trim(option) == 'scale_TQU') then

       write(*,*) ''
       write(*,*) '   The scale_TQU operation scales T,Q, and U for an input map by a factor'
       write(*,*) '   according to the three inputs.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor scale_TQU [infile] [outfile] [T scale] [Q scale] [U scale]'
       write(*,*) ''

    else if (trim(option) == 'ud_grade') then

       write(*,*) ''
       write(*,*) '   The ud_grade operation returns a map with a new nside resolution.'
       write(*,*) '   This operation should not be used on masks as it smooths pixel values,'
       write(*,*) '   which can cause the new mask to have pixel values not equal to 0 or 1.'
       write(*,*) '   Instead use ud_grade_mask to change the nside of a mask.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor ud_grade [input map] [nside_out] [output map]'
       write(*,*) ''

    else if (trim(option) == 'ud_grade_mask') then

       write(*,*) ''
       write(*,*) '   The ud_grade_mask operation returns a mask with a new nside resolution.'
       write(*,*) '   The threshold argument is the value at which all downgraded pixels with'
       write(*,*) '   a value less than [threshold] are set to zero.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor ud_grade_mask [input map] [nside_out] [output map]'
       write(*,*) '                                   [threshold]'
       write(*,*) ''


    else if (trim(option) == 'ud_grade_rms') then

       write(*,*) ''
       write(*,*) '   The ud_grade_rms operation returns a rms map  with a new nside resolution.'
       write(*,*) '   The map is ud_graded in variance (rms^2) and the variance is scaled with'
       write(*,*) '   respect to the difference in npix.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor ud_grade_rms [input map] [nside_out] [output map]'
       write(*,*) ''


    else if (trim(option) == 'mask2misspix') then

       write(*,*) ''
       write(*,*) '   The mask2misspix operation returns a map with pixel values equal to'
       write(*,*) '   the missing pixel value where to input mask is less than 0.5.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor mask2misspix [input map] [mask] [output map]'
       write(*,*) ''


    else if (trim(option) == 'rms_half_ring_diff') then
       write(*,*) ''
       write(*,*) 'The rms_half_ring_diff operation rescales a (full) rms map from the wcov files'
       write(*,*) 'based on the half-ring difference of the same map and their corresponding'
       write(*,*) 'rms maps (also from wcov). The mask is to mask out the galactic plane.'
       write(*,*) 'The npol option is to specify how many polarizations to calculate the'
       write(*,*) 'scaling factor for. If less than nmpas, the last maps get the samescaling'
       write(*,*) 'as map1. Default is the same number as input nmaps.'
       write(*,*) ''
       write(*,*) "Usage: map_editor rms_half_ring_diff [half_ring1] [half_ring2] [rms1] [rms2]"
       write(*,*) "      [mask] [full rms map] [outmap] [npol (1/3), optional]"
       write(*,*) ''

       
    else if (trim(option) == 'apply_mask') then

       write(*,*) ''
       write(*,*) '   The apply_mask operation applies a mask to an input map, and creates'
       write(*,*) '   a new masked map. '
       write(*,*) '   Usage: map_editor apply_mask [input map] [mask] [output map]'
       write(*,*) ''

    else if (trim(option) == 'maskcount') then

       write(*,*) ''
       write(*,*) '   The maskcount operation returns the number of masked and unmasked'
       write(*,*) '   pixels in a mask file as well as the percentage of the map masked'
       write(*,*) '   and unmasked.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor maskcount [maskfile]'
       write(*,*) ''

    else if (trim(option) == 'make_procmask') then

       write(*,*) ''
       write(*,*) '   The operation takes in a map and a threshold, and returns a mask'
       write(*,*) '   where the pixels that have a value difference wrt. the average of'
       write(*,*) '   the neighbouring pixels that is greater than the threshold have'
       write(*,*) '   been masked out.'
       write(*,*) ''
       write(*,*) 'Usage: map_editor make_procmap [map] [threshold] [mask out] '
       write(*,*) ''

    else if (trim(option) == 'expand_mask_neighbours') then

       write(*,*) ''
       write(*,*) '   The operation takes in a mask and a threshold (optional), and'
       write(*,*) '   returns a copy of the input mask where also pixels that have '
       write(*,*) '   more than or equal to "threshold" neighbouring pixels that are'
       write(*,*) '   less than zero have been masked out. (default threshold = 3)'
       write(*,*) ''
       write(*,*) '   Usage: map_editor expand_mask_neighbours [mask in] [mask out] '
       write(*,*) '                     [n_neigh (optional)] '
       write(*,*) ''

    else if (trim(option) == 'weighted_sum') then

       write(*,*) ''
       write(*,*) '   The weighted_sum operation creates a new map that is the weighted sum of'
       write(*,*) '   a group of input maps. See help infofile for details on the input file.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor weighted_sum [infofile] [output map]'
       write(*,*) ''

    else if (trim(option) == 'infofile') then

       write(*,*) ''
       write(*,*) '   Layout of infofile for weighted_sum.'
       write(*,*) '   ------------------------------------'
       write(*,*) ''
       write(*,*) '   nside:'
       write(*,*) '   nmaps (TQU)'
       write(*,*) '   mapname1, weight1'
       write(*,*) '   mapname2, weight2'
       write(*,*) '   etc....'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) '   64'
       write(*,*) '   1'
       write(*,*) '   mymap1.fits 4'
       write(*,*) '   mymap2.fits 1'
       write(*,*) '   mymap3.fits 2'
       write(*,*) '   mymap4.fits 1'
       write(*,*) ''

    else if (trim(option) == 'compute_mean_stddev') then
       
       write(*,*) ''
       write(*,*) '   The compute_mean_stddev operation computes the mean and standard deviation'
       write(*,*) '   of a group of maps and then outputs two uniform maps: a mean map and a '
       write(*,*) '   standard deviation map. I.e. maps where all pixels are the mean and stddev'
       write(*,*) '   values named outprefix_mean.fits, outprefix_stddev.fits.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor compute_mean_stddev [outprefix] [mask] [map1]'
       write(*,*) '                             [map2](optional) [map3](optional) ...'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '   map_editor compute_mean_stddev 3maps mask.fits map1.fits map2.fits map3.fits'
       write(*,*) ''

    else if (trim(option) == 'subtract_mono_dipole') then
       
       write(*,*) ''
       write(*,*) '   The subtract_mono_dipole operation calcualtes and removes the monopole'
       write(*,*) '   and dipole signals from the input map, creating a residual map. The user'
       write(*,*) '   can input their own monopole and dipole values to subtract. The mono-'
       write(*,*) '   pole and dipole values are returned in console in the following form:'
       write(*,*) ''
       write(*,*) '   Monopole:    Dipole-x:    Dipole-y:    Dipole-z: '
       write(*,*) ''
       write(*,*) '   Usage: map_editor subtract_mono_dipole [input map] [mask] [output map]'
       write(*,*) '                               (optional) [mono] [dipx] [dipy] [dipz]'
       write(*,*) ''

    else if (trim(option) == 'smooth') then

       write(*,*) ''
       write(*,*) '   The smoothing operation reads a map, computes its spherical '
       write(*,*) '   harmonics transform, deconvolve a given beam (and pixel window),'
       write(*,*) '   convolve with a new given beam (and pixel window), and finally'
       write(*,*) '   outputs resulting the inverse spherical harmonics transform.'
       write(*,*) ''
       write(*,*) '   Both the input and output beams can be given on two formats,'
       write(*,*) '   either in the form of the FWHM of a Gaussian beam, or as a'
       write(*,*) '   FITS ascii table.'
       write(*,*) ''
       write(*,*) '   The beam id must be one of the following strings: '
       write(*,*) ''
       write(*,*) '        f2f -> Input and output beams are read from files'
       write(*,*) '        f2g -> Input beam is read from file, output beam is Gaussian'
       write(*,*) '        g2f -> Input beam is Gaussian, output beam is read from file'
       write(*,*) '        g2g -> Input and output beams are Gaussian'
       write(*,*) ''
       write(*,*) '   If a file is requested, then [input/output beam] must be a filename.'
       write(*,*) '   If a Gaussian is requested, then it must be a floating point number,'
       write(*,*) '   giving the FWHM of the Gaussian beam in arcmin.'
       write(*,*) ''
       write(*,*) '   The variables lmin, lmax and Nside may formally be choosen freely '
       write(*,*) '   (independent of the input map), but a good rule of thumb is'
       write(*,*) '   lmax < 3*Nside_in, and preferrably even lmax < 2*Nside_in.'
       write(*,*) ''
       write(*,*) '   The command line format is as follows:'
       write(*,*) ''
       write(*,*) '        map_editor smooth [beam id] [input filename] [lmin] [lmax]'
       write(*,*) '               [nside_out] [input beam] [output beam]'
       write(*,*) '               [output filename] [radius fill in pixels; optional]'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor smooth f2g map.fits 2 1024 512 beam.fits 60. smoothed.fits'
       write(*,*) ''

    else if (trim(option) == 'smooth_zerospin') then

       write(*,*) ''
       write(*,*) '   The smooth_zerospin operation smooths maps like the "smooth"' 
       write(*,*) '   operation, but maps in Q and U are smoothed like maps in T,'
       write(*,*) '   i.e. zero spin maps. For more help see the "smooth" operation'
       write(*,*) ''
       write(*,*) '   The command line format is as follows:'
       write(*,*) ''
       write(*,*) '        map_editor smooth_zerospin [beam id] [input filename] [lmin]'
       write(*,*) '               [lmax] [nside_out] [input beam] [output beam]'
       write(*,*) '               [output filename] [radius fill in pixels; optional]'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor smooth_zerospin g2g map.fits 2 1024 512 15. 60. smoothed.fits'
       write(*,*) ''

    else if (trim(option) == 'smooth_rms_quick') then
       write(*,*) ''
       write(*,*) '   The smoothing operation reads a RMS map and simultes 5 white'
       write(*,*) '   noise maps with the input RMS. Then it computes the spherical '
       write(*,*) '   harmonics transform, deconvolve a given beam (and pixel window),'
       write(*,*) '   convolve with a new given beam (and pixel window), and finally'
       write(*,*) '   outputs resulting the inverse spherical harmonics transform of'
       write(*,*) '   the white noise maps. Finally it computes the RMS of the smoothed'
       write(*,*) '   white noise maps and scales a smoothed map of the input to match'
       write(*,*) '   the simulated RMS amplitude.'
       write(*,*) '   '
       write(*,*) '   This is the most accurate (and correct) way of smoothing RMS maps.'
       write(*,*) '   To get a more accurate result one have to use many simulations,'
       write(*,*) "   see 'smooth_rms' or 'smooth_rms_degrade'"
       write(*,*) ''
       write(*,*) '   Both the input and output beams can be given on two formats,'
       write(*,*) '   either in the form of the FWHM of a Gaussian beam, or as a'
       write(*,*) '   FITS ascii table.'
       write(*,*) ''
       write(*,*) '   The beam id must be one of the following strings: '
       write(*,*) ''
       write(*,*) '        f2f -> Input and output beams are read from files'
       write(*,*) '        f2g -> Input beam is read from file, output beam is Gaussian'
       write(*,*) '        g2f -> Input beam is Gaussian, output beam is read from file'
       write(*,*) '        g2g -> Input and output beams are Gaussian'
       write(*,*) ''
       write(*,*) '   If a file is requested, then [input/output beam] must be a filename.'
       write(*,*) '   If a Gaussian is requested, then it must be a floating point number,'
       write(*,*) '   giving the FWHM of the Gaussian beam in arcmin.'
       write(*,*) ''
       write(*,*) '   The variables lmin, lmax, seed and Nside may formally be choosen freely '
       write(*,*) '   (independent of the input map), but a good rule of thumb is'
       write(*,*) '   lmax < 3*Nside_in, and preferrably even lmax < 2*Nside_in.'
       write(*,*) ''
       write(*,*) '   The command line format is as follows:'
       write(*,*) ''
       write(*,*) '        map_editor smooth_rms_quick [beam id] [input filename] [lmin] [lmax]'
       write(*,*) '               [nside_out] [input beam] [output beam] [output filename]'
       write(*,*) '               [seed] [r_fill; pixel radius for bad pixels, optional]'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor smooth_rms_quick g2g rms_map.fits 0 1024 512 5. 60. rms_smoothed.fits 46261 '
       write(*,*) ''

    else if (trim(option) == 'smooth_rms' .or. trim(operation) == 'smooth_rms_degrade' ) then

       write(*,*) ''
       write(*,*) '   The smoothing operation reads a RMS map and simultes N white'
       write(*,*) '   noise maps with the input RMS. Then it computes the spherical '
       write(*,*) '   harmonics transform, deconvolve a given beam (and pixel window),'
       write(*,*) '   convolve with a new given beam (and pixel window), and finally'
       write(*,*) '   outputs resulting the inverse spherical harmonics transform of'
       write(*,*) '   the white noise maps. Finally it computes the RMS of the smoothed'
       write(*,*) '   white noise maps and scales a smoothed map of the input to match'
       write(*,*) '   the simulated RMS amplitude. Note: the input RMS map is smoothed'
       write(*,*) '   in variance-domain (i.e. RMS^2)'
       write(*,*) '   '
       write(*,*) '   This is the most accurate (and correct) way of smoothing RMS maps.'
       write(*,*) '   To get a more accurate result one have to use MORE simulations.'
       write(*,*) '   '
       write(*,*) "   The 'smooth_rms_degrade' option will degrade the input RMS map to"
       write(*,*) '   the output Nside before computing the white noise maps. This is also'
       write(*,*) '   done in variance-domain.'
       write(*,*) ''
       write(*,*) '   Both the input and output beams can be given on two formats,'
       write(*,*) '   either in the form of the FWHM of a Gaussian beam, or as a'
       write(*,*) '   FITS ascii table.'
       write(*,*) ''
       write(*,*) '   The beam id must be one of the following strings: '
       write(*,*) ''
       write(*,*) '        f2f -> Input and output beams are read from files'
       write(*,*) '        f2g -> Input beam is read from file, output beam is Gaussian'
       write(*,*) '        g2f -> Input beam is Gaussian, output beam is read from file'
       write(*,*) '        g2g -> Input and output beams are Gaussian'
       write(*,*) ''
       write(*,*) '   If a file is requested, then [input/output beam] must be a filename.'
       write(*,*) '   If a Gaussian is requested, then it must be a floating point number,'
       write(*,*) '   giving the FWHM of the Gaussian beam in arcmin.'
       write(*,*) ''
       write(*,*) '   The variables lmin, lmax, seed and Nside may formally be choosen freely '
       write(*,*) '   (independent of the input map), but a good rule of thumb is'
       write(*,*) '   lmax < 3*Nside_in, and preferrably even lmax < 2*Nside_in.'
       write(*,*) ''
       write(*,*) '   The command line format is as follows:'
       write(*,*) ''
       write(*,*) '   Usage:  map_editor {smooth_rms, smooth_rms_degrade} [beam id] '
       write(*,*) '           [input filename] [lmin] [lmax] [nside_out] [input beam] '
       write(*,*) '           [output beam] [output filename] [seed] [N_sims] [OPTIONS]'
       write(*,*) ''        
       write(*,*) '   Options: -verbose <integer>   : add more output to terminal (1: per sim, 2:debug) '        
       write(*,*) '            -rfill <integer>     : radius of pixels to fill in bad pixels'        
       write(*,*) ''        
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor smooth_rms g2g rms_map.fits 0 1024 512 5. 60. '
       write(*,*) '             rms_smoothed.fits 46261 100'
       write(*,*) ''
       
    else if (trim(option) == 'add_gaussian_noise') then

       write(*,*) ''
       write(*,*) '   Gaussian noise can be added to an existing map by calling the program'
       write(*,*) '   with the "add_gaussian_noise" option. This mode of operation assumes'
       write(*,*) '   an input filename, an RMS/Nobs filename, a seed and an output filename.'
       write(*,*) '   If only four arguments are provided, the second file is assumed to '
       write(*,*) '   contain the noise RMS. If five is provided (the fifth being '
       write(*,*) '   the standard deviation of a single detector hit, sigma_0), it is'
       write(*,*) '   assumed to contain the number of observations for each pixel, Nobs.'
       write(*,*) ''
       write(*,*) '   Example usage:'
       write(*,*) ''
       write(*,*) '        map_editor add_gaussian_noise map.fits type rms.fits -128341 noisymap.fits '
       write(*,*) '        map_editor add_gaussian_noise map.fits type nobs.fits -128341 noisymap.fits 100.'
       write(*,*) ''

    else if (trim(option) == 'convert_beam') then

       write(*,*) ''
       write(*,*) '   Convert beamfiles between formats'
       write(*,*) ''
       write(*,*) '   Usage: map_editor convert_beam [fwhm] [lmax] [nmaps] [outfile] [infile]'
       write(*,*) ''
       write(*,*) '   Note that infile is optional.'
       write(*,*) ''

    else if (trim(option) == 'merge_maps') then

       write(*,*) ''
       write(*,*) '   The merge_maps operation merges two maps into one map. This operation'
       write(*,*) '   makes a new map compensating for different multipole values and beam widths.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor merge_maps [map_lowl] [beam_lowl] [map_highl] [beam_highl]'
       write(*,*) '                   [l_trans_start] [l_trans_stop] [threshold] [outmap] [prefix]'

    else if (trim(option) == 'fit_gain_offset') then

       write(*,*) ''
       write(*,*) '   The fit_gain_offset operation compares two maps for calibration. Estimates'
       write(*,*) '   for the monopole and gain differences are calculated. A residual map is'
       write(*,*) '   returned.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor fit_gain_offset [map1] [map2] [gain mask] [offset mask]'
       write(*,*) '                                     [output residual filename]'
       write(*,*) ''

    else if (trim(option) == 'fit_gain_offset_dipole') then

       write(*,*) ''
       write(*,*) '   The fit_gain_offset operation compares two maps for calibration. Estimates'
       write(*,*) '   for the monopole, gain and dipole differences are calculated. A residual'
       write(*,*) '   map is returned.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor fit_gain_offset [map1] [map2] [gain mask] [offset mask]'
       write(*,*) '                                     [output residual filename]'
       write(*,*) ''

    else if (trim(option) == 'fill_zero_mask_udgrade') then
       write(*,*) ''
       write(*,*) 'The fill_zero_mask_udgrade operation fills in new values in the input map'
       write(*,*) 'where the mask < 0.5 based on ud-grading. The mask is downgraded until there' 
       write(*,*) 'are no more "zero" pixels or nside = 1. The map is ud-graded down to the '
       write(*,*) 'given nside, excluding the pixels where mask < 0.5, and back up again to the'
       write(*,*) 'original resolution. The output map have the values where mask < 0.5 '
       write(*,*) "replaced by the ud-graded map's pixels."
       write(*,*) ''
       write(*,*) 'Usage: map_editor fill_zero_mask_udgrade [input map] [mask] [output map]'
       write(*,*) ''


    else if (trim(option) == 'print_stats' .or. trim(option) == 'stats') then
       
       write(*,*) ''
       write(*,*) '   The stats (print_stats) operation returns map statistics. A mask may be '
       write(*,*) '   applied to the input map. Statistics are printed as follows:'
       write(*,*) ''
       write(*,*) '   -------------------------'
       write(*,*) '   nside'
       write(*,*) '   nmaps'
       write(*,*) '   ordering'
       write(*,*) '   -------------------------'
       write(*,*) ''
       write(*,*) '   Polarization (1 = T, 2 = Q, 3 = U)'
       write(*,*) '   Mean'
       write(*,*) '   RMS'
       write(*,*) '   Min'
       write(*,*) '   Max'
       write(*,*) ''
       write(*,*) '   Usage: map_editor stats [input map] [mask](optional)'
       write(*,*) ''

    else if (trim(option) == 'size_info') then
       
       write(*,*) ''
       write(*,*) '   The size_info operation returns a short map size information. '
       write(*,*) '   It is similar to the stats (print_stats) operation, but without'
       write(*,*) '   statistics. The operation reads and output the following values:'
       write(*,*) ''
       write(*,*) '   nside'
       write(*,*) '   nmaps'
       write(*,*) '   ordering'
       write(*,*) ''
       write(*,*) '   Usage: map_editor size_info [input map]'
       write(*,*) ''

    else if (trim(option) == 'print_isolatitude') then
       
       write(*,*) ''
       write(*,*) '   The print_iso_latitude operation computes the mean and standard deviation'
       write(*,*) '   of the input map along iso_latitude lines of a map.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor print_isolatitude [mapfile] [mask] [outfile]'
       write(*,*) ''

    else if (trim(option) == 'print_map_to_ascii') then
       
       write(*,*) ''
       write(*,*) '   The print_map_to_ascii reads a map and prints the pixel values in ascii'
       write(*,*) '   format to the output file.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor print_map_to_ascii [mapfile in] [mask] [mapfile out]'
       write(*,*) ''

    else if (trim(option) == 'print_two_maps_to_ascii') then
       
       write(*,*) ''
       write(*,*) '   The print_two_maps_to_ascii operation reads a map and prints the '
       write(*,*) '   pixel values in ascii format to the output file.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor print_two_maps_to_ascii [mapfile1] [mapfile2]'
       write(*,*) '                                              [mask] [mapfile out]'
       write(*,*) ''
       
    else if (trim(option) == 'firas_calib') then

       write(*,*) ''
       write(*,*) '   The firas_calib operation creates a new map which is a composite of'
       write(*,*) '   absolutely calibrated FIRAS maps by weighting FIRAS maps according to the'
       write(*,*) '   input bandpass file. The output map can be used for calibration of the'
       write(*,*) '   gain and monopole of other maps with bandpasses within the 68-2911GHz range.'
       write(*,*) ''
       write(*,*) ' '
       write(*,*) '   Usage: map_editor firas_calib [bandpass file] [map weights (output)] [output map]'
       write(*,*) ' '

    else if (trim(option) == 'zero_count') then

       write(*,*) ' '
       write(*,*) '   The zero_count operation takes an input map and counts the number of '
       write(*,*) '   zero values in each map and prints it to the terminal'
       write(*,*) ' '
       write(*,*) '   Usage: map_editor zero_count [inmap]'
       write(*,*) ' '

    else if (trim(option) == 'value_count') then

       write(*,*) ' '
       write(*,*) '   The value_count operation takes an input map and counts the number of '
       write(*,*) '   pixels of the same value as the input value in each map (relative '
       write(*,*) '   difference < 1d-6). The result is printed to screen. If [value] = 0, '
       write(*,*) '   zero_count will be called instead.'
       write(*,*) ' '
       write(*,*) '   Usage: map_editor value_count [inmap] [value]'
       write(*,*) ' '

    else if (trim(option) == 'misspix_count') then

       write(*,*) ' '
       write(*,*) '   The misspix_count operation takes an input map and counts the number of '
       write(*,*) '   missing pixels in the map. The result is printed to screen.'
       write(*,*) ' '
       write(*,*) '   Usage: map_editor misspix_count [inmap]'
       write(*,*) ' '

    else if (trim(option) == 'copy_T_to_QU') then

       write(*,*) ''
       write(*,*) '   The copy_T_to_QU operation reads an input map and outputs a map of'
       write(*,*) '   nmaps = 3, where map1 is a copy of map1 (T polarization) of the input'
       write(*,*) '   map, and map 2 and 3 is copies of map1'
       write(*,*) ''
       write(*,*) '   Usage: map_editor copy_T_to_QU [input map] [output map]'
       write(*,*) ''

    else if (trim(option) == 'copy_pol_to_T') then

       write(*,*) ''
       write(*,*) '   The copy_pol_to_T operation reads an input map and outputs a map of'
       write(*,*) '   nmaps = 3, where map1 is a copy of one polarization of the input map.'
       write(*,*) ''
       write(*,*) '   Usage: map_editor copy_pol_to_T [input map] [pol] [output map]'
       write(*,*) ''

    else if (trim(option) == 'copy_pol_to_pol') then
       write(*,*) ''
       write(*,*) '   The copy_pol_to_pol operation reads two input mapa and outputs a map'
       write(*,*) '   equal to the second inputmap, where the first polarization/map-index'
       write(*,*) '   of the first inputmap is copied to the second polarization/map-index '
       write(*,*) '   of the second inputmap. The new map is written to the given output '
       write(*,*) '   mapfile.'    
       write(*,*) ''
       write(*,*) '   Usage: map_editor copy_pol_to_pol [map1] [pol1] [map2] [pol2]'
       write(*,*) '                                     [output map]'
       write(*,*) ''


    else if (trim(option) == 'median_filter_source_holes') then

       write(*,*) ''
       write(*,*) '   The median_filter_source_holes operation searches through the map for'
       write(*,*) '   pixels with significantly lower value (1e-4) than the mean within a '
       write(*,*) '   radius of [radius] arcmins of the pixel. The [treshold] value sets a '
       write(*,*) '   lower limit of the mean in order to count the "hole".'
       write(*,*) '   The operation then procedes to make a mask of 2*[radius] around the'
       write(*,*) '   pixels and creating a map where all pixels within the mask is set to'
       write(*,*) '   the median of the pixles within a radius of 2*[radius].'
       write(*,*) '   The mask is then smoothed with 30 arcmins and the output map is set to:'
       write(*,*) '      out = (1-mask)*inmap + mask*median_map'
       write(*,*) '   '

       write(*,*) ''
       write(*,*) '   Usage: map_editor median_filter_source_holes [inmap] [radius (arcmins)]'
       write(*,*) '                    [threshold] [outmap] (options)'
       write(*,*) ''
       write(*,*) '   Options:'
       write(*,*) '    -mask [mask]   (a mask to filter inside, specific filtering)'
       write(*,*) '    -output_mask   (also output the two-radii mask and the apodization mask)'
       write(*,*) ''


    else if (trim(option) == 'rescale_dust_amplitude_map') then
       write(*,*) ''
       write(*,*) 'Usage: map_editor rescale_dust_amplitude_map [ampl_inmap] [Td_map] [beta_map]'
       write(*,*) '             [old nu_ref (GHz)] [new nu_ref (GHz)] [ampl_outmap]'

       write(*,*) '   The rescale_dust_amplitude_map operation rescales the dust amplitude'
       write(*,*) '   map pixel by pixel from one reference frequency to another. '
       write(*,*) ''


    else if (trim(option) == 'median_filter_specific_values') then

       write(*,*) ''
       write(*,*) '   The median_filter_specific_values operation searches through the map'
       write(*,*) '   for pixels with a relative difference less than 1e-6 of the specified'
       write(*,*) '   value. All pixels found is then set to the median value of the pixles'
       write(*,*) '   within a radius of [radius] arcmins.'
       write(*,*) '   '
       write(*,*) 'Usage: map_editor median_filter_specific_value [inmap] [radius (arcmin)]'
       write(*,*) '                                               [value] [outmap]'
       write(*,*) ''

    else if (trim(option) == 'median_filter_misspix') then

       write(*,*) ''
       write(*,*) '   The median_filter_misspix operation searches through the map for pixels'
       write(*,*) '   with a relative difference less than 1e-5 of the missing pixel value.'
       write(*,*) '   All missing pixels found are set to the median value of the non-'
       write(*,*) '   missing pixels within a radius <= [radius] pixels.'
       write(*,*) '   It does so by checking from a radius of 1 pixel up to [radius] pixels'
       write(*,*) '   if there are more than 4 real values inside the radius. If more than 4'
       write(*,*) '   is found, the pixel value is set to be the median of these.'
       write(*,*) '   The [[dpc]] parameter is optional (logical, default .false.) used for'
       write(*,*) '   referencing to missing pixels in other maps if nmaps = 3 or 10.'
       write(*,*) '  '
       write(*,*) '  Usage: map_editor median_filter_misspix [map_in] [radius in pixels]'
       write(*,*) '                                          [map_out] [[dpc]]' 

    else if (trim(option) == 'median_filter_misspix_neighbor') then

       write(*,*) ''
       write(*,*) '   The median_filter_misspix_neighbor operation searches through the '
       write(*,*) '   map for pixels with a relative difference less than 1e-5 of the'
       write(*,*) '   missing pixel value. All missing pixels found are set to the median '
       write(*,*) '   value of the non-missing pixels of its neighbors.'
       write(*,*) '  '
       write(*,*) '  Usage: map_editor median_filter_misspix_neighbor [map_in] [map_out]'
       write(*,*) '' 

    else if (trim(option) == 'help') then

       write(*,*) ''
       write(*,*) '   Execute the program as follows to get more information about an operation:'
       write(*,*) ''
       write(*,*) '         map_editor help [operation]'
       write(*,*) ''
       write(*,*) '   For a list of all available operation, run the program without any'
       write(*,*) '   arguments.'
       write(*,*) ''


    else

       write(*,*) ''
       write(*,*) '   Execute the program as follows to get more information about an operation:'
       write(*,*) ''
       write(*,*) '         map_editor help [operation]'
       write(*,*) ''
       write(*,*) '   For a list of all available operation, run the program without any'
       write(*,*) '   arguments.'
       write(*,*) ''

    end if

  end subroutine print_help

end program map_editor
