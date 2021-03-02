cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(Cosmotools
	LANGUAGES C CXX Fortran
	)

# Looking for MPI
find_package(MPI COMPONENTS Fortran)
# Looking for OpenMP
find_package(OpenMP)

# TODO: Fix FindHEALPix.cmake file to correctly link HEALPix library and modules
set(HEALPIX_LIBRARIES
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/lib/libhealpix.a"	
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/lib/libsharp.a"
	)
set(HEALPIX_INCLUDES
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/include"
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/include/libsharp"
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/build/mod"
	)
# TODO: Include Find CFitsIO module
set(CFITSIO_LIBRARIES
	"/mn/stornext/u3/maksymb/commander/quiet/build/install/lib/libcfitsio.a"	
	)

# setting the directory where to output all .mod and .o files
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/mod)
#================================================================================
# Setting the Fortran source directory location
set(COSMOTOOLS_Fortran_SRC						"${CMAKE_SOURCE_DIR}/src/f90")
set(COSMOTOOLS_Fortran_SRC_INCLUDE		"${COSMOTOOLS_Fortran_SRC}/include")
set(COSMOTOOLS_Fortran_SRC_MAP_EDITOR "${COSMOTOOLS_Fortran_SRC}/map_editor")
set(COSMOTOOLS_Fortran_SRC_MAPTOOL		"${COSMOTOOLS_Fortran_SRC}/maptool")
set(COSMOTOOLS_Fortran_SRC_POSTMAP		"${COSMOTOOLS_Fortran_SRC}/postmap")
set(COSMOTOOLS_Fortran_SRC_UTILS			"${COSMOTOOLS_Fortran_SRC}/utils")
# Setting the C++ source directory location
set(COSMOTOOLS_CXX_SRC			"${CMAKE_SOURCE_DIR}/src/cpp")
set(COSMOTOOLS_Python_SRC		"${CMAKE_SOURCE_DIR}/src/python")
# Tempita
set(TEMPITA_DIR "${COSMOTOOLS_Python_SRC}")
#================================================================================
# TEMPITA
#================================================================================
# Processing *.in files with Tempita language to obtain *.f90 files
#add_custom_target(tempita ALL "")
list(APPEND tempita_in_files
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_covfile_mod.f90.in"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_hdf_mod.f90.in"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mapfile_mod.f90.in"
	)
list(APPEND tempita_out_files
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_covfile_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_hdf_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mapfile_mod.f90"
	)
foreach(in_file out_file IN ZIP_LISTS tempita_in_files tempita_out_files)
	execute_process(
		COMMAND ${TEMPITA_DIR}/tempita_proc.py < $< > $@ 
		INPUT_FILE "${in_file}"
		OUTPUT_FILE "${out_file}"
		)
	#message("${in_file}")
	#message("${out_file}")
endforeach()
#================================================================================
# LIBQUIET -- compile everything from include into QUIET library 
#================================================================================
# Compile math_tools first -- serves as one of the main inputs to other codes
add_library(quiet SHARED
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/math_tools.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/d1mach.f"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/drc3jj.f"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/fake_mpi_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/l1_read_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/locate_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/minmax_tools.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/powell_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quasi_newton_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_acceptlist_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_angcorr_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_apex_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_assembly_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_ces_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_constrained_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_defs.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_detector_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_ephem_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_fft_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_fileutils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_filter_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_gain_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_glitchy_noise_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_healpix_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_i2qu_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_l1_defs.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_Lx_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mask_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_module_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mpi_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mpi_utils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_noise_estimation_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_nr_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_patch_detect_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_patch_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_pixspace_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_pmac_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_pointing_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_pointing_mod2.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_postutils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_shared_output_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_sidelobe_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_stat_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_status_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_system_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_target_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_task_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_todsim_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_typeb_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_utils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/scalautils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/scalawrap.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/sdeOOF.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/sort_utils.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/spline_1D_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/spline_2D_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/typeb_correction_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/wmap_beam_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/ziggurat.f90"
	#
	#"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_ephem_backend.cpp"
	#"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_system_backend.cpp"
	# Tempita processed files
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_covfile_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_hdf_mod.f90"
	"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mapfile_mod.f90"
	#"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_covfile_mod.f90.in"
	#"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_hdf_mod.f90.in"
	#"${COSMOTOOLS_Fortran_SRC_INCLUDE}/quiet_mapfile_mod.f90.in"
	)
# Including HEALPix Fortran *.mod and *.o files
target_include_directories(quiet
	PUBLIC
	"${HEALPIX_INCLUDES}"
	)

target_link_libraries(quiet
	# Including MPI
	MPI::MPI_Fortran
	# Including OpenMP
	OpenMP::OpenMP_Fortran
	# Including HEALPix + Sharp2
	"${HEALPIX_LIBRARIES}"
	# Including CFitsIO (as HEALPix needs it)
	#"${CFITSIO_LIBRARIES}"	
	)
install(TARGETS quiet ARCHIVE DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#================================================================================

# Adding map_editor targets
#add_executable(map_editor
#	"${COSMOTOOLS_Fortran_SRC_MAP_EDITOR}/map_editor.f90"
#	"${COSMOTOOLS_Fortran_SRC_MAP_EDITOR}/map_editor_complex_ops_mod.f90"
#	"${COSMOTOOLS_Fortran_SRC_MAP_EDITOR}/map_editor_simple_ops_mod.f90"
#	"${COSMOTOOLS_Fortran_SRC_MAP_EDITOR}/map_editor_utils.f90"
#	)
#
## Including HEALPix Fortran *.mod and *.o files
#target_include_directories(map_editor
#	PUBLIC
#	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/include"
#	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/include/libsharp"
#	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/build/mod"
#	# Including math_tools.f90
#	"/mn/stornext/u3/maksymb/cmake_projects/Cosmotools/src/f90/include"
#	)
#
#target_link_libraries(map_editor
#	# Including MPI
#	MPI::MPI_Fortran
#	# Including OpenMP
#	OpenMP::OpenMP_Fortran
#	# Including HEALPix + Sharp2
#	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/lib/libhealpix.a"	
#	"/mn/stornext/u3/maksymb/commander/quiet/build/install/healpix/lib/libsharp.a"
#	# Including CFitsIO (as HEALPix needs it)
#	#"/mn/stornext/u3/maksymb/commander/quiet/build/install/lib/libcfitsio.a"	
#	#
#	math_tools
#	)
