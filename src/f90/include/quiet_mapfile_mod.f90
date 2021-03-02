! We were starting to get a lot of different functions for reading and
! writing maps, especially with the hdf support, so I pulled it out
! of quiet_fileutils and into this module. I also templatized it to
! make it shorter, but even so, it is pretty long. However, the interface
! is very simple. There are just two functions you need to care about here:
! read_map and write_map. Each of these have two main variants:
!
!  (read|write)_map(map,                ordering, filename)
!  (read|write)_map(map, pixels, nside, ordering, filename)
!
! Both of these handle both full and sparse maps. The difference
! is just how the information is presented to the user. If you
! want to work with full-sky maps, use the upper variant.
!
! Both FITS and HDF formats are supported. Which to use is
! determined by the file extension, or by the optional
! argument "type", which can be "hdf" or "fits".
!
! map can be 1, 2 or 3-dimensional, but the last case is
! only supported by hdf (for fits, only the first slice in
! the third direction will be written).
!
! map can be integer(i4b) or real(sp) or real(dp).
!
! Working with hdf maps is much faster than fits maps, by
! more than a factor of 10 in my tests. It is possible
! to make the fits part faster (but still slower than hdf)
! by calling read_bintab, but I did not think this was
! worth it, as it would make things less elegant, and
! we will want to use hdf in any case.
module quiet_mapfile_mod
  use quiet_hdf_mod
  use fitstools
  use head_fits
  use quiet_utils
  implicit none

  interface write_map
     module procedure write_sparse_map_1d_dp
     module procedure write_sparse_map_2d_dp
     module procedure write_sparse_map_3d_dp
     module procedure write_sparse_map_1d_sp
     module procedure write_sparse_map_2d_sp
     module procedure write_sparse_map_3d_sp
     module procedure write_sparse_map_1d_int
     module procedure write_sparse_map_2d_int
     module procedure write_sparse_map_3d_int
     module procedure write_full_map_1d_dp
     module procedure write_full_map_2d_dp
     module procedure write_full_map_3d_dp
     module procedure write_full_map_1d_sp
     module procedure write_full_map_2d_sp
     module procedure write_full_map_3d_sp
     module procedure write_full_map_1d_int
     module procedure write_full_map_2d_int
     module procedure write_full_map_3d_int
  end interface

  interface read_map
     module procedure read_sparse_map_1d_dp
     module procedure read_sparse_map_2d_dp
     module procedure read_sparse_map_3d_dp
     module procedure read_sparse_map_1d_sp
     module procedure read_sparse_map_2d_sp
     module procedure read_sparse_map_3d_sp
     module procedure read_sparse_map_1d_int
     module procedure read_sparse_map_2d_int
     module procedure read_sparse_map_3d_int
     module procedure read_full_map_1d_dp
     module procedure read_full_map_2d_dp
     module procedure read_full_map_3d_dp
     module procedure read_full_map_1d_sp
     module procedure read_full_map_2d_sp
     module procedure read_full_map_3d_sp
     module procedure read_full_map_1d_int
     module procedure read_full_map_2d_int
     module procedure read_full_map_3d_int
  end interface


contains

  subroutine write_sparse_map_1d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_1d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_1d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_2d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_2d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_2d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_3d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_3d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_3d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_1d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_1d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_1d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_2d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_2d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_2d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_3d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_3d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_3d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_1d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_1d_int(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_1d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_2d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_2d_int(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_2d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_sparse_map_3d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_sparse_hdf_map_3d_int(map, pixels, nside, ordering, filename)
       case("fits"); call write_sparse_fits_map_3d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_1d_dp(map, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_1d_dp(map, ordering, filename)
       case("fits"); call write_full_fits_map_1d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_2d_dp(map, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_2d_dp(map, ordering, filename)
       case("fits"); call write_full_fits_map_2d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_3d_dp(map, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_3d_dp(map, ordering, filename)
       case("fits"); call write_full_fits_map_3d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_1d_sp(map, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_1d_sp(map, ordering, filename)
       case("fits"); call write_full_fits_map_1d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_2d_sp(map, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_2d_sp(map, ordering, filename)
       case("fits"); call write_full_fits_map_2d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_3d_sp(map, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_3d_sp(map, ordering, filename)
       case("fits"); call write_full_fits_map_3d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_1d_int(map, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_1d_int(map, ordering, filename)
       case("fits"); call write_full_fits_map_1d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_2d_int(map, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_2d_int(map, ordering, filename)
       case("fits"); call write_full_fits_map_2d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine write_full_map_3d_int(map, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call write_full_hdf_map_3d_int(map, ordering, filename)
       case("fits"); call write_full_fits_map_3d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_1d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_1d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_1d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_2d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_2d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_2d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_3d_dp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_3d_dp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_3d_dp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_1d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_1d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_1d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_2d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_2d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_2d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_3d_sp(map, pixels, nside, ordering, filename, type)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_3d_sp(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_3d_sp(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_1d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_1d_int(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_1d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_2d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_2d_int(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_2d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_sparse_map_3d_int(map, pixels, nside, ordering, filename, type)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable :: pixels(:)
    integer(i4b)      :: nside, ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    select case(get_map_type(filename, type))
       case("hdf");  call read_sparse_hdf_map_3d_int(map, pixels, nside, ordering, filename)
       case("fits"); call read_sparse_fits_map_3d_int(map, pixels, nside, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
  end subroutine
  subroutine read_full_map_1d_dp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(dp) ,dimension(:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_1d_dp(map, ordering, filename)
       case("fits"); call read_full_fits_map_1d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = 1
  end subroutine
  subroutine read_full_map_2d_dp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_2d_dp(map, ordering, filename)
       case("fits"); call read_full_fits_map_2d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine
  subroutine read_full_map_3d_dp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_3d_dp(map, ordering, filename)
       case("fits"); call read_full_fits_map_3d_dp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine
  subroutine read_full_map_1d_sp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(sp) ,dimension(:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_1d_sp(map, ordering, filename)
       case("fits"); call read_full_fits_map_1d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = 1
  end subroutine
  subroutine read_full_map_2d_sp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_2d_sp(map, ordering, filename)
       case("fits"); call read_full_fits_map_2d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine
  subroutine read_full_map_3d_sp(map, ordering, filename, type, nside, nmap)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_3d_sp(map, ordering, filename)
       case("fits"); call read_full_fits_map_3d_sp(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine
  subroutine read_full_map_1d_int(map, ordering, filename, type, nside, nmap)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_1d_int(map, ordering, filename)
       case("fits"); call read_full_fits_map_1d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = 1
  end subroutine
  subroutine read_full_map_2d_int(map, ordering, filename, type, nside, nmap)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_2d_int(map, ordering, filename)
       case("fits"); call read_full_fits_map_2d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine
  subroutine read_full_map_3d_int(map, ordering, filename, type, nside, nmap)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map
    integer(i4b)      :: ordering
    character(len=*)  :: filename
    character(len=*), optional :: type
    integer(i4b),     optional :: nside, nmap
    select case(get_map_type(filename, type))
       case("hdf");  call read_full_hdf_map_3d_int(map, ordering, filename)
       case("fits"); call read_full_fits_map_3d_int(map, ordering, filename)
       case default; call assert(.false., "Unknown type " // get_map_type(filename,type))
    end select
    if(present(nside)) nside = npix2nside(size(map,1))
    if(present(nmap)) nmap = size(map,2)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !             Support routines                 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_sparse_fits_map_1d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = 1
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcld(unit, i+1, 1, 1, size(map,1), map(:), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_2d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcld(unit, i+1, 1, 1, size(map,1), map(:,i), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_3d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcld(unit, i+1, 1, 1, size(map,1), map(:,i,1), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_1d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = 1
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcle(unit, i+1, 1, 1, size(map,1), map(:), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_2d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcle(unit, i+1, 1, 1, size(map,1), map(:,i), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_3d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpcle(unit, i+1, 1, 1, size(map,1), map(:,i,1), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_1d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = 1
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpclj(unit, i+1, 1, 1, size(map,1), map(:), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_2d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpclj(unit, i+1, 1, 1, size(map,1), map(:,i), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_sparse_fits_map_3d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, status, unit, ordering, i, nmap
    character(len=*)  :: filename
    nmap = size(map,2)
    call write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    status = 0
    do i = 1, nmap
       call ftpclj(unit, i+1, 1, 1, size(map,1), map(:,i,1), status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine write_full_fits_map_1d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = 1
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:), [size(map,1),1])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_2d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:,:), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_3d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:,:,1), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_1d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = 1
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:), [size(map,1),1])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_2d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:,:), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_3d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab((reshape(map(:,:,1), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_1d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = 1
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab(dble(reshape(map(:), [size(map,1),1])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_2d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab(dble(reshape(map(:,:), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_full_fits_map_3d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:),     intent(in) :: map
    integer(i4b),                       intent(in) :: ordering
    character(len=*),                   intent(in) :: filename

    integer(i4b)     :: npix, nmaps, i, nside
    logical(lgt)     :: exist, polarization
    character(len=6) :: order
    character(len=80), dimension(1:120)    :: header

    npix         = size(map,1)
    nside        = npix2nside(npix)
    nmaps = size(map,2)
    polarization = (nmaps == 3)
    order = 'RING'; if(ordering == 2) order = 'NESTED'
    header = ""
    call add_card(header,"COMMENT","*************************************")
    ! start putting information relative to this code and run
    call write_minimal_header(header, 'map', append=.true., &
         nside = nside, ordering = order, &
         polar = polarization )
    call add_card(header,"COMMENT","*************************************")
    call write_bintab(dble(reshape(map(:,:,1), [size(map,1),size(map,2)])), &
         & npix, nmaps, header, size(header), "!" // trim(filename))
  end subroutine
  subroutine write_sparse_hdf_map_1d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_2d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_3d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_1d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_2d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_3d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_1d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_2d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_sparse_hdf_map_3d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:) :: map
    integer(i4b)      :: pixels(:), nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 0)
    call write_hdf(file, "pixels", pixels)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_1d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_2d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_3d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_1d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_2d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_3d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_1d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:), [size(map,1),1,1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_2d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:), [size(map,1),size(map,2),1]))
    call close_hdf_file(file)
  end subroutine
  subroutine write_full_hdf_map_3d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:) :: map
    integer(i4b)      :: nside, ordering, i, nmap
    character(len=*)  :: filename
    type(hdf_file)    :: file
    nside = npix2nside(size(map,1))
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "nside", nside)
    call write_hdf(file, "ordering", ordering)
    call write_hdf(file, "fullsky", 1)
    call write_hdf(file, "maps", reshape(map(:,:,:), [size(map,1),size(map,2),size(map,3)]))
    call close_hdf_file(file)
  end subroutine
  subroutine read_sparse_fits_map_1d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n))
    nmap = 1
    do i = 1, nmap
       call ftgcvd(unit, i+off, 1, 1, n, hpx_dbadval, map(:), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_2d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap))
    
    do i = 1, nmap
       call ftgcvd(unit, i+off, 1, 1, n, hpx_dbadval, map(:,i), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_3d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap,1))
    
    do i = 1, nmap
       call ftgcvd(unit, i+off, 1, 1, n, hpx_dbadval, map(:,i,1), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_1d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n))
    nmap = 1
    do i = 1, nmap
       call ftgcve(unit, i+off, 1, 1, n, hpx_sbadval, map(:), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_2d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap))
    
    do i = 1, nmap
       call ftgcve(unit, i+off, 1, 1, n, hpx_sbadval, map(:,i), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_3d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap,1))
    
    do i = 1, nmap
       call ftgcve(unit, i+off, 1, 1, n, hpx_sbadval, map(:,i,1), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_1d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n))
    nmap = 1
    do i = 1, nmap
       call ftgcvj(unit, i+off, 1, 1, n, -1, map(:), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_2d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap))
    
    do i = 1, nmap
       call ftgcvj(unit, i+off, 1, 1, n, -1, map(:,i), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_sparse_fits_map_3d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, unit, nmap, n, ordering, status
    integer(i4b) :: i, off
    logical(lgt) :: anynull
    character(len=*)  :: filename
    call read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    n = size(pixels)
    status = 0
    allocate(map(n,nmap,1))
    
    do i = 1, nmap
       call ftgcvj(unit, i+off, 1, 1, n, -1, map(:,i,1), anynull, status)
       if (status > 0) call printerror(status)
    end do
    call ftclos(unit, status)
    if (status > 0) call printerror(status)
  end subroutine
  subroutine read_full_fits_map_1d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_1d_dp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1))
    map = hpx_dbadval
    map(pixels) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_2d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_2d_dp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2)))
    map = hpx_dbadval
    map(pixels,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_3d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_3d_dp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2),size(tmap,3)))
    map = hpx_dbadval
    map(pixels,:,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_1d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_1d_sp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1))
    map = hpx_sbadval
    map(pixels) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_2d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_2d_sp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2)))
    map = hpx_sbadval
    map(pixels,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_3d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_3d_sp(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2),size(tmap,3)))
    map = hpx_sbadval
    map(pixels,:,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_1d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_1d_int(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1))
    map = -1
    map(pixels) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_2d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_2d_int(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2)))
    map = -1
    map(pixels,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_full_fits_map_3d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map, tmap
    integer(i4b), allocatable, dimension(:) :: pixels
    integer(i4b) :: nside, unit, nmap, ordering, npix, i
    character(len=*) :: filename
    call read_sparse_fits_map_3d_int(tmap, pixels, nside, ordering, filename)
    npix = 12*nside**2
    allocate(map(0:npix-1,size(tmap,2),size(tmap,3)))
    map = -1
    map(pixels,:,:) = tmap
    deallocate(tmap,pixels)
  end subroutine
  subroutine read_sparse_hdf_map_1d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_2d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_3d_dp(map, pixels, nside, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2),ext(3)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_1d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_2d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_3d_sp(map, pixels, nside, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2),ext(3)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_1d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_2d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_sparse_hdf_map_3d_int(map, pixels, nside, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    if(fullsky == 1) then
       npix = 12*nside**2
       allocate(pixels(npix))
       pixels = irange(npix)-1
    else
       call get_size_hdf(file, "pixels", ext); npix = ext(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
    end if
    call get_size_hdf(file, "maps", ext)
    allocate(map(ext(1),ext(2),ext(3)))
    call read_hdf(file, "maps", map)
  end subroutine
  subroutine read_full_hdf_map_1d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:), allocatable :: map
    real(dp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_dbadval
       call read_hdf(file, "maps", tmap)
       map(pixels) = reshape(tmap(:,1,1), [size(tmap,1)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_2d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:), allocatable :: map
    real(dp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_dbadval
       call read_hdf(file, "maps", tmap)
       map(pixels,:) = reshape(tmap(:,:,1), [size(tmap,1),size(tmap,2)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_3d_dp(map, ordering, filename)
    implicit none
    real(dp) ,dimension(:,:,:), allocatable :: map
    real(dp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2),ext(3)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2),ext(3)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_dbadval
       call read_hdf(file, "maps", tmap)
       map(pixels,:,:) = reshape(tmap(:,:,:), [size(tmap,1),size(tmap,2),size(tmap,3)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_1d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:), allocatable :: map
    real(sp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_sbadval
       call read_hdf(file, "maps", tmap)
       map(pixels) = reshape(tmap(:,1,1), [size(tmap,1)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_2d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:), allocatable :: map
    real(sp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_sbadval
       call read_hdf(file, "maps", tmap)
       map(pixels,:) = reshape(tmap(:,:,1), [size(tmap,1),size(tmap,2)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_3d_sp(map, ordering, filename)
    implicit none
    real(sp) ,dimension(:,:,:), allocatable :: map
    real(sp), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2),ext(3)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2),ext(3)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = hpx_sbadval
       call read_hdf(file, "maps", tmap)
       map(pixels,:,:) = reshape(tmap(:,:,:), [size(tmap,1),size(tmap,2),size(tmap,3)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_1d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:), allocatable :: map
    integer(i4b), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = -1
       call read_hdf(file, "maps", tmap)
       map(pixels) = reshape(tmap(:,1,1), [size(tmap,1)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_2d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:), allocatable :: map
    integer(i4b), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = -1
       call read_hdf(file, "maps", tmap)
       map(pixels,:) = reshape(tmap(:,:,1), [size(tmap,1),size(tmap,2)])
       deallocate(pixels, tmap)
    end if
  end subroutine
  subroutine read_full_hdf_map_3d_int(map, ordering, filename)
    implicit none
    integer(i4b) ,dimension(:,:,:), allocatable :: map
    integer(i4b), dimension(:,:,:), allocatable :: tmap
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b)      :: nside, nmap, n, ordering, fullsky, npix, ext(3), ext2(1)
    character(len=*)  :: filename
    type(hdf_file)    :: file
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside", nside)
    call read_hdf(file, "ordering", ordering)
    call read_hdf(file, "fullsky", fullsky)
    call get_size_hdf(file, "maps", ext)
    if(fullsky == 1) then
       allocate(map(0:ext(1)-1,ext(2),ext(3)))
       call read_hdf(file, "maps", map)
    else
       call get_size_hdf(file, "pixels", ext2); npix = ext2(1)
       allocate(pixels(npix))
       call read_hdf(file, "pixels", pixels)
       allocate(map(0:12*nside**2-1,ext(2),ext(3)))
       allocate(tmap(ext(1),ext(2),ext(3)))
       map = -1
       call read_hdf(file, "maps", tmap)
       map(pixels,:,:) = reshape(tmap(:,:,:), [size(tmap,1),size(tmap,2),size(tmap,3)])
       deallocate(pixels, tmap)
    end if
  end subroutine


  !!!!!!!!!!!!!!!!!!!
  ! Helper routines !
  !!!!!!!!!!!!!!!!!!!

  subroutine write_sparse_helper(nmap, pixels, nside, ordering, filename, unit)
    implicit none
    integer(i4b) :: pixels(:), nside, status, unit, nmap, n, tfields, ordering
    integer(i4b) :: i, itn
    character(len=20), allocatable, dimension(:) :: ttype, tform, tunit
    character(len=80) :: header(120), comment, card
    character(len=6)  :: order
    character(len=2)  :: stn
    character(Len=*)  :: filename

    n    = size(pixels)
    allocate(ttype(nmap+1), tform(nmap+1), tunit(nmap+1))

    if(ordering == 1) then; order = "RING"; else; order = "NESTED"; end if
    tfields  = nmap+1 ! An extra field for the indices
    ttype  = 'unknown'; tunit  = ''; tform(1) = '1J'; tform(2:) = '1D'
    header = ''
    status = 0
    unit = getlun()

    call ftinit(unit,"!"//trim(filename),1,status)
    call ftphps(unit,32,0,0,status)
    call ftpdat(unit,status) ! format (yyyy-mm-dd)
    call ftcrhd(unit,status)
    if (status > 0) call printerror(status)
    call ftphbn(unit, n, nmap+1, ttype, tform, tunit, '', 0, status)
    if (status > 0) call printerror(status)

    call add_card(header,"COMMENT","*************************************")
    call write_minimal_header(header, 'cutmap', append=.true., &
         nside = nside, ordering = order, &
         polar = nmap==3)
    call add_card(header,"COMMENT","*************************************")

    comment = ''
    ! Hack: merge header with the information we already put there
    do i=1,size(header)
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPE1 is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          if (itn > tfields) cycle
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          call putrec(unit,header(i), status)    ! write new TTYPE1
          call ftpkys(unit,'TFORM'//stn,tform(1),comment,status) ! and write new TFORM1 right after
         if (status > 0) call printerror(status)
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
          if (status > 0) call printerror(status)
       endif
    enddo

    call ftpclj(unit, 1, 1, 1, n, pixels, status)
    deallocate(ttype, tform, tunit)
  end subroutine

  subroutine read_sparse_helper(nmap, pixels, nside, ordering, filename, unit, off)
    implicit none
    integer(i4b), allocatable, dimension(:)   :: pixels
    integer(i4b) :: nside, status, unit, nmap, n, tfields, ordering
    integer(i4b) :: i, off, nfield, repeat, width, datacode
    logical(lgt) :: anynull, sparse
    character(len=80) :: comment
    character(len=20) :: indexing, order, tform
    character(Len=*)  :: filename

    unit = getlun()
    status = 0
    call ftnopn(unit,trim(filename) // "[1]",0,status)
    if (status > 0) call printerror(status)

    call ftgkyj(unit,"NSIDE",   nside,   comment,status)
    if (status > 0) call printerror(status)
    call ftgkyj(unit,"NAXIS2",  n,       comment,status)
    if (status > 0) call printerror(status)
    call ftgkyj(unit,"TFIELDS", nfield,  comment,status)
    if (status > 0) call printerror(status)
    call ftgkys(unit,"INDXSCHM",indexing,comment,status)
	if (status > 0) then; indexing = 'IMPLICIT'; status = 0; end if
    call ftgkys(unit,"ORDERING",order,   comment,status)
    if (status > 0) call printerror(status)

    if(indexing == 'IMPLICIT') then
       sparse = .false.
       nmap   = nfield
       off    = 0
    elseif(indexing == 'EXPLICIT') then
       sparse = .true.
       nmap   = nfield-1
       off    = 1
    else
       write(*,*) "read_sparse_map got unrecognized indexing sheme: " // trim(indexing) // "!"
       stop
    end if
    if(order == "RING") then; ordering = 1; else; ordering = 2; end if

    ! The full format uses stupid repeats, so we need to check gat
    call ftgkys(unit,"TFORM" // trim(itoa(off+1)), tform, comment,status)
    if (status > 0) call printerror(status)
    call ftbnfm(tform, datacode, repeat, width, status)
    if (status > 0) call printerror(status)
    n = n*repeat

    ! Ok, time to read the actual data.
    allocate(pixels(n))

    if(sparse) then
       call ftgcvj(unit, 1, 1, 1, n, 0, pixels, anynull, status)
       if (status > 0) call printerror(status)
    else
       do i = 1, n; pixels(i) = i-1; end do
    end if
  end subroutine

  function get_map_type(filename, itype) result(type)
    implicit none
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: itype
    character(len=512)                     :: type
    integer(i4b)                           :: n, i
    if(present(itype)) then
       type = itype
    else
       n = len_trim(filename)
       do i = n,1,-1
          if(filename(i:i) == ".") exit
       end do
       type = filename(i+1:n)
    end if
  end function

end module

