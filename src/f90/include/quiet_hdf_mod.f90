! Due to the extreme amount of redundancy in hdf_mod, it is now generated
! by using the handy tempita template language. All the machinery for
! doing this is included in the repository, so this should just work.
module quiet_hdf_mod
  use healpix_types
  use quiet_utils
  use hdf5
  implicit none

  type hdf_file
     character(len=512) :: filename, setname
     integer(hid_t)     :: filehandle, sethandle, status
  end type hdf_file

  interface read_hdf
     module procedure read_hdf_0d_dp
     module procedure read_hdf_0d_sp
     module procedure read_hdf_0d_int
     module procedure read_hdf_1d_dp
     module procedure read_hdf_1d_sp
     module procedure read_hdf_1d_int
     module procedure read_hdf_2d_dp
     module procedure read_hdf_2d_sp
     module procedure read_hdf_2d_int
     module procedure read_hdf_3d_dp
     module procedure read_hdf_3d_sp
     module procedure read_hdf_3d_int
     module procedure read_hdf_4d_dp
     module procedure read_hdf_4d_sp
     module procedure read_hdf_4d_int
     module procedure read_hdf_5d_dp
     module procedure read_hdf_5d_sp
     module procedure read_hdf_5d_int
     module procedure read_hdf_6d_dp
     module procedure read_hdf_6d_sp
     module procedure read_hdf_6d_int
     module procedure read_hdf_7d_dp
     module procedure read_hdf_7d_sp
     module procedure read_hdf_7d_int
     module procedure read_hdf_slice_0d_dp
     module procedure read_hdf_slice_0d_sp
     module procedure read_hdf_slice_0d_int
     module procedure read_hdf_slice_1d_dp
     module procedure read_hdf_slice_1d_sp
     module procedure read_hdf_slice_1d_int
     module procedure read_hdf_slice_2d_dp
     module procedure read_hdf_slice_2d_sp
     module procedure read_hdf_slice_2d_int
     module procedure read_hdf_slice_3d_dp
     module procedure read_hdf_slice_3d_sp
     module procedure read_hdf_slice_3d_int
     module procedure read_hdf_slice_4d_dp
     module procedure read_hdf_slice_4d_sp
     module procedure read_hdf_slice_4d_int
     module procedure read_hdf_slice_5d_dp
     module procedure read_hdf_slice_5d_sp
     module procedure read_hdf_slice_5d_int
     module procedure read_hdf_slice_6d_dp
     module procedure read_hdf_slice_6d_sp
     module procedure read_hdf_slice_6d_int
     module procedure read_hdf_slice_7d_dp
     module procedure read_hdf_slice_7d_sp
     module procedure read_hdf_slice_7d_int
  end interface

  interface read_alloc_hdf
     module procedure read_alloc_hdf_1d_dp
     module procedure read_alloc_hdf_1d_sp
     module procedure read_alloc_hdf_1d_int
     module procedure read_alloc_hdf_2d_dp
     module procedure read_alloc_hdf_2d_sp
     module procedure read_alloc_hdf_2d_int
     module procedure read_alloc_hdf_3d_dp
     module procedure read_alloc_hdf_3d_sp
     module procedure read_alloc_hdf_3d_int
     module procedure read_alloc_hdf_4d_dp
     module procedure read_alloc_hdf_4d_sp
     module procedure read_alloc_hdf_4d_int
     module procedure read_alloc_hdf_5d_dp
     module procedure read_alloc_hdf_5d_sp
     module procedure read_alloc_hdf_5d_int
     module procedure read_alloc_hdf_6d_dp
     module procedure read_alloc_hdf_6d_sp
     module procedure read_alloc_hdf_6d_int
     module procedure read_alloc_hdf_7d_dp
     module procedure read_alloc_hdf_7d_sp
     module procedure read_alloc_hdf_7d_int
  end interface

  interface write_hdf
     module procedure write_hdf_0d_dp
     module procedure write_hdf_0d_sp
     module procedure write_hdf_0d_int
     module procedure write_hdf_1d_dp
     module procedure write_hdf_1d_sp
     module procedure write_hdf_1d_int
     module procedure write_hdf_2d_dp
     module procedure write_hdf_2d_sp
     module procedure write_hdf_2d_int
     module procedure write_hdf_3d_dp
     module procedure write_hdf_3d_sp
     module procedure write_hdf_3d_int
     module procedure write_hdf_4d_dp
     module procedure write_hdf_4d_sp
     module procedure write_hdf_4d_int
     module procedure write_hdf_5d_dp
     module procedure write_hdf_5d_sp
     module procedure write_hdf_5d_int
     module procedure write_hdf_6d_dp
     module procedure write_hdf_6d_sp
     module procedure write_hdf_6d_int
     module procedure write_hdf_7d_dp
     module procedure write_hdf_7d_sp
     module procedure write_hdf_7d_int
     module procedure write_hdf_slice_0d_dp
     module procedure write_hdf_slice_0d_sp
     module procedure write_hdf_slice_0d_int
     module procedure write_hdf_slice_1d_dp
     module procedure write_hdf_slice_1d_sp
     module procedure write_hdf_slice_1d_int
     module procedure write_hdf_slice_2d_dp
     module procedure write_hdf_slice_2d_sp
     module procedure write_hdf_slice_2d_int
     module procedure write_hdf_slice_3d_dp
     module procedure write_hdf_slice_3d_sp
     module procedure write_hdf_slice_3d_int
     module procedure write_hdf_slice_4d_dp
     module procedure write_hdf_slice_4d_sp
     module procedure write_hdf_slice_4d_int
     module procedure write_hdf_slice_5d_dp
     module procedure write_hdf_slice_5d_sp
     module procedure write_hdf_slice_5d_int
     module procedure write_hdf_slice_6d_dp
     module procedure write_hdf_slice_6d_sp
     module procedure write_hdf_slice_6d_int
     module procedure write_hdf_slice_7d_dp
     module procedure write_hdf_slice_7d_sp
     module procedure write_hdf_slice_7d_int
  end interface

  interface slice
     module procedure slice_0d
     module procedure slice_1d
     module procedure slice_2d
     module procedure slice_3d
     module procedure slice_4d
     module procedure slice_5d
     module procedure slice_6d
     module procedure slice_7d
  end interface

contains

  ! *****************************************************
  ! Initialization and cleanup routines
  ! *****************************************************
  subroutine initialize_quiet_hdf_mod
    implicit none
    logical(lgt), save :: initialized = .false.
    integer(i4b)       :: status
    if(initialized) return
    call h5open_f(status)
    call assert(status==0, 'quiet_hdf_mod: Could not initialize hdf module')
    initialized = .true.
  end subroutine initialize_quiet_hdf_mod

  subroutine cleanup_quiet_hdf_mod
    implicit none
    integer(i4b) :: status
    call h5close_f(status)
    call assert(status==0, 'quiet_hdf_mod: Could not close hdf module')
  end subroutine cleanup_quiet_hdf_mod

  ! *****************************************************
  ! Basic file open and close routines
  ! *****************************************************
  subroutine open_hdf_file(filename, file, mode)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=1), intent(in) :: mode
    type(hdf_file)               :: file

    ! Initialize
    call initialize_quiet_hdf_mod

    ! Open file in either read or write mode
    file%filename = filename
    if (mode == 'r') then
       call h5fopen_f(file%filename, H5F_ACC_RDONLY_F, file%filehandle, file%status)
    else if (mode == 'w') then
       call h5fcreate_f(file%filename, H5F_ACC_TRUNC_F, file%filehandle, file%status)
    else
       write(*,*) 'quiet_hdf_mod: Unknown hdf file mode =', mode
       stop
    end if

    ! Initalize sethandle to empty value
    file%setname   = ''
    file%sethandle = -1
  end subroutine open_hdf_file

  subroutine close_hdf_file(file)
    implicit none
    type(hdf_file) :: file
    call close_hdf_set(file)
    call h5fclose_f(file%filehandle, file%status)
    call assert(file%status>=0, 'quiet_hdf_mod: Could not close file')
  end subroutine close_hdf_file

  subroutine open_hdf_set(file, setname)
    implicit none
    type(hdf_file)               :: file
    character(len=*), intent(in) :: setname
    if (trim(file%setname) == trim(setname)) return
    call close_hdf_set(file)
    file%setname = setname
    call h5dopen_f(file%filehandle, file%setname, file%sethandle, file%status)
  end subroutine open_hdf_set

  subroutine close_hdf_set(file)
    implicit none
    type(hdf_file) :: file
    if (file%sethandle == -1) return
    call h5dclose_f(file%sethandle, file%status)
    call assert(file%status>=0, 'quiet_hdf_mod: Could not close set')
    file%sethandle = -1
    file%setname   = ''
  end subroutine close_hdf_set

  ! *****************************************************
  ! Query operations
  ! *****************************************************
  function get_rank_hdf(file, setname) result(rank)
    implicit none
    type(hdf_file)                :: file
    character(len=*), intent(in)  :: setname
    integer(i4b)                  :: rank
    integer(hid_t)                :: space
    call open_hdf_set(file, setname)
    call h5dget_space_f(file%sethandle, space, file%status)
    call h5sget_simple_extent_ndims_f(space, rank, file%status)
    call h5sclose_f(space, file%status)
  end function

  subroutine get_size_hdf(file, setname, ext)
    implicit none
    type(hdf_file)                  :: file
    character(len=*),   intent(in)  :: setname
    integer(i4b),       intent(out) :: ext(:)
    integer(i4b)                    :: rank
    integer(hid_t)                  :: space, n
    integer(hsize_t), allocatable, dimension(:) :: ext_hdf, mext_hdf
    call open_hdf_set(file, setname)
    call h5dget_space_f(file%sethandle, space, file%status)
    call h5sget_simple_extent_ndims_f(space, rank, file%status)
    allocate(ext_hdf(rank), mext_hdf(rank))
    call h5sget_simple_extent_dims_f(space, ext_hdf, mext_hdf, file%status)
    call h5sclose_f(space, file%status)
    n = min(size(ext),rank)
    ext(:n) = ext_hdf(:n)
    deallocate(ext_hdf, mext_hdf)
  end subroutine get_size_hdf

  ! *****************************************************
  ! Set read operations
  ! *****************************************************

  subroutine read_hdf_0d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) , intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_0d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) , intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_0d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) , intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_1d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_1d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_1d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_2d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_2d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_2d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_3d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_3d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_3d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_4d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_4d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_4d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_5d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_5d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_5d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_6d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_6d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_6d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_7d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_7d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_hdf_7d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:,:), intent(out) :: val
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_1d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:), allocatable, intent(out) :: val
    integer(i4b) :: n(1)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_1d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:), allocatable, intent(out) :: val
    integer(i4b) :: n(1)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_1d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:), allocatable, intent(out) :: val
    integer(i4b) :: n(1)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_2d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(2)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_2d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(2)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_2d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(2)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_3d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(3)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_3d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(3)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_3d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(3)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_4d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(4)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_4d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(4)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_4d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(4)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_5d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(5)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_5d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(5)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_5d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(5)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_6d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(6)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_6d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(6)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_6d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(6)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_7d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(dp) ,dimension(:,:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(7)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6),n(7)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_7d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    real(sp) ,dimension(:,:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(7)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6),n(7)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine

  subroutine read_alloc_hdf_7d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in)  :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:,:), allocatable, intent(out) :: val
    integer(i4b) :: n(7)
    if(allocated(val)) deallocate(val)
    call get_size_hdf(file, setname, n)
    allocate(val(n(1),n(2),n(3),n(4),n(5),n(6),n(7)))
    call open_hdf_set(file, setname)
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot read data from hdf set")
  end subroutine


  ! *****************************************************
  ! Set write operations
  ! *****************************************************

  subroutine write_hdf_0d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) , intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_0d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) , intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_0d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) , intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_1d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_1d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_1d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_2d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_2d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_2d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_3d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_3d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_3d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_4d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_4d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_4d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_5d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_5d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_5d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_6d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_6d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_6d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_7d_dp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(dp) ,dimension(:,:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F64LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_7d_sp(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    real(sp) ,dimension(:,:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_IEEE_F32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine

  subroutine write_hdf_7d_int(file, setname, val)
    implicit none
    type(hdf_file) :: file
    character(len=*), intent(in) :: setname
    integer(i4b) ,dimension(:,:,:,:,:,:,:), intent(in) :: val
    call create_hdf_set(file, setname, shape(val), H5T_STD_I32LE)
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, val, int(shape(val),hsize_t), file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot write data set")
  end subroutine


  ! *****************************************************
  ! Sliced set operations.
  !  These are like read/write, but the dataset is
  !  indexed with a slice. Note that the dataset must
  !  exist beforehand. Use crate_hdf_set for this.
  ! *****************************************************

  subroutine read_hdf_slice_0d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) , intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_0d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) , intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_0d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) , intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_1d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_1d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_1d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_2d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_2d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_2d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_3d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_3d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_3d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_4d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_4d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_4d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_5d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_5d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_5d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_6d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_6d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_6d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_7d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_7d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine read_hdf_slice_7d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:,:,:), intent(out) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dread_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_0d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) , intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_0d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) , intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_0d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) , intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_1d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_1d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_1d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_2d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_2d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_2d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_3d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_3d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_3d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_4d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_4d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_4d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_5d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_5d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_5d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_6d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_6d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_6d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_7d_dp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(dp) ,dimension(:,:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_DOUBLE, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_7d_sp(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    real(sp) ,dimension(:,:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_REAL, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine

  subroutine write_hdf_slice_7d_int(file, setname, slice, arr)
    implicit none
    type(hdf_file) :: file
    character(len=*),  intent(in) :: setname
    integer(i4b),      intent(in) :: slice(:,:)
    integer(i4b) ,dimension(:,:,:,:,:,:,:), intent(in) :: arr
    integer(hid_t)                :: dspace, mspace
    integer(i4b),     allocatable :: ext(:)
    integer(hsize_t)              :: hslice(3,size(slice,2))
    ! Set up data spaces for memory and disk
    call h5screate_simple_f(size(shape(arr)),int(shape(arr),hsize_t), mspace, file%status)
    call open_hdf_set(file, setname)
    allocate(ext(get_rank_hdf(file, setname)))
    call get_size_hdf(file, setname, ext)
    call h5screate_simple_f(size(ext), int(ext,hsize_t), dspace, file%status)
    ! Specify the slice
    hslice = int(parse_hdf_slice(slice, ext),hsize_t)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, hslice(1,:), hslice(2,:), &
     & file%status, stride=hslice(3,:))
    call h5dwrite_f(file%sethandle, H5T_NATIVE_INTEGER, arr, int(shape(arr),hsize_t), &
     & file%status, file_space_id=dspace, mem_space_id=mspace)
    call h5sclose_f(dspace, file%status)
    call h5sclose_f(mspace, file%status)
    deallocate(ext)
  end subroutine


  ! *****************************************************
  ! Dataset creation operation
  ! *****************************************************
  subroutine create_hdf_set(file, setname, ext, type_id)
    implicit none

    type(hdf_file)                               :: file
    character(len=*),                 intent(in) :: setname
    integer(i4b),     dimension(:),   intent(in) :: ext
    integer(i4b)                                 :: type_id
    integer(hid_t) :: space
    if (trim(file%setname) /= trim(setname)) call close_hdf_set(file)
    file%setname = setname
    call h5screate_simple_f(size(ext), int(ext,hsize_t), space, file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot create data space")
    call h5dcreate_f(file%filehandle, file%setname, type_id, space, file%sethandle, file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot create data set")
    call h5sclose_f(space, file%status)
    call assert(file%status>=0, "quiet_hdf_mod: Cannot close data space")
  end subroutine create_hdf_set

  ! Group creation. Once created, they can be used by specifying "group/dset" instaed
  ! of just "dset".
  subroutine create_hdf_group(file, group)
    implicit none
    type(hdf_file)   :: file
    character(len=*) :: group
    integer(i4b)     :: gid
    call h5gcreate_f(file%filehandle, group, gid, file%status)
    call h5gclose_f(gid, file%status)
  end subroutine

  ! **********************
  ! Helper functions
  ! **********************

  function slice_0d() result(res)
    implicit none
    integer(i4b) :: res(3,0)
    res = 0
  end function

  function slice_1d(s0) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0(:)
    integer(i4b)               :: res(3,1)
    select case(size(s0))
       case(0);  res(:,1) = [1,-1,1]
       case(1);  res(:,1) = [s0(1),s0(1),1]
       case(2);  res(:,1) = [s0(1),s0(2),1]
       case(3:); res(:,1) = s0
    end select
  end function

  function slice_2d(s0,s1) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1
    integer(i4b)               :: res(3,2)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
  end function

  function slice_3d(s0,s1,s2) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1,s2
    integer(i4b)               :: res(3,3)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
    res(:,3:3) = slice_1d(s2)
  end function

  function slice_4d(s0,s1,s2,s3) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1,s2,s3
    integer(i4b)               :: res(3,4)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
    res(:,3:3) = slice_1d(s2)
    res(:,4:4) = slice_1d(s3)
  end function

  function slice_5d(s0,s1,s2,s3,s4) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1,s2,s3,s4
    integer(i4b)               :: res(3,5)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
    res(:,3:3) = slice_1d(s2)
    res(:,4:4) = slice_1d(s3)
    res(:,5:5) = slice_1d(s4)
  end function

  function slice_6d(s0,s1,s2,s3,s4,s5) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1,s2,s3,s4,s5
    integer(i4b)               :: res(3,6)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
    res(:,3:3) = slice_1d(s2)
    res(:,4:4) = slice_1d(s3)
    res(:,5:5) = slice_1d(s4)
    res(:,6:6) = slice_1d(s5)
  end function

  function slice_7d(s0,s1,s2,s3,s4,s5,s6) result(res)
    implicit none
    integer(i4b), dimension(:) :: s0,s1,s2,s3,s4,s5,s6
    integer(i4b)               :: res(3,7)
    res(:,1:1) = slice_1d(s0)
    res(:,2:2) = slice_1d(s1)
    res(:,3:3) = slice_1d(s2)
    res(:,4:4) = slice_1d(s3)
    res(:,5:5) = slice_1d(s4)
    res(:,6:6) = slice_1d(s5)
    res(:,7:7) = slice_1d(s6)
  end function

  function parse_hdf_slice(slice, ext) result(hslice)
    implicit none
    integer(i4b), intent(in) :: slice(:,:), ext(:)
    integer(i4b)             :: hslice(3,size(slice,2)), i
    hslice = slice
    ! Negative indices count from the end, with -1 being the last valid index
    where(hslice([1,2],:) < 0) hslice([1,2],:) = hslice([1,2],:) + spread(ext,1,2) + 1
    ! We need to translate "to" into "count"
    hslice(2,:) = (hslice(2,:)-hslice(1,:)+hslice(3,:))/hslice(3,:)
    ! 0 based
    hslice(1,:) = hslice(1,:) - 1
  end function

end module quiet_hdf_mod

