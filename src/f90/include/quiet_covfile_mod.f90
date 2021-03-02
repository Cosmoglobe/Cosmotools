module quiet_covfile_mod
  use healpix_types
  use quiet_hdf_mod
  implicit none

  interface read_covmat
     module procedure read_covmat_2
     module procedure read_covmat_4
  end interface

  interface write_covmat
     module procedure write_covmat_2
     module procedure write_covmat_4
  end interface


contains

  ! I want to be able to read/write covmatrices/equation sets
  ! with a consistent interface no matter if they are unf or hdf.
  ! So I will implement RW for hdf, and then wrapper functions which
  ! choose the appropriate one.

  subroutine read_covmat_2(covmat, pixels, nside, order, ncomp, fname, inv, type, verbose)
    implicit none
    real(dp), allocatable         :: covmat(:,:)
    integer(i4b), allocatable     :: pixels(:)
    integer(i4b)                  :: ncomp, nside, order
    logical(lgt),     optional    :: inv, verbose
    character(len=*)              :: fname
    character(len=*), optional    :: type
    character(len=32)             :: typ_
    typ_ = get_extension(fname, type)
    select case(typ_)
       case("hdf"); call read_covmat_2_hdf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
       case("unf"); call read_covmat_2_unf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
       case default; stop "Unknown covmat filetype."
    end select
  end subroutine

  subroutine read_covmat_4(covmat, pixels, nside, order,  fname, inv, type, verbose)
    implicit none
    real(dp), allocatable         :: covmat(:,:,:,:)
    integer(i4b), allocatable     :: pixels(:)
    integer(i4b)                  ::  nside, order
    logical(lgt),     optional    :: inv, verbose
    character(len=*)              :: fname
    character(len=*), optional    :: type
    character(len=32)             :: typ_
    typ_ = get_extension(fname, type)
    select case(typ_)
       case("hdf"); call read_covmat_4_hdf(covmat, pixels, nside, order,  fname, inv, verbose)
       case("unf"); call read_covmat_4_unf(covmat, pixels, nside, order,  fname, inv, verbose)
       case default; stop "Unknown covmat filetype."
    end select
  end subroutine

  subroutine write_covmat_2(covmat, pixels, nside, order, ncomp, fname, inv, type, verbose)
    implicit none
    real(dp)         :: covmat(:,:)
    integer(i4b)     :: pixels(:)
    integer(i4b)                  :: ncomp, nside, order
    logical(lgt),     optional    :: inv, verbose
    character(len=*)              :: fname
    character(len=*), optional    :: type
    character(len=32)             :: typ_
    typ_ = get_extension(fname, type)
    select case(typ_)
       case("hdf"); call write_covmat_2_hdf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
       case("unf"); call write_covmat_2_unf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
       case default; stop "Unknown covmat filetype."
    end select
  end subroutine

  subroutine write_covmat_4(covmat, pixels, nside, order,  fname, inv, type, verbose)
    implicit none
    real(dp)         :: covmat(:,:,:,:)
    integer(i4b)     :: pixels(:)
    integer(i4b)                  ::  nside, order
    logical(lgt),     optional    :: inv, verbose
    character(len=*)              :: fname
    character(len=*), optional    :: type
    character(len=32)             :: typ_
    typ_ = get_extension(fname, type)
    select case(typ_)
       case("hdf"); call write_covmat_4_hdf(covmat, pixels, nside, order,  fname, inv, verbose)
       case("unf"); call write_covmat_4_unf(covmat, pixels, nside, order,  fname, inv, verbose)
       case default; stop "Unknown covmat filetype."
    end select
  end subroutine

  subroutine read_covmat_2_hdf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
    implicit none
    real(dp), allocatable         :: covmat(:,:)
    integer(i4b), allocatable     :: pixels(:)
    integer(i4b)                  :: ncomp, nside, order, ext(4), inv_
    character(len=*)              :: fname
    logical(lgt),     optional    :: inv, verbose
	logical(lgt)                  :: verb
    type(hdf_file)                :: hfile
    verb = .false.; if(present(verbose)) verb = verbose
    if(allocated(covmat)) deallocate(covmat)
    if(allocated(pixels)) deallocate(pixels)
    call open_hdf_file(fname, hfile, "r")
    call get_size_hdf(hfile, "cov", ext)
    call open_hdf_file(fname, hfile, "r")
    call read_hdf(hfile, "nside",    nside)
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    call read_hdf(hfile, "ordering", order)
    if(verb) write(*,'(a10,i5)') "order: ", order
    call read_hdf(hfile, "inverse", inv_)
    if(verb) write(*,'(a10,l)') "inverse: ", inv_
    allocate(covmat(ext(1)*ext(2),ext(3)*ext(4)), pixels(ext(1)))
    if(verb) write(*,'(a10,4i5)') "shape: ", ext
    ncomp = size(covmat,1)/size(pixels)
    call open_hdf_set(hfile, "cov")
    call h5dread_f(hfile%sethandle, H5T_NATIVE_DOUBLE, covmat, int(ext,hsize_t), hfile%status)
    call read_hdf(hfile, "pixels",   pixels)
    call close_hdf_file(hfile)
    if(present(inv)) inv = inv_ /= 0
  end subroutine

  subroutine read_covmat_4_hdf(covmat, pixels, nside, order,  fname, inv, verbose)
    implicit none
    real(dp), allocatable         :: covmat(:,:,:,:)
    integer(i4b), allocatable     :: pixels(:)
    integer(i4b)                  ::  nside, order, ext(4), inv_
    character(len=*)              :: fname
    logical(lgt),     optional    :: inv, verbose
	logical(lgt)                  :: verb
    type(hdf_file)                :: hfile
    verb = .false.; if(present(verbose)) verb = verbose
    if(allocated(covmat)) deallocate(covmat)
    if(allocated(pixels)) deallocate(pixels)
    call open_hdf_file(fname, hfile, "r")
    call get_size_hdf(hfile, "cov", ext)
    call open_hdf_file(fname, hfile, "r")
    call read_hdf(hfile, "nside",    nside)
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    call read_hdf(hfile, "ordering", order)
    if(verb) write(*,'(a10,i5)') "order: ", order
    call read_hdf(hfile, "inverse", inv_)
    if(verb) write(*,'(a10,l)') "inverse: ", inv_
    allocate(covmat(ext(1),ext(2),ext(3),ext(4)), pixels(ext(1)))
    if(verb) write(*,'(a10,4i5)') "shape: ", ext
    call read_hdf(hfile, "cov",      covmat)
    call read_hdf(hfile, "pixels",   pixels)
    call close_hdf_file(hfile)
    if(present(inv)) inv = inv_ /= 0
  end subroutine

  subroutine write_covmat_2_hdf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
    implicit none
    real(dp)         :: covmat(:,:)
    integer(i4b)     :: pixels(:)
    integer(i4b)                  :: ncomp, nside, order, ext(4), inv_
    character(len=*)              :: fname
    logical(lgt),     optional    :: inv, verbose
	logical(lgt)                  :: verb
    type(hdf_file)                :: hfile
    verb = .false.; if(present(verbose)) verb = verbose
    ext = [size(covmat,1)/ncomp,ncomp,size(covmat,2)/ncomp,ncomp]
    call open_hdf_file(fname, hfile, "w")
    call write_hdf(hfile, "nside",    nside)
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    call write_hdf(hfile, "ordering", order)
    if(verb) write(*,'(a10,i5)') "order: ", order
    call write_hdf(hfile, "inverse", inv_)
    if(verb) write(*,'(a10,l)') "inverse: ", inv_
    if(verb) write(*,'(a10,4i5)') "shape: ", ext
    call create_hdf_set(hfile, "cov", ext, H5T_IEEE_F64LE)
    call open_hdf_set(hfile, "cov")
    call h5dwrite_f(hfile%sethandle, H5T_NATIVE_DOUBLE, covmat, int(ext,hsize_t), hfile%status)
    call write_hdf(hfile, "pixels",   pixels)
    call close_hdf_file(hfile)
    if(present(inv)) inv = inv_ /= 0
  end subroutine

  subroutine write_covmat_4_hdf(covmat, pixels, nside, order,  fname, inv, verbose)
    implicit none
    real(dp)         :: covmat(:,:,:,:)
    integer(i4b)     :: pixels(:)
    integer(i4b)                  ::  nside, order, ext(4), inv_
    character(len=*)              :: fname
    logical(lgt),     optional    :: inv, verbose
	logical(lgt)                  :: verb
    type(hdf_file)                :: hfile
    verb = .false.; if(present(verbose)) verb = verbose
    ext = shape(covmat)
    call open_hdf_file(fname, hfile, "w")
    call write_hdf(hfile, "nside",    nside)
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    call write_hdf(hfile, "ordering", order)
    if(verb) write(*,'(a10,i5)') "order: ", order
    call write_hdf(hfile, "inverse", inv_)
    if(verb) write(*,'(a10,l)') "inverse: ", inv_
    if(verb) write(*,'(a10,4i5)') "shape: ", ext
    call write_hdf(hfile, "cov",      covmat)
    call write_hdf(hfile, "pixels",   pixels)
    call close_hdf_file(hfile)
    if(present(inv)) inv = inv_ /= 0
  end subroutine

  subroutine read_covmat_2_unf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
    implicit none
    character(len=*)              :: fname
    real(dp), allocatable         :: covmat(:,:)
    integer(i4b), allocatable     :: pixels(:)
    logical(lgt),     optional    :: inv, verbose
    integer(i4b) :: nside, order, i, j, ncomp, ntot, n, unit
    logical(lgt) :: inverse, verb

    verb = .false.; if(present(verbose)) verb = verbose
    unit = getlun()
    if(allocated(covmat)) deallocate(covmat)
    if(allocated(pixels)) deallocate(pixels)
    open(unit, file=fname, form='unformatted', action="read", status="old")
	
        read(unit) ntot
    if(verb) write(*,'(a10,i5)') "ntot: ", ntot
    read(unit) order
    if(verb) write(*,'(a10,i5)') "order: ", order

    
    read(unit) ncomp
    if(verb) write(*,'(a10,i5)') "ncomp: ", ncomp
    allocate(covmat(ntot,ntot), pixels(ntot/ncomp))
    do i = 1, ntot
       read(unit) covmat(:,i)
    end do
    read(unit) inverse
    if(verb) write(*,'(a10,l)') "inverse: ", inverse
    read(unit) nside
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    read(unit) pixels
    close(unit)
    if(present(inv)) inv = inverse
  end subroutine

  subroutine read_covmat_4_unf(covmat, pixels, nside, order,  fname, inv, verbose)
    implicit none
    character(len=*)              :: fname
    real(dp), allocatable         :: covmat(:,:,:,:)
    integer(i4b), allocatable     :: pixels(:)
    logical(lgt),     optional    :: inv, verbose
    integer(i4b) :: nside, order, i, j, ncomp, ntot, n, unit
    logical(lgt) :: inverse, verb

    verb = .false.; if(present(verbose)) verb = verbose
    unit = getlun()
    if(allocated(covmat)) deallocate(covmat)
    if(allocated(pixels)) deallocate(pixels)
    open(unit, file=fname, form='unformatted', action="read", status="old")
	
        read(unit) ntot
    if(verb) write(*,'(a10,i5)') "ntot: ", ntot
    read(unit) order
    if(verb) write(*,'(a10,i5)') "order: ", order

    
    read(unit) ncomp
    if(verb) write(*,'(a10,i5)') "ncomp: ", ncomp
    n = ntot/ncomp
    allocate(covmat(n,ncomp,n,ncomp), pixels(n))
    do j = 1, ncomp
       do i = 1, n
          read(unit) covmat(:,:,i,j)
       end do
    end do
    read(unit) inverse
    if(verb) write(*,'(a10,l)') "inverse: ", inverse
    read(unit) nside
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    read(unit) pixels
    close(unit)
    if(present(inv)) inv = inverse
  end subroutine

  subroutine write_covmat_2_unf(covmat, pixels, nside, order, ncomp, fname, inv, verbose)
    implicit none
    character(len=*)              :: fname
    real(dp)         :: covmat(:,:)
    integer(i4b)     :: pixels(:)
    logical(lgt),     optional    :: inv, verbose
    integer(i4b) :: nside, order, i, j, ncomp, ntot, n, unit
    logical(lgt) :: inverse, verb

    verb = .false.; if(present(verbose)) verb = verbose
    unit = getlun()
    open(unit, file=fname, form='unformatted')
	n = size(covmat,1)
    ntot = size(covmat,1)
        write(unit) ntot
    if(verb) write(*,'(a10,i5)') "ntot: ", ntot
    write(unit) order
    if(verb) write(*,'(a10,i5)') "order: ", order

    
    write(unit) ncomp
    if(verb) write(*,'(a10,i5)') "ncomp: ", ncomp
    do i = 1, ntot
       write(unit) covmat(:,i)
    end do
    inverse = .false.; if(present(inv)) inverse = inv
    write(unit) inverse
    if(verb) write(*,'(a10,l)') "inverse: ", inverse
    write(unit) nside
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    write(unit) pixels
    close(unit)
    
  end subroutine

  subroutine write_covmat_4_unf(covmat, pixels, nside, order,  fname, inv, verbose)
    implicit none
    character(len=*)              :: fname
    real(dp)         :: covmat(:,:,:,:)
    integer(i4b)     :: pixels(:)
    logical(lgt),     optional    :: inv, verbose
    integer(i4b) :: nside, order, i, j, ncomp, ntot, n, unit
    logical(lgt) :: inverse, verb

    verb = .false.; if(present(verbose)) verb = verbose
    unit = getlun()
    open(unit, file=fname, form='unformatted')
	n = size(covmat,1)
    ntot = size(covmat,1)*size(covmat,2)
    write(unit) ntot
    if(verb) write(*,'(a10,i5)') "ntot: ", ntot
    write(unit) order
    if(verb) write(*,'(a10,i5)') "order: ", order

    
    ncomp = size(covmat,2)
    write(unit) ncomp
    do j = 1, ncomp
       do i = 1, n
          write(unit) covmat(:,:,i,j)
       end do
    end do
    inverse = .false.; if(present(inv)) inverse = inv
    write(unit) inverse
    if(verb) write(*,'(a10,l)') "inverse: ", inverse
    write(unit) nside
    if(verb) write(*,'(a10,i5)') "nside: ", nside
    write(unit) pixels
    close(unit)
    
  end subroutine


  function get_extension(filename, override) result(ext)
    implicit none
    character(len=*)                  :: filename
    character(len=*),    optional     :: override
    character(len=len_trim(filename)) :: ext
    integer(i4b)                      :: n, i
    if(present(override)) then
       ext = override
       return
    end if
    n = len_trim(filename)
    do i = n,1,-1
       if(filename(i:i) == ".") exit
    end do
    ext = filename(i+1:n)
  end function

end module

