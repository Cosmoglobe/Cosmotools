module quick_qpoint
	implicit none
	integer, parameter :: dp=8
	real(dp), parameter :: PI = 3.1415926535897931
	real(dp), parameter :: TWOPI = 6.2831853071795862
	real(dp), parameter :: PI_BY_4 = 0.78539816339744828d0
  	real(dp), parameter :: QP_DTOR  = 1.7453292519943295769d-2 ! Conversion from degrees to radians
  	real(dp), parameter :: QP_RTOD  = 1.d0 / 1.7453292519943295769d-2 ! Conversion from radians to degrees
	real(dp), parameter :: CBI_GEODETIC_LONGITUDE  = -67.76166667d0  ! 67 degrees, 45 minutes, 42.0 seconds WEST
	real(dp), parameter :: CBI_GEODETIC_LATITUDE   =  -23.02822222d0 ! 23 degrees, 1 minute, 41.6 seconds SOUTHcontains

	contains

function fortran_getazimuth(phi,theta,time) result(az)
	real(dp) :: phi,theta,az,lst,el,time
	real(dp), dimension(3) :: vec
	real(dp), dimension(3,3) :: rotmat

	call ang2vec(theta, phi, vec)
	lst    = QP_mjd2lst(time, CBI_GEODETIC_LONGITUDE)
	rotmat = QP_matrix_equ2hor(lst, CBI_GEODETIC_LATITUDE)
	vec    = matmul(rotmat, vec)
	call vec2ang(vec, el, az)

end function fortran_getazimuth

function fortran_getalpha(phi,theta,psi) result(alpha)

	implicit none
!	logical galactic
	real(dp):: phi,theta,psi,alpha
	real(dp), dimension(3):: p0
    real(dp), dimension(3,1) :: s0
	real(dp), dimension(3,3) :: rotmat_equ2gal
!	write(*,*) phi,theta,psi
	!rotmat_equ2gal  = QP_matrix_equ2gal()

	call QP_p0_equ(QP_RTOD*phi,QP_RTOD*theta, p0)
	call QP_s0_equ(QP_RTOD*phi,QP_RTOD*theta,QP_RTOD*psi, s0)
	alpha = PI_BY_4 - QP_DTOR * QP_Psi(p0,s0)
end function fortran_getalpha


!subroutine fortran_getalphaarray(n,phi,theta,psi,alpha)
!	real(dp), dimension(0:n-1) :: phi,theta,psi,alpha
!	
!	do i=0,n-1
!		alpha(i) = fortran_getalpha(phi(i),theta(i),psi(i),.false.)
!	enddo
!end subroutine fortran_getalphaarray


subroutine QP_p0_equ(ra, dec, vec) 
  implicit none
  
  real(dp), intent(in) :: ra, dec
  real(dp), dimension(3), intent(out) :: vec
  
  vec(1) = cos(QP_DTOR * ra) * cos(QP_DTOR * dec)
  vec(2) = sin(QP_DTOR * ra) * cos(QP_DTOR * dec)
  vec(3) = sin(QP_DTOR * dec)
  
end subroutine QP_p0_equ

subroutine QP_s0_equ(ra, dec, par, vec)
  implicit none
  
  real(dp),               intent(in)  :: ra, dec, par
  real(dp), dimension(3), intent(out) :: vec
  
  vec(1) = -1.d0 * cos(QP_DTOR * ra) * sin(QP_DTOR * dec) * cos(QP_DTOR * par) - &
       & sin(QP_DTOR * ra) * sin(QP_DTOR * par)
  vec(2) = -1.d0 * sin(QP_DTOR * ra) * sin(QP_DTOR * dec) * cos(QP_DTOR * par) + &
       & cos(QP_DTOR * ra) * sin(QP_DTOR * par)
  vec(3) = cos(QP_DTOR * dec) * cos(QP_DTOR * par)
  
end subroutine QP_s0_equ



function QP_psi(p_equ, s_equ) 
  implicit none

  real(dp), dimension(3), intent(in) :: p_equ, s_equ
  real(dp)                           :: QP_psi

  real(dp) :: psi

  psi = 180.d0 - QP_par(p_equ, s_equ)

  if (psi < 180.d0) then
     QP_psi = psi
  else
     QP_psi = psi - 360.d0
  end if

end function QP_psi

function QP_par(p_equ, s_equ) 
  implicit none

  real(dp), dimension(3), intent(in) :: p_equ, s_equ
  real(dp)                           :: QP_par

  real(dp) :: par
  real(dp), dimension(3) :: z

  z    = 0.d0
  z(3) = 1.d0    

  par = atan2(dot_product(s_equ, crossprod(z, p_equ)), &
       & dot_product(s_equ, crossprod(p_equ, crossprod(z, p_equ)))) / QP_DTOR

  if (par == 180.d0) then
     QP_par = -180.d0
  else
     QP_par =  par
  end if

end function QP_par

function crossprod(u, v) RESULT(w)
  real(dp), dimension(3), intent(in) :: u, v
  real(dp), dimension(3)             :: w
  
  w(1) = u(2) * v(3) - u(3) * v(2)
  w(2) = u(3) * v(1) - u(1) * v(3)
  w(3) = u(1) * v(2) - u(2) * v(1)
end function crossprod

function QP_mjd2lst(mjd, longitude)
  implicit none

  real(dp), intent(in) :: mjd, longitude
  real(dp)             :: QP_mjd2lst

  real(dp) :: s, du, tu, gmst, lst

  s  = 86400.d0 * (mjd - floor(mjd))   ! seconds since midnight UTC
  du = mjd - 51544.5d0                 ! days since J2000 epoch
  tu = du / 36525.d0                   ! convert du to centuries

  ! Greenwich Mean Sidereal Time
  ! Formula from Astrophysical Formulae, Volume 2 (3Ed) by Kenneth Lang, 1999.
  ! Conforms to IAU 1976 System of Astronomical Constants, 1980 IAU Theory of Nutation, FK5 catalog
  gmst = s + 24110.54841d0 + tu * (8640184.812866d0 + tu * (0.093104d0 - tu * 0.0000062d0)); 

  gmst = gmst * (360.d0 / 86400.d0);

  ! add telescope longitude, fix to the range [0,360)
  lst = gmst + longitude
  do while (lst < 0.d0)    
     lst = lst + 360.d0
  end do
  do while (lst >= 360.d0) 
     lst = lst - 360.d0
  end do

  QP_mjd2lst = lst

end function QP_mjd2lst


subroutine ang2vec(theta, phi, vector)
  !=======================================================================
  !     renders the vector (x,y,z) corresponding to angles
  !     theta (co-latitude measured from North pole, in [0,Pi] radians)
  !     and phi (longitude measured eastward, in radians)
  !     North pole is (x,y,z)=(0,0,1)
  !     added by EH, Feb 2000
  !=======================================================================
  REAL(dp), INTENT(IN) :: theta, phi
  REAL(dp), INTENT(OUT), dimension(1:) :: vector

  REAL(dp) :: sintheta
  !=======================================================================

  if (theta<0.0_dp .or. theta>pi)  then
     stop "ANG2VEC: theta is out of range [0, Pi]"
  endif
  sintheta = SIN(theta)

  vector(1) = sintheta * COS(phi)
  vector(2) = sintheta * SIN(phi)
  vector(3) = COS(theta)

  return
end subroutine ang2vec
!=======================================================================
subroutine vec2ang(vector, theta, phi)
  !=======================================================================
  !     renders the angles theta, phi corresponding to vector (x,y,z)
  !     theta (co-latitude measured from North pole, in [0,Pi] radians)
  !     and phi (longitude measured eastward, in [0,2Pi[ radians)
  !     North pole is (x,y,z)=(0,0,1)
  !     added by EH, Feb 2000
  !=======================================================================
  REAL(dp), INTENT(IN), dimension(1:) :: vector
  REAL(dp), INTENT(OUT) :: theta, phi

  REAL(dp) :: dnorm, z
  !=======================================================================

  dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)

  z = vector(3) / dnorm
  theta = ACOS(z)

  phi = 0.0_dp
  if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
       &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
  if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

  return
end subroutine vec2ang

function QP_matrix_hor2equ(lst, lat) result(matrix)
  implicit none

  real(dp), intent(in)     :: lst, lat
  real(dp), dimension(3,3) :: matrix

  call compute_euler_matrix_zyz(QP_DTOR * (lst - 180.d0), QP_DTOR * (lat - 90.d0), 0.d0, matrix)
  
end function QP_matrix_hor2equ

function QP_matrix_equ2hor(lst, lat) result(matrix)
  implicit none

  real(dp), intent(in)     :: lst, lat
  real(dp), dimension(3,3) :: matrix

  matrix = QP_matrix_hor2equ(lst, lat)
  matrix = transpose(matrix)
  
end function QP_matrix_equ2hor

! Convention: First psi around z, then theta around y, then phi around z
subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
  implicit none
                                                                              
  real(dp),                 intent(in)  :: phi, theta, psi
  real(dp), dimension(3,3), intent(out) :: euler_matrix

  real(dp) :: sphi, cphi, sth, cth, spsi, cpsi
   
  sphi = sin(phi)
  cphi = cos(phi)
   
  sth  = sin(theta)
  cth  = cos(theta)

  spsi = sin(psi)
  cpsi = cos(psi)

  euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
  euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
  euler_matrix(1,3) =                sth * cphi
  euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
  euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
  euler_matrix(2,3) =                sth * sphi
  euler_matrix(3,1) =              - sth * cpsi
  euler_matrix(3,2) =                sth * spsi
  euler_matrix(3,3) =                cth

end subroutine compute_euler_matrix_zyz


end module
