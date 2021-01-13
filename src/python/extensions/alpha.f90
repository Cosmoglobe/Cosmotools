function fortran_getalpha(phi,theta,psi) result(alpha)
	use quick_qpoint
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

function fortran_getazimuth(phi,theta,time) result(az)
  use quick_qpoint
  implicit none
  
  real(dp) :: phi,theta,az,lst,el,time
  real(dp), dimension(3) :: vec
  real(dp), dimension(3,3) :: rotmat

  call ang2vec(theta, phi, vec)
  lst    = QP_mjd2lst(time, CBI_GEODETIC_LONGITUDE)
  rotmat = QP_matrix_equ2hor(lst, CBI_GEODETIC_LATITUDE)
  vec    = matmul(rotmat, vec)
  call vec2ang(vec, el, az)

end function fortran_getazimuth
