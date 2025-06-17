module utils
  use params
    implicit none
    contains

      ! Tri-linear interpolation
        ! @param[in] pos position vector
        ! @param[out] ii zonal lower grid point for interoplation
        ! @param[out] jj meridional lower grid point for interpolation
        ! @param[out] kk vertical lower grid point for interpolation
        ! @param[out] ww interpolation weights
      !This util should only do one parcel at a time,
      !grid2par and oar2grid is where we will loop over the parcel arrays
    pure subroutine trilinear(pos,ii, jj, kk, ww)
    double precision, intent(in)  :: pos(3)

    integer,          intent(out) :: ii, jj, kk
    double precision, intent(out) :: ww(0:1, 0:1, 0:1)
    double precision              :: xyz(3)
    double precision              :: w00, w10, w01, w11
    double precision              :: px, py, pz, pxc, pyc, pzc


    ! (i, j, k)
    xyz = (pos - lower) * dxi
    ii = floor(xyz(1))
    jj = floor(xyz(2))
    kk = floor(xyz(3))

    px = xyz(1) - dble(ii)
    pxc = 1 - px

    py = xyz(2) - dble(jj)
    pyc = 1 - py

    pz = xyz(3) - dble(kk)
    pzc = 1 - pz

    w00 = pyc * pxc
    w10 = pyc * px
    w01 = py * pxc
    w11 = py * px

    ! Note order of indices is k,j,i
    ww(0, 0, 0) = pzc * w00
    ww(0, 0, 1) = pzc * w10
    ww(0, 1, 0) = pzc * w01
    ww(0, 1, 1) = pzc * w11
    ww(1, 0, 0) = pz * w00
    ww(1, 0, 1) = pz * w10
    ww(1, 1, 0) = pz * w01
    ww(1, 1, 1) = pz * w11

end subroutine trilinear



pure subroutine set_atmos(xi,yi,zi,theta_array,qv_array)

integer, intent(in) :: xi,yi,zi
double precision, intent(out),dimension(:,:,:) :: theta_array, qv_array
double precision :: epsilon,  RH, theta,press,exn,temp,ws
integer :: k,j,i
epsilon = 0.622
      ! reference pressure in Pa (1000 hPa)     ! Rd / cp
RH = 0.7            ! example: 80% relative humidity
theta = 300
do i = 1, xi
  do j = 1, yi
    do k = 1, zi
    press = surf_press*exp(-zi*dx(3)/pressure_scale_height)
    exn=(press/ref_press)**(r_d/c_p)
    temp=theta*exn
    ws = 3.8/(press*e**(-17.2693882*(temp-273.15)/(temp-35.86))-6.109)
    theta_array(k,j,i) = theta
    qv_array(k,j,i) = ws*RH
    end do
  end do
end do


end subroutine 
end module