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



subroutine set_atmos(xi,yi,zi,theta_array,qv_array,RH,theta)

integer, intent(in) :: xi,yi,zi
double precision, intent(out),dimension(:,:,:) :: theta_array, qv_array
double precision, dimension(zi,yi,xi) :: press_array, temp_array, ws_array
double precision, intent(in) :: RH, theta
double precision :: epsilon,press,exn,temp,ws
integer :: k,j,i
epsilon = 0.622
      ! reference pressure in Pa (1000 hPa)     ! Rd / cp

do i = 1, xi
  do j = 1, yi
    do k = 1, zi
    press = surf_press*exp(-k*dx(3)/pressure_scale_height)
    exn=(press/ref_press)**(r_d/c_p)
    temp=theta*exn
    ws = 3.8/((0.01*press)*exp(-17.2693882*(temp-273.15)/(temp-35.86))-6.109)
    theta_array(k,j,i) = theta
    qv_array(k,j,i) = ws*RH
    press_array(k,j,i) = press
    temp_array(k,j,i) = temp
    ws_array(k,j,i) = ws
    end do
  end do
end do
! Plot the k component of qv_array against k*dx(3)
open(unit=10, file='atmos_profile.dat', status='replace')
write(10,*) '# Height(m)   Pressure(Pa)   Temperature(K)   Saturation_specific_humidity(kg/kg)  &
 &Potential_temperature(K)   Specific_humdity(kg/kg)   Relative_humidity'
do k = 1, zi
  write(10,*) k*dx(3), press_array(k,1,1), temp_array(k,1,1), ws_array(k,1,1), &
  theta_array(k,1,1), qv_array(k,1,1), qv_array(k,1,1)/ws_array(k,1,1)
end do
close(10)

end subroutine 
end module