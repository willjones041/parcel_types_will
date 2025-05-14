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

end module