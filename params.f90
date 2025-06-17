module params
    use constants
    implicit none

    ! Parameters
    !Abel and shipway fall speed constants
    double precision, parameter :: a1 = 4854.0
    double precision, parameter :: a2 = 446.0
    double precision, parameter :: b1 = 1.0
    double precision, parameter :: b2 = 0.782
    double precision, parameter :: f1 = 195
    double precision, parameter :: f2 = 4085.35
    !Euler's number
    double precision, parameter :: e = 2.718281828459045
    !SL air density
    double precision, parameter :: ro_0 = 1.2256
    !Density of liquid water
    double precision, parameter :: ro_w = 1000
    
    !Grid inputs to find dx and 
    double precision,  parameter ::lower(3) = [0,0,0]
    double precision,   parameter :: extent(3) = [5000,12800,12800]
    integer, parameter :: nx = 5, ny =5, nz=5
    double precision,   parameter:: dx(3) = extent/dble((/nz,ny,nx/))
    double precision, parameter :: dxi(3) = one/dx
    !Parameters for the atmosphere
    double precision, parameter :: pressure_scale_height = 8619 !From Bomex, please change (placeholder)
    double precision, parameter :: surf_press = 100000
    double precision, parameter :: ref_press = 100000
   
    ! gas constant of dry air
    double precision, protected :: r_d=287.05
    ![J/(kg*K)] specific heat at constant pressure and dry air:
    double precision, protected :: c_p = 1005.7d0
    ![J/kg] latent heat of vaporization:
    double precision, protected :: l_v = 2.501e6
    
    !diffusivity of water on air
    double precision, parameter :: diffus = 2.42e-5
    !kinematic vscocity of air at sea level
    double precision, parameter :: visc = 1.14e-5
    double precision, protected :: sc = visc/diffus
    !rainwater density
    double precision, parameter :: ro_r = 1000.0
    !Thermal conductivity of air 
    double precision, parameter :: k_a =0.02623
    !Gas constant for water vapour
    double precision, parameter :: r_v = 461.52
    !shape parameter for the rainwater DSD
    double precision, parameter :: shape = 2.5
    ! Switch to turn on rain on rain self collection within parcel DSDs
    logical, parameter :: aggregation_on = .true.
    !Option to switch off evaporation
    logical, parameter :: evaporation_off = .false.
    ! Collection efficiency rain on rain
    double precision, parameter :: err = 1.0


end module params