 module parcel_types

    use parcel_ellipsoid
    use params
    use constants
    use grid_container
    use utils
    implicit none

    type, extends(ellipsoid_parcel_type) :: idealised_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: humidity
        double precision, allocatable, dimension(:) :: buoyancy

        contains
            procedure :: alloc => idealised_parcel_alloc
            procedure :: dealloc=> idealised_parcel_dealloc
            procedure :: resize => idealised_parcel_resize

            ! get_buoyancy added here
    end type

    type, extends(ellipsoid_parcel_type) :: realistic_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: qv
        double precision, allocatable, dimension(:) :: ql
        double precision, allocatable, dimension(:) :: theta
        double precision, allocatable, dimension(:) :: Nl ! optional droplet number
        logical :: has_droplets = .false.

        contains
            procedure :: alloc => realistic_parcel_alloc
            procedure :: dealloc=> realistic_parcel_dealloc
            procedure :: resize => realistic_parcel_resize

            ! get_buoyancy added here
    end type

    type, extends(base_parcel_type) :: prec_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: volume
        double precision, allocatable, dimension(:) :: qr
        double precision, allocatable, dimension(:) :: Nr ! droplet number
        double precision, allocatable, dimension(:) :: vterm
        double precision, allocatable, dimension(:) :: prevp
        double precision, allocatable, dimension(:) :: nrevp


        contains
            procedure :: alloc => prec_parcel_alloc
            procedure :: dealloc=> prec_parcel_dealloc
            procedure :: resize => prec_parcel_resize
            procedure :: sedimentation
            procedure ::  evaporation
            procedure goners
            

            ! get_buoyancy added here
    end type

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_alloc(this, num)
            class(idealised_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%ellipsoid_parcel_type%alloc(num)

            allocate(this%humidity(num))
            allocate(this%buoyancy(num))

            call this%register_attribute(this%humidity, "humidity", "kg/kg")
            call this%register_attribute(this%buoyancy, "buoyancy", "m/s^2")

        end subroutine idealised_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_dealloc(this)
            class(idealised_parcel_type), intent(inout) :: this

            call try_deallocate(this%humidity)
            call try_deallocate(this%buoyancy)
            call this%ellipsoid_parcel_type%dealloc

        end subroutine idealised_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_resize(this, new_size)
            class(idealised_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%ellipsoid_parcel_type%resize(new_size)

            call resize_array(this%humidity, new_size, this%local_num)
            call resize_array(this%buoyancy, new_size, this%local_num)

            call this%reset_attribute(this%humidity, "humidity")
            call this%reset_attribute(this%buoyancy, "buoyancy")

        end subroutine idealised_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_alloc(this, num)
            class(realistic_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%ellipsoid_parcel_type%alloc(num)

            allocate(this%theta(num))
            allocate(this%qv(num))
            allocate(this%ql(num))

            call this%register_attribute(this%theta, "theta", "K")
            call this%register_attribute(this%qv, "qv", "kg/kg")
            call this%register_attribute(this%ql, "ql", "kg/kg")

            if(this%has_droplets) then
               allocate(this%Nl(num))
               call this%register_attribute(this%Nl, "Nl", "/m^3")
            endif

        end subroutine realistic_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_dealloc(this)
            class(realistic_parcel_type), intent(inout) :: this

            call try_deallocate(this%theta)
            call try_deallocate(this%qv)
            call try_deallocate(this%ql)
            call try_deallocate(this%Nl)

            call this%ellipsoid_parcel_type%dealloc

        end subroutine realistic_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_resize(this, new_size)
            class(realistic_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%ellipsoid_parcel_type%resize(new_size)

            call resize_array(this%theta, new_size, this%local_num)
            call resize_array(this%qv, new_size, this%local_num)
            call resize_array(this%ql, new_size, this%local_num)

            call this%reset_attribute(this%theta, "theta")
            call this%reset_attribute(this%qv, "qv")
            call this%reset_attribute(this%ql, "ql")

            if(this%has_droplets) then
                call resize_array(this%Nl, new_size, this%local_num)
                call this%reset_attribute(this%Nl, "Nl")
            endif

        end subroutine realistic_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_alloc(this, num)
            class(prec_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%base_alloc(num)

            allocate(this%volume(num))
            allocate(this%qr(num))
            allocate(this%nr(num))
            allocate(this%vterm(num))
            allocate(this%prevp(num))
            allocate(this%nrevp(num))            

            call this%register_attribute(this%volume, "volume", "m^3")
            call this%register_attribute(this%qr, "qr", "kg/kg")
            call this%register_attribute(this%nr, "nr", "/m^3")
            call this%register_attribute(this%vterm, "vterm", "ms^-1")
            call this%register_attribute(this%prevp, "prevp", "kg/kg")
            call this%register_attribute(this%nrevp, "nrevp", "/m^3")

        end subroutine prec_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_dealloc(this)
            class(prec_parcel_type), intent(inout) :: this

            call try_deallocate(this%volume)
            call try_deallocate(this%qr)
            call try_deallocate(this%nr)
            call try_deallocate(this%vterm)
            call try_deallocate(this%prevp)
            call try_deallocate(this%nrevp)


            call this%base_dealloc

        end subroutine prec_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_resize(this, new_size)
            class(prec_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%base_resize(new_size)

            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%qr, new_size, this%local_num)
            call resize_array(this%nr, new_size, this%local_num)
            call resize_array(this%vterm, new_size, this%local_num)
            call resize_array(this%prevp, new_size, this%local_num)
            call resize_array(this%nrevp, new_size, this%local_num)

            call this%reset_attribute(this%volume, "volume")
            call this%reset_attribute(this%qr, "qr")
            call this%reset_attribute(this%nr, "nr")
            call this%reset_attribute(this%vterm, "vterm")
            call this%reset_attribute(this%prevp, "prevp")
            call this%reset_attribute(this%nrevp, "nrevp")
           

        end subroutine prec_parcel_resize

       


    subroutine sedimentation(this,mesh)

        class(prec_parcel_type), intent(inout) :: this
        class(Grid), intent(in) :: mesh
        double precision :: D, press, exn,theta, qv, temp, ro_air
        integer :: n

            parcel_loop: do n = 1, this%local_num

                !Here I believe we will call trilinear and par2grid to get 
                !the theta and qv values we need to make the parcel fall
                !This Could be wrapped up in grid2par

            !! ----------------------grid2par---------------------------
                call mesh%grid2par(position=this%position(:, n), theta=theta, qv=qv)
                press=surf_press*exp(-this%position(1, n)/pressure_scale_height)
                exn=(press/ref_press)**(r_d/c_p)
                temp=theta*exn
                ro_air = press/(r_d*temp)
                !! ----------------------------------------------------
                
                !! Approximate diameter scaling from Rooney 2025
                D = ((ro_air/ro_w)*(this%qr(n)/this%nr(n)))**(f13)
                
               
                this%vterm(n) =  (a1*(D**(b1))*(exp(-f1*D)))+a2*(D**(b2)) &
                &* (exp(-f2*D))*(ro_0/ro_air)**(f12)
               
             
                
                D = 0.0
                press = 0.0
                exn = 0.0
                theta = 0.0
                qv = 0.0
                temp = 0.0
                ro_air = 0.0




            end do parcel_loop

    end subroutine sedimentation

    subroutine evaporation(this,mesh)
        class(prec_parcel_type), intent(inout) :: this
        class(Grid), intent(inout) :: mesh
        double precision :: press, exn,theta, qv, temp, ro_air, slope, vent_r,abliq,ws
        double precision :: evap_mass, evap_heat
        integer :: n



        parcel_loop: do n =1, this%local_num
                        
                !Here I believe we will call trilinear and par2grid to get 
                !the theta and qv values we need to make the parcel fall
                !This Could be wrapped up in grid2par
                call mesh%grid2par(position=this%position(:, n), theta=theta, qv=qv)
                press=surf_press*exp(-this%position(1, n)/pressure_scale_height)
                exn=(press/ref_press)**(r_d/c_p)
                
                temp=theta*exn
                ro_air = press/(r_d*temp)
                ws = 3.8/(0.01*press*e**(-17.2693882*(temp-273.15)/(temp-35.86))-6.109)

                !Calculating slope coefficient 
                slope = ((pi/6)*(ro_r/ro_air)*(this%nr(n)/this%qr(n))*(shape+1)*(shape+2)*(shape+3))**((f13))
                
                vent_r = two*pi*(this%nr(n))*ro_air* &
                    (0.78*((one+shape)/(slope)) &
                    +  0.31*((a1*ro_air/visc)**(f12))*(Sc**(f13))*((ro_0/ro_air)**(f14))  &
                    * (gamma((f12*b1 +shape +f52))/gamma((1+shape))) &
                    *((one + (f12*f1)/slope)**(-(f12*b1 + shape + f52))) &
                    *((slope)**(-f12*b1 -f32)))

                !Thermodynamic coefficient
                abliq = ((l_v**2)/(k_a*r_v* (temp**2))) + (1/(ro_air*ws*diffus))
                
                !Sink term for rainwater mass mixing ratio due to evaporation
                evap_mass = -(((qv/ws)-1)/(ro_air*ABliq))*vent_r
                evap_heat = evap_mass*(l_v/(c_p*exn))
                !call mesh%par2grid(position=this%position(:,n),evap_mass=evap_mass,evap_heat=evap_heat)
                if (evaporation_off .eqv. .true.) then
                    this%prevp = 0
                else 
                    this%prevp(n) = (((qv/ws)-1)/(ro_air*ABliq))*vent_r
                end if 
                !Sink term for rainwater number concentration due to evaporation
                !!Here I have added a switch for rain on rain aggregation 
                if (aggregation_on .eqv. .true.) then 
                    this%nrevp(n) = this%prevp(n)*(((this%nr(n))/ro_air)/(this%qr(n))) - 8*err*this%qr(n)*(this%nr(n))
                else 
                    this%nrevp(n) = this%prevp(n)*(((this%nr(n))/ro_air)/(this%qr(n)))
                end if 
                
        end do parcel_loop 

    end subroutine evaporation
    
    subroutine goners(this,rainfall)
        class(prec_parcel_type), intent(inout) :: this
        double precision, intent(out) :: rainfall
        integer, allocatable :: pid(:)  ! Declare pid as an allocatable array
        integer :: n_del
        integer :: n
    
        n_del = 0
        allocate(pid(this%local_num))  ! Allocate pid with the size of local_num
    
        do n = 1, this%local_num
            if (this%position(1, n) <= 0) then
                n_del = n_del + 1
                pid(n_del) = n
                print *, "Ping! Parcel impinged"
                !! Check me I may be wrong
                rainfall = rainfall + 522*((ro_0/ro_r)*(this%qr(n)/this%nr(n))**(3*f13))
                
                cycle
            else if (this%qr(n) <= 0) then
                n_del = n_del + 1
                pid(n_del) = n
                print *, "Poof! Parcel evaporated"
                cycle
            else if (this%nr(n) <= 0) then
                n_del = n_del + 1
                pid(n_del) = n
                print *, "Poof! Parcel evaporated"
                cycle
            end if
        end do
        
        if (n_del > 0) then
            pid = pid(1:n_del)  ! Resize pid to only include the first n_del elements
            call this%pdelete(pid=pid, n_del=n_del)
        end if
        deallocate(pid)  ! Deallocate pid to free memory
    end subroutine goners


end module
