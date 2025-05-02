 module parcel_types

    use parcel_ellipsoid
    use params
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

        contains
            procedure :: alloc => prec_parcel_alloc
            procedure :: dealloc=> prec_parcel_dealloc
            procedure :: resize => prec_parcel_resize
            procedure :: fall => prec_parcel_fall
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
            allocate(this%Nr(num))

            call this%register_attribute(this%volume, "volume", "m^3")
            call this%register_attribute(this%qr, "qr", "kg/kg")
            call this%register_attribute(this%Nr, "Nr", "/m^3")

        end subroutine prec_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_dealloc(this)
            class(prec_parcel_type), intent(inout) :: this

            call try_deallocate(this%volume)
            call try_deallocate(this%qr)
            call try_deallocate(this%Nr)

            call this%base_dealloc

        end subroutine prec_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_resize(this, new_size)
            class(prec_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%base_resize(new_size)

            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%qr, new_size, this%local_num)
            call resize_array(this%Nr, new_size, this%local_num)

            call this%reset_attribute(this%volume, "volume")
            call this%reset_attribute(this%qr, "qr")
            call this%reset_attribute(this%Nr, "Nr")

        end subroutine prec_parcel_resize

    subroutine prec_parcel_fall(this,ro_air)

        class(prec_parcel_type), intent(inout) :: this
        double precision :: ro_air, D
        integer :: i 
        
        do i = 1, size(this%delta_pos(3,:))
            D = ((ro_air/ro_w)*(this%qr(i)/this%Nr(i)))
            this%delta_pos(3,i) =  (a1*(D**(b1))*(e**(-f1*D)))+a2*(D**(b2))*(e**(-f2*D))*(ro_0/ro_air)**(0.5)
        end do 

    end subroutine

    
end module
