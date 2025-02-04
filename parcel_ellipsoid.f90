 module parcel_ellipsoid

    use parcel_container
    implicit none

    ! Adding the ellipsoid geomerty to the dynamics
    type, extends(dynamic_parcel_type) :: ellipsoid_parcel_type ! add procedures for ellipsoid below
        double precision, allocatable, dimension(:,:) :: B
        double precision, allocatable, dimension(:,:) :: delta_B
        double precision, allocatable, dimension(:,:) :: strain
        double precision, allocatable, dimension(:,:) :: Vetas
        double precision, allocatable, dimension(:,:) :: Vtaus
        character(len=16)   :: shape_type ! (e.g. "ellipsoid5", for B with 5 elements, strain with 8)

        contains
            procedure :: alloc => ellipsoid_parcel_alloc
            procedure :: dealloc => ellipsoid_parcel_dealloc
            procedure :: resize => ellipsoid_parcel_resize

            ! Other ellipsoid procedures to go here
    end type

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_alloc(this, num)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%dynamic_parcel_type%alloc(num)

            if (this%shape_type == "ellipsoid5") then
                allocate(this%B(5, num))
                allocate(this%delta_B(5, num))
                allocate(this%strain(8, num))
                allocate(this%Vetas(3, num))
                allocate(this%Vtaus(3, num))

                call this%register_attribute(this%B(1, :), "B11", "m^2")
                call this%register_attribute(this%B(2, :), "B12", "m^2")
                call this%register_attribute(this%B(3, :), "B13", "m^2")
                call this%register_attribute(this%B(4, :), "B22", "m^2")
                call this%register_attribute(this%B(5, :), "B23", "m^2")
                call this%register_attribute(this%delta_B(1, :), "B11_rk_tendency", "m^2/s")
                call this%register_attribute(this%delta_B(2, :), "B12_rk_tendency", "m^2/s")
                call this%register_attribute(this%delta_B(3, :), "B13_rk_tendency", "m^2/s")
                call this%register_attribute(this%delta_B(4, :), "B22_rk_tendency", "m^2/s")
                call this%register_attribute(this%delta_B(5, :), "B23_rk_tendency", "m^2/s")
                call this%register_attribute(this%strain(1, :), "DUDX", "1/s")
                call this%register_attribute(this%strain(2, :), "DUDY", "1/s")
                call this%register_attribute(this%strain(3, :), "DUDZ", "1/s")
                call this%register_attribute(this%strain(4, :), "DVDX", "1/s")
                call this%register_attribute(this%strain(5, :), "DVDY", "1/s")
                call this%register_attribute(this%strain(6, :), "DVDZ", "1/s")
                call this%register_attribute(this%strain(7, :), "DWDX", "1/s")
                call this%register_attribute(this%strain(8, :), "DWDY", "1/s")
                call this%register_attribute(this%Vetas(1, :), "Veta1", "m")
                call this%register_attribute(this%Vetas(2, :), "Veta2", "m")
                call this%register_attribute(this%Vetas(3, :), "Veta3", "m")
                call this%register_attribute(this%Vtaus(1, :), "Vtau1", "m")
                call this%register_attribute(this%Vtaus(2, :), "Vtau2", "m")
                call this%register_attribute(this%Vtaus(3, :), "Vtau3", "m")
            else
                print *, "Ellipsoid type not known."
                stop
            endif

        end subroutine ellipsoid_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_dealloc(this)
            class(ellipsoid_parcel_type), intent(inout) :: this

            call try_deallocate(this%B)
            call try_deallocate(this%delta_B)
            call try_deallocate(this%strain)
            call try_deallocate(this%Vetas)
            call try_deallocate(this%Vtaus)

            call this%dynamic_parcel_type%dealloc

        end subroutine ellipsoid_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_resize(this, new_size)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%dynamic_parcel_type%resize(new_size)

            call resize_array(this%B, new_size, this%local_num)
            call resize_array(this%delta_B, new_size, this%local_num)
            call resize_array(this%strain, new_size, this%local_num)
            call resize_array(this%Vetas, new_size, this%local_num)
            call resize_array(this%Vtaus, new_size, this%local_num)

            if (this%shape_type == "ellipsoid5") then
                call this%reset_attribute(this%B(1, :), "B11")
                call this%reset_attribute(this%B(2, :), "B12")
                call this%reset_attribute(this%B(3, :), "B13")
                call this%reset_attribute(this%B(4, :), "B22")
                call this%reset_attribute(this%B(5, :), "B23")
                call this%reset_attribute(this%delta_B(1, :), "B11_rk_tendency")
                call this%reset_attribute(this%delta_B(2, :), "B12_rk_tendency")
                call this%reset_attribute(this%delta_B(3, :), "B13_rk_tendency")
                call this%reset_attribute(this%delta_B(4, :), "B22_rk_tendency")
                call this%reset_attribute(this%delta_B(5, :), "B23_rk_tendency")
                call this%reset_attribute(this%strain(1, :), "DUDX")
                call this%reset_attribute(this%strain(2, :), "DUDY")
                call this%reset_attribute(this%strain(3, :), "DUDZ")
                call this%reset_attribute(this%strain(4, :), "DVDX")
                call this%reset_attribute(this%strain(5, :), "DVDY")
                call this%reset_attribute(this%strain(6, :), "DVDZ")
                call this%reset_attribute(this%strain(7, :), "DWDX")
                call this%reset_attribute(this%strain(8, :), "DWDY")
                call this%reset_attribute(this%Vetas(1, :), "Veta1")
                call this%reset_attribute(this%Vetas(2, :), "Veta2")
                call this%reset_attribute(this%Vetas(3, :), "Veta3")
                call this%reset_attribute(this%Vtaus(1, :), "Vtau1")
                call this%reset_attribute(this%Vtaus(2, :), "Vtau2")
                call this%reset_attribute(this%Vtaus(3, :), "Vtau3")
            else
                print *, "Ellipsoid type not known."
                stop
            endif

        end subroutine ellipsoid_parcel_resize

end module
