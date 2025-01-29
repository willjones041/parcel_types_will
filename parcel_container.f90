module parcel_container
    use datatypes, only : int64
    use armanip, only : resize_array
    implicit none

    type scalar_ptr
        double precision, pointer :: sptr(:)
    end type scalar_ptr

    type vector_ptr
        double precision, pointer :: vptr(:,:)
    end type vector_ptr

    type, abstract :: position_type
        double precision, allocatable, dimension(:,:) :: position
        integer             :: attr_num     ! number of parcel attributes for serialisation
        integer             :: local_num    ! local number of parcels
        integer(kind=int64) :: total_num    ! global number of parcels (over all MPI ranks)
        integer             :: max_num      ! capacity per attribute, i.e. maximum number of parcels

        contains
            procedure :: position_allocate
            procedure :: position_deallocate
            !procedure :: position_serialize
            !procedure :: position_deserialize
            ! Basic procedures common to all derived classes:
            !procedure :: pack                 => parcel_pack
            !procedure :: unpack               => parcel_unpack
            !procedure :: delete               => parcel_delete
            ! Pure virtual procedures:
            ! procedure(parcel_allocate),       deferred :: allocate
            ! procedure(parcel_deallocate),     deferred :: deallocate
            ! procedure(parcel_serialize),      deferred :: serialize
            ! procedure(parcel_deserialize),    deferred :: deserialize
    end type

    type, extends(position_type) :: base_parcel_type
        type(scalar_ptr), allocatable, dimension(:) :: scalar_attribs
        type(vector_ptr), allocatable, dimension(:) :: vector_attribs
        character(len=16), allocatable, dimension(:) :: scalar_names !for reverse lookup
        character(len=16), allocatable, dimension(:) :: vector_names !for reverse lookup
        character(len=16), allocatable, dimension(:) :: scalar_units !for reverse lookup
        character(len=16), allocatable, dimension(:) :: vector_units !for reverse lookup
        integer :: n_scalars
        integer :: n_vectors

        contains
            procedure :: base_parcel_allocate
            procedure :: base_parcel_deallocate
            procedure :: update_attr_num
            procedure :: print_me
    end type

    type, extends(base_parcel_type) :: idealised_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: humidity
        double precision, allocatable, dimension(:) :: buoyancy

        contains
            procedure :: idealised_parcel_allocate
            procedure :: idealised_parcel_deallocate
    end type

    type, extends(base_parcel_type) :: realistic_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: qv
        double precision, allocatable, dimension(:) :: ql
        double precision, allocatable, dimension(:) :: theta

        contains
            procedure :: realistic_parcel_allocate
            procedure :: realistic_parcel_deallocate
    end type

    type, extends(base_parcel_type) :: prec_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: qr
        double precision, allocatable, dimension(:) :: Nr

        contains
            procedure :: prec_parcel_allocate
            procedure :: prec_parcel_deallocate
    end type

    ! Lower type
    type, extends(realistic_parcel_type) :: realistic_ellipsoid_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: volume
        double precision, allocatable, dimension(:,:) :: B

        contains
           procedure          :: realistic_ellipsoid_parcel_allocate
           procedure          :: realistic_ellipsoid_parcel_deallocate
    end type

    ! Lower type
    type, extends(idealised_parcel_type) :: idealised_ellipsoid_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: volume
        double precision, allocatable, dimension(:,:) :: B

        contains
           procedure          :: idealised_ellipsoid_parcel_allocate
           procedure          :: idealised_ellipsoid_parcel_deallocate
    end type

    ! Lower type
    type, extends(prec_parcel_type) ::prec_sphere_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: volume
        ! procedure to add the volume
        contains
            procedure          :: prec_sphere_parcel_allocate
            procedure          :: prec_sphere_parcel_deallocate
    end type

    contains


       !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


        subroutine update_attr_num(this)
            class(base_parcel_type), intent(inout), target :: this
            integer :: j
            this%attr_num = this%n_scalars


            do j = 1, this%n_vectors
                this%attr_num = this%attr_num + size(this%vector_attribs(j)%vptr, dim=1)
            end do
        end subroutine

        subroutine print_me(this)
            class(base_parcel_type), intent(inout), target :: this
            integer :: j

            do j = 1, this%n_scalars
                write(*,*) this%scalar_names(j)
                write(*,*) size(this%scalar_attribs(j)%sptr, dim=1)
            end do


            do j = 1, this%n_vectors
                write(*,*) this%vector_names(j)
                write(*,*) size(this%vector_attribs(j)%vptr, dim=1)
                write(*,*) size(this%vector_attribs(j)%vptr, dim=2)
            end do

            write(*,*) "this%attr_num"
            write(*,*) this%attr_num
        end subroutine


        ! Allocate parcel memory
        ! ATTENTION: Extended types must allocate additional parcel attributes
        !            in their own routine.
        ! @param[in] num number of parcels
        ! @param[in] n_pos number of spatial dimensions
        subroutine position_allocate(this, num, n_pos)
            class(position_type), intent(inout) :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos

            allocate(this%position(n_pos, num))
            this%max_num = num

        end subroutine position_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate parcel memory
        ! ATTENTION: Extended types must deallocate additional parcel attributes
        !            in their own routine.
        subroutine position_deallocate(this)
            class(position_type), intent(inout) :: this

            if (.not. allocated(this%position)) then
                return
            endif

            this%local_num = 0
            this%total_num = 0
            this%max_num   = 0

            deallocate(this%position)

        end subroutine position_deallocate

        subroutine base_parcel_allocate(this, num, n_pos)
            class(base_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer :: n_vectors_temp

            n_vectors_temp = this%n_vectors
            this%n_vectors = this%n_vectors + 1 ! Bookkeeping before allocation

            call this%position_allocate(num, n_pos)
            allocate(this%scalar_attribs(this%n_scalars))
            allocate(this%vector_attribs(this%n_vectors))
            allocate(this%scalar_names(this%n_scalars))
            allocate(this%vector_names(this%n_vectors))
            this%vector_names(n_vectors_temp+1) = 'position'
            this%vector_attribs(n_vectors_temp+1)%vptr => this%position

        end subroutine base_parcel_allocate

        subroutine base_parcel_deallocate(this)
            class(base_parcel_type), intent(inout) :: this
            call this%position_deallocate()
            deallocate(this%scalar_attribs)
            deallocate(this%vector_attribs)

        end subroutine base_parcel_deallocate

        subroutine idealised_parcel_allocate(this, num, n_pos)
            class(idealised_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer :: n_scalars_temp

            n_scalars_temp = this%n_scalars
            this%n_scalars=this%n_scalars + 2 ! Bookkeeping before allocation
            call this%base_parcel_allocate(num, n_pos)
            allocate(this%humidity(num))
            allocate(this%buoyancy(num))
            this%scalar_names(n_scalars_temp+1) = 'humidity'
            this%scalar_attribs(n_scalars_temp)%sptr => this%humidity
            this%scalar_names(n_scalars_temp+2) = 'buoyancy'
            this%scalar_attribs(n_scalars_temp)%sptr => this%buoyancy

        end subroutine idealised_parcel_allocate

        subroutine idealised_parcel_deallocate(this)
            class(idealised_parcel_type), intent(inout) :: this
            call this%base_parcel_deallocate()
            deallocate(this%humidity)
            deallocate(this%buoyancy)
        end subroutine idealised_parcel_deallocate

        subroutine realistic_parcel_allocate(this, num, n_pos)
            class(realistic_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer :: n_scalars_temp

            n_scalars_temp = this%n_scalars
            this%n_scalars=this%n_scalars + 3 ! Bookkeeping before allocation
            call this%base_parcel_allocate(num, n_pos)
            allocate(this%qv(num))
            allocate(this%ql(num))
            allocate(this%theta(num))
            this%scalar_names(n_scalars_temp+1) = 'qv'
            this%scalar_attribs(n_scalars_temp+1)%sptr => this%qv
            this%scalar_names(n_scalars_temp+2) = 'ql'
            this%scalar_attribs(n_scalars_temp+2)%sptr => this%ql
            this%scalar_names(n_scalars_temp+3) = 'theta'
            this%scalar_attribs(n_scalars_temp+3)%sptr => this%theta
        end subroutine realistic_parcel_allocate

        subroutine realistic_parcel_deallocate(this)
            class(realistic_parcel_type), intent(inout) :: this
            call this%base_parcel_deallocate()
            deallocate(this%qv)
            deallocate(this%ql)
            deallocate(this%theta)
        end subroutine realistic_parcel_deallocate

        subroutine prec_parcel_allocate(this, num, n_pos)
            class(prec_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer :: n_scalars_temp

            n_scalars_temp = this%n_scalars
            this%n_scalars = this%n_scalars + 2 ! Bookkeeping before allocation
            call this%base_parcel_allocate(num, n_pos)
            allocate(this%qr(num))
            allocate(this%Nr(num))
            this%scalar_names(n_scalars_temp+1) = 'qr'
            this%scalar_attribs(n_scalars_temp+1)%sptr => this%qr
            this%scalar_names(n_scalars_temp+2) = 'Nr'
            this%scalar_attribs(n_scalars_temp+2)%sptr => this%Nr

        end subroutine prec_parcel_allocate

        subroutine prec_parcel_deallocate(this)
            class(prec_parcel_type), intent(inout) :: this
            call this%base_parcel_deallocate()
            deallocate(this%qr)
            deallocate(this%Nr)
        end subroutine prec_parcel_deallocate

        subroutine prec_sphere_parcel_allocate(this, num, n_pos)
            class(prec_sphere_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos

            this%n_scalars= 1 ! Bookkeeping before allocation
            this%n_vectors= 0
            call this%prec_parcel_allocate(num, n_pos)
            allocate(this%volume(num))
            this%scalar_names(1) = 'volume'
            this%scalar_attribs(1)%sptr => this%volume
            call this%update_attr_num

        end subroutine prec_sphere_parcel_allocate

        subroutine prec_sphere_parcel_deallocate(this)
            class(prec_sphere_parcel_type), intent(inout) :: this
            call this%prec_parcel_deallocate()
            deallocate(this%volume)
            this%n_scalars=0 ! Reset base_parcel_type_properties
            this%n_vectors=0 ! Reset base_parcel_type_properties
        end subroutine prec_sphere_parcel_deallocate

        subroutine realistic_ellipsoid_parcel_allocate(this, num, n_pos, n_shape)
            class(realistic_ellipsoid_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer,        intent(in)    :: n_shape

            this%n_scalars= 1 ! Bookkeeping before allocation
            this%n_vectors= 1 ! Bookkeeping before allocation
            call this%realistic_parcel_allocate(num, n_pos)
            allocate(this%volume(num))
            allocate(this%B(n_shape, num))
            this%scalar_names(1) = 'volume'
            this%scalar_attribs(1)%sptr => this%volume
            this%vector_names(1) = 'B'
            this%vector_attribs(1)%vptr => this%B
            call this%update_attr_num

        end subroutine realistic_ellipsoid_parcel_allocate

        subroutine realistic_ellipsoid_parcel_deallocate(this)
            class(realistic_ellipsoid_parcel_type), intent(inout) :: this
            call this%realistic_parcel_deallocate()
            deallocate(this%volume)
            deallocate(this%B)
            this%n_scalars=0 ! Reset base_parcel_type_properties
            this%n_vectors=0 ! Reset base_parcel_type_properties
        end subroutine realistic_ellipsoid_parcel_deallocate

        subroutine idealised_ellipsoid_parcel_allocate(this, num, n_pos, n_shape)
            class(idealised_ellipsoid_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos
            integer,        intent(in)    :: n_shape

            this%n_scalars= 1 ! Bookkeeping before allocation
            this%n_vectors= 1 ! Bookkeeping before allocation
            call this%idealised_parcel_allocate(num, n_pos)
            allocate(this%volume(num))
            allocate(this%B(n_shape, num))
            this%scalar_names(1) = 'volume'
            this%scalar_attribs(1)%sptr => this%volume
            this%vector_names(1) = 'B'
            this%vector_attribs(1)%vptr => this%B
            call this%update_attr_num

        end subroutine idealised_ellipsoid_parcel_allocate

        subroutine idealised_ellipsoid_parcel_deallocate(this)
            class(idealised_ellipsoid_parcel_type), intent(inout) :: this
            call this%idealised_parcel_deallocate()
            deallocate(this%volume)
            deallocate(this%B)
            this%n_scalars=0 ! Reset base_parcel_type_properties
            this%n_vectors=0 ! Reset base_parcel_type_properties
        end subroutine idealised_ellipsoid_parcel_deallocate

end module
