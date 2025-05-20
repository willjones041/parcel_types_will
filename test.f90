program test
use parcel_types
use grid_container
use params
implicit none
    

        type(prec_parcel_type) :: prec_parcels
        type(Grid) :: testgrid

        double precision, allocatable, dimension(:,:,:) :: theta_init, qv_init
        double precision :: RH = 0.8
        integer :: i, j, k, tlim, t, delt, n
        
        t = 0
        delt = 1.0
        tlim = 10

        

        ! Allocate and initialize the grid
        call testgrid%alloc()

        ! Allocate and set initial values for theta and qv
        allocate(theta_init(nx, ny, nz))
        allocate(qv_init(nx, ny, nz))
        theta_init = 300.0
        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    qv_init(i,j,k) = RH * exp(17.27 * (theta_init(i,j,k) - 273.15) / (theta_init(i,j,k) - 35.86))
                end do
            end do
        end do
        

        ! Set the fields in the grid
        call testgrid%set_fields(theta_init, qv_init)

        ! Print the grid values for verification
        ! call testgrid%print_me()


        call prec_parcels%set_dimension(3)
        call prec_parcels%alloc(1)

        ! Set the initial position of the parcels
        
        ! Set the parcel microphysical properties
        prec_parcels%qr = 0.001
        prec_parcels%nr = 100000
        prec_parcels%position(1,:) = 50
        

        do while (t <= tlim)
            call prec_parcels%evaporation(testgrid)
            call prec_parcels%sedimentation(testgrid)
            
            parcel_loop: do n =1, prec_parcels%local_num
            
            
            prec_parcels%position(1,n) = prec_parcels%position(1,n) + prec_parcels%vterm(n) * delt
            prec_parcels%qr(n) = prec_parcels%qr(n) + prec_parcels%prevp(n) * delt
            prec_parcels%nr(n) = prec_parcels%nr(n) + prec_parcels%nrevp(n) * delt
            
            
            end do parcel_loop
            print *, "time = ",t,"position =", prec_parcels%position(1,1), "Nr =", prec_parcels%nr(n), "qr =", prec_parcels%qr(n), &
            "vterm =", prec_parcels%vterm(n), "prevp =", prec_parcels%prevp(n), "nrevp =", prec_parcels%nrevp(n)
            
            t = t+delt
        end do 

        call prec_parcels%dealloc
    
end program

