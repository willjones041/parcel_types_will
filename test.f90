program test
use parcel_types
use grid_container
use params
use utils
implicit none
    

        type(prec_parcel_type) :: prec_parcels
        type(Grid) :: testgrid
        
        double precision, allocatable, dimension(:,:,:) :: theta_init, qv_init
        integer ::  tlim, t, delt, n,i
        double precision :: rainfall
        
        t = 0
        delt = 1
        tlim = 20*60

        

        ! Allocate and initialize the grid
        call testgrid%alloc()

        ! Allocate and set initial values for theta and qv
        allocate(theta_init(nz, ny, nx))
        allocate(qv_init(nz, ny, nx))
        !---------------------------------------------------------------------------
        call set_atmos(xi=nx,yi=ny,zi=nz,theta_array = theta_init,qv_array=qv_init)
        !--------------------------This needs to be a call to a subroutine ---------
        ! Set the fields in the grid
        call testgrid%set_fields(theta_init, qv_init)
       
        ! Print the grid values for verification
        call testgrid%print_me()


        call prec_parcels%set_dimension(3)
        call prec_parcels%alloc(100)
       

        ! Set the initial position of the parcels
        
        ! Set the parcel microphysical properties
        prec_parcels%qr = 0.01
        prec_parcels%nr = 1000
        
        prec_parcels%position(1,:) = [(200 + (i-1)*100, i=1, 90)]
        prec_parcels%position(2,:) = extent(2)/2
        prec_parcels%position(3,:) = extent(2)/2
        rainfall = 0

        do while (t <= tlim)
            call prec_parcels%evaporation(testgrid)
            call prec_parcels%sedimentation(testgrid)
            
            parcel_loop: do n =1, prec_parcels%local_num
            
            
            prec_parcels%position(1,n) = prec_parcels%position(1,n) - prec_parcels%vterm(n) * delt
            prec_parcels%qr(n) = prec_parcels%qr(n) + prec_parcels%prevp(n) * delt
            prec_parcels%nr(n) = prec_parcels%nr(n) + prec_parcels%nrevp(n) * delt
            
            print *, "time = ",t,"position =", prec_parcels%position(1,n), "Nr =", prec_parcels%nr(n), "qr =", prec_parcels%qr(n), &
             "vterm =", prec_parcels%vterm(n), "prevp =", prec_parcels%prevp(n), "nrevp =", prec_parcels%nrevp(n)
            
            
            
        
            end do parcel_loop

            !Removing impinged and evaporated parcels
            call prec_parcels%goners(rainfall)
            print *, "Rainfall rate = ", (rainfall), "kg"
            t = t+delt
        end do 

        call prec_parcels%dealloc
    
end program

