module grid_container
    use params
    use utils
    implicit none

    type :: Grid
        private
        !mesh spacing
        double precision :: dx(3)
        !inverse meshspacing
        double precision:: dxi
        !Declare specific humidity field
        double precision, allocatable, dimension(:,:,:) :: qv
        !Declare potential temperature field
        double precision, allocatable, dimension(:,:,:) :: theta
        !Declare liquid water mixing ratio
        double precision, allocatable, dimension(:,:,:) :: ql
        contains
            procedure :: alloc => prec_grid_alloc
            procedure :: set_fields => set_prec_fields
            procedure :: print_me => print_me_grid
            procedure :: get_attribs => get_grid_attribs
            procedure :: set_attribs => set_grid_attribs
            procedure :: grid2par
            procedure :: par2grid

       end type Grid
    
contains

        subroutine prec_grid_alloc(this)
        !!Here is where the user input is taken to allocate the fields
        ! and initialise the fixed parameters
            class(Grid), intent(inout) :: this
            this%dx = extent/dble([nz,ny,nx])
            
            allocate(this%qv(nz,ny,nx))
            allocate(this%theta(nz,ny,nz))

        end subroutine prec_grid_alloc

        subroutine set_prec_fields(this, theta_val, qv_val)
            !! Set the values of theta and qv fields
            class(Grid), intent(inout) :: this
            double precision, intent(in), dimension(:,:,:) :: theta_val, qv_val
    
            if (size(theta_val) /= size(this%theta)) then
                print*, "Error: theta_val size does not match grid dimensions"
                stop
            end if
    
            if (size(qv_val) /= size(this%qv)) then
                print*, "Error: qv_val size does not match grid dimensions"
                stop
            end if
    
            this%theta = theta_val
            this%qv = qv_val
        end subroutine set_prec_fields

        subroutine print_me_grid(this)
        !! This subroutine prints the gridded values 
        !! So I can check them
            class(Grid),intent(in) :: this
            print*, this%theta
            print*, this%qv
        end subroutine

        subroutine get_grid_attribs(this,ii,jj,kk,thetag_val,qvg_val) 
            !!This is the getter function for theta and qv
            !! It takes in the indexes of the grid box and outputs the values for each
            class(Grid), intent(in) :: this
            double precision, intent(out),dimension(:,:,:) :: thetag_val, qvg_val
            integer, intent(in) :: ii,jj,kk

            thetag_val = this%theta(ii:ii+1,jj:jj+1,kk:kk+1)
            qvg_val = this%qv(ii:ii+1,jj:jj+1,kk:kk+1)
            
        end subroutine get_grid_attribs

        subroutine set_grid_attribs(this, ii, jj, kk, thetag_val, qvg_val)
            !! This is the setter function for theta and qv
            !! It takes in the indexes of the grid box and updates the values for each
            class(Grid), intent(inout) :: this
            double precision, intent(in), dimension(:,:,:) :: thetag_val, qvg_val
            integer, intent(in) :: ii, jj, kk
        
            ! Update the grid's theta and qv subarrays
            this%theta(ii:ii+1, jj:jj+1, kk:kk+1) = thetag_val
            this%qv(ii:ii+1, jj:jj+1, kk:kk+1) = qvg_val
        end subroutine set_grid_attribs
        
        
        
        subroutine grid2par(this, position, theta, qv)
            class(Grid), intent(in) :: this
            double precision, intent(in) :: position(3)
            double precision, intent(out) :: theta, qv
            double precision, dimension(0:1, 0:1, 0:1) :: weights, theta_subarray, qv_subarray
            integer :: is, js, ks
            ! Call trilinear to compute weights and indices
            call trilinear(pos=position, ii=is, jj=js, kk=ks, ww=weights)
            ! Retrieve subarray values from the grid
            call this%get_attribs(ii=is, jj=js, kk=ks, thetag_val=theta_subarray, qvg_val=qv_subarray)
            
            ! Compute theta and qv using weights
            theta = sum(theta_subarray * weights)
            qv = sum(qv_subarray * weights)
        end subroutine grid2par


        subroutine par2grid(this, position, evap_mass,evap_heat)
            class(Grid), intent(inout) :: this
            double precision, intent(in) :: position(3)
            double precision, intent(in) :: evap_mass,evap_heat
            double precision, dimension(0:1, 0:1, 0:1) :: weights
            integer :: is, js, ks
            double precision, dimension(0:1, 0:1, 0:1) :: theta_subarray, qv_subarray
        
            ! Call trilinear to compute weights and indices
            call trilinear(pos=position, ii=is, jj=js, kk=ks, ww=weights)
        
            ! Retrieve the current subarray values from the grid
            call this%get_attribs(ii=is, jj=js, kk=ks, thetag_val=theta_subarray, qvg_val=qv_subarray)
        
            ! Update the subarray values using weights
            theta_subarray = theta_subarray + weights *evap_heat 
            qv_subarray = qv_subarray + weights *evap_mass
        
            ! Write the updated values back to the grid
            call this%set_attribs(ii=is, jj=js, kk=ks, thetag_val=theta_subarray, qvg_val=qv_subarray)
        end subroutine par2grid

end module grid_container