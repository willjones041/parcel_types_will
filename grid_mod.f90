module grid_container
    use params
    
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
        contains
            procedure :: alloc => prec_grid_alloc
            procedure :: set_fields => set_prec_fields
            procedure :: print_me => print_me_grid
            procedure :: get_attribs => get_prec_attribs
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
            print*, this%theta(1,:,:)
            print*, this%qv(1,:,:)
        end subroutine 

        subroutine get_prec_attribs(this,ii,jj,kk,thetag_val,qvg_val) 
            !!This is the getter function for theta and qv
            !! It takes in the indexes of the grid box and outputs the values for each
            class(Grid), intent(in) :: this
            double precision, intent(out),dimension(:,:,:) :: thetag_val, qvg_val
            integer, intent(in) :: ii,jj,kk

            thetag_val = this%theta(ii:ii+1,jj:jj+1,kk:kk+1)
            qvg_val = this%qv(ii:ii+1,jj:jj+1,kk:kk+1)

        end subroutine get_prec_attribs


end module grid_container