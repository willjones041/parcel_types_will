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
            procedure :: fill => fill_grid
            procedure :: print_me => print_me_grid
            procedure :: get_attribs
       end type Grid
    
contains

        subroutine prec_grid_alloc(this)
        !!Here is where the user input is taken to allocate the fields
        ! and initialise the fixed parameters
            class(Grid), intent(inout) :: this
            this%dx = extent/dble([nx,ny,nz])
            
            allocate(this%qv(nx,ny,nz))
            allocate(this%theta(nx,ny,nz))

        end subroutine prec_grid_alloc

        subroutine fill_grid(this)
            !This is a placeholder subroutine that allows me to fill
            ! the grid without needing a netcdf4 file
            !REPLACE WITH READ NETCDF SUBROUTINE
            class(Grid), intent(inout) :: this
            integer :: i,j,k

            zloop: do k = 1, nx
                yloop: do j = 1, ny
                    xloop: do i = 1, nz
                    !Setting the isothermal case up first
                                this%theta(k,j,i) = 300
                                this%qv(k,j,i) = 0.001
                        end do xloop
                    end do yloop
                end do zloop
        end subroutine

        subroutine print_me_grid(this)
        !! This subroutine prints the gridded values 
        !! So I can check them
            class(Grid),intent(in) :: this
            print*, this%theta
        end subroutine 

        subroutine get_attribs(this,ii,jj,kk,thetag_val,qvg_val) 
            !!This is the getter function for theta and qv
            !! It takes in the indexes of the grid box and outputs the values for each
            class(Grid), intent(in) :: this
            double precision, intent(out),dimension(:,:,:) :: thetag_val, qvg_val
            integer, intent(in) :: ii,jj,kk

            thetag_val = this%theta(kk:kk+1,jj:jj+1,ii:ii+1)
            qvg_val = this%qv(kk:kk+1,jj:jj+1,ii:ii+1)
        end subroutine get_attribs


end module grid_container