program test
use parcel_types
use grid_container
implicit none
    type(Grid) :: testgrid
    type(realistic_parcel_type) :: parcels
    type(prec_parcel_type) :: prec_parcels
    integer :: i
    double precision :: ro_air



call testgrid%alloc()
call testgrid%fill()


call parcels%set_dimension(3)
parcels%has_droplets=.true.
parcels%has_labels=.true.
parcels%shape_type="ellipsoid5"
call prec_parcels%set_dimension(3)

call parcels%alloc(4)
call prec_parcels%alloc(4)

parcels%volume=1.0
parcels%ql=2.0
parcels%qv=3.0
parcels%theta=4.0
parcels%position=5.0
parcels%B=5.0
parcels%Nl=0.0
parcels%label=0

prec_parcels%volume=1.0


! call parcels%print_me
! call prec_parcels%print_me

call parcels%resize(7)
call prec_parcels%resize(7)

do  i = 1, size(prec_parcels%qr)

    prec_parcels%qr(i)=i
    prec_parcels%Nr(i)=i*2
end do
! call parcels%print_me
! call prec_parcels%print_me

!I have set up a use case of updating the z_position deltas 
!based off the terminal velocity equation

ro_air = 1.11
!Testing the fall subroutine
call prec_parcels%fall(ro_air)

! do i = 1, size(prec_parcels%delta_pos(3,:))
!     print *, prec_parcels%delta_pos(3,i)
! end do

!call prec_parcels%print_me

!call parcels%dealloc
call prec_parcels%dealloc




end program

