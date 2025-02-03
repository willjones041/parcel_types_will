program test
use parcel_container
implicit none

    type(realistic_parcel_type) :: parcels
    type(prec_parcel_type) :: prec_parcels

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
prec_parcels%qr=2.0
prec_parcels%Nr=3.0

call parcels%print_me
call prec_parcels%print_me

call parcels%dealloc
call prec_parcels%dealloc

end program

