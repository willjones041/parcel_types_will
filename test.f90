program test
use parcel_container
implicit none

    type(realistic_ellipsoid_parcel_type) :: parcels
    type(prec_sphere_parcel_type) :: prec_parcels

call parcels%realistic_ellipsoid_parcel_allocate(4, 3, 5)
call prec_parcels%prec_sphere_parcel_allocate(4, 3)

parcels%volume=1.0
parcels%ql=2.0
parcels%qv=3.0
parcels%theta=4.0
parcels%position=5.0
parcels%B=5.0

prec_parcels%volume=1.0
prec_parcels%qr=2.0
prec_parcels%Nr=3.0

call parcels%print_me
call prec_parcels%print_me

call parcels%realistic_ellipsoid_parcel_deallocate
call prec_parcels%prec_sphere_parcel_deallocate

end program

