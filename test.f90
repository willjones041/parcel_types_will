program test
use parcel_container
implicit none

    type(realistic_ellipsoid_parcel_type) :: parcels
    type(prec_sphere_parcel_type) :: prec_parcels
    integer :: j
call parcels%realistic_ellipsoid_parcel_allocate(4, 3, 5)
call prec_parcels%prec_sphere_parcel_allocate(4, 3)

parcels%volume=1.0
parcels%ql=2.0
parcels%qv=3.0
parcels%theta=4.0
parcels%position=5.0
parcels%B=5.0

do j = 1, parcels%n_scalars
    write(*,*) parcels%scalar_names(j)
    write(*,*) size(parcels%scalar_attribs(j)%sptr, dim=1)
end do


do j = 1, parcels%n_vectors
    write(*,*) parcels%vector_names(j)
    write(*,*) size(parcels%vector_attribs(j)%vptr, dim=1)
    write(*,*) size(parcels%vector_attribs(j)%vptr, dim=2)
end do

call parcels%realistic_ellipsoid_parcel_deallocate
call prec_parcels%prec_sphere_parcel_deallocate

end program

