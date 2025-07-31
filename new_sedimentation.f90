 subroutine sedimentation(this)

        class(prec_parcel_type), intent(inout) :: this
        double precision :: asr1,asr2,slope
        double precision :: D


        integer :: n

        !$omp parallel do default(shared) private(n,D)
        do n = 1, this%local_num

                !The D* implemtation in case we want to keep it
                !D = ((rho_air/rho_w)*(this%qr(n)/this%nr(n)))**(f13)
                !this%delta_pos(this%z_dim, n) = this%delta_pos(this%z_dim, n) - (a1*(D**(b1))*(exp(-f1*D)))+a2*(D**(b2)) &
                !&* (exp(-f2*D))*(rho_ref/rho_air)**(f12)

                ! This is the slope parameter of the distribution
                slope = ((pi/6)*(ro_0/ro_air)*(this%nr(n)/this%qr(n))*(shape+1)*(shape+2)*(shape+3))**((f13))

                !These are the mass weighted integrals for abel and shipway terminal velocity
                asr1 = a1*((ro_0/ro_air)**(f12))*(slope**(1+shape)*(slope+f1)**(-(1+shape+b1))) &
                *(gamma(1+shape+b1)/gamma(1+shape)) 
                asr2 = a2*((ro_0/ro_air)**(f12))*(slope**(1+shape)*(slope+f2)**(-(1+shape+b2))) &
                *(gamma(1+shape+b2)/gamma(1+shape)) 
                this%delta_pos = asr1 + asr2
        end do
        !$omp end parallel do

    end subroutine sedimentation