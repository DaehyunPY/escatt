module inner
    use kind_type 
    use global 
    implicit none
    complex(dp),       allocatable, private, protected :: inner_a(:)
    complex(dp), save, allocatable, private, protected :: inner_u(:, :)
contains 
    

subroutine inner_coeff(l)
    use math_const, only: i => math_i 
    use gsl_special, only: gsl_sf_bessel_jsl, gsl_sf_bessel_ysl
    integer(i4), intent(in) :: l 
    real   (dp) :: ka, sb_j, sb_y, tmp2 
    complex(dp) :: tmp1 
    integer(i4) :: j

    ka = (2.d0*Scatt)**0.5d0*ra
    sb_j = gsl_sf_bessel_jsl(l, ka)
    sb_y = gsl_sf_bessel_ysl(l, ka)
    
    tmp1 = A(l)*(sb_j -K(l)*sb_y)
    do j = 1, N 
        tmp2       = H(j, N)/(2.d0*(E(j) -Scatt))
        inner_a(j) = tmp1*tmp2/R(l)
    end do     
end subroutine inner_coeff









! ==================================================
! PROCESS
! ==================================================
! inner achive -------------------------------------
subroutine PROC_inner_achive(l)
    integer(i4), intent(in) :: l 
    complex(qp) :: sum 
    integer(i4) :: i, j

    if(allocated(inner_u) == .false.) allocate(inner_u(0:L, 1:N))
    if(allocated(inner_a) == .false.) allocate(inner_a(1:N))
    call inner_coeff(l) 
    do i = 1, N 
        sum = 0.d0 
        do j = 1, N 
            sum = sum +inner_a(j)*H(j, i)
        end do 
        inner_u(l, i) = sum 
    end do 
    if(allocated(inner_a) == .true.) deallocate(inner_a)
end subroutine PROC_inner_achive
! end inner achive ---------------------------------
! inner plot ---------------------------------------
subroutine PROC_inner_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use hamiltonian, only: coord_r, coord_theta
    use gsl_special, only: gsl_sf_legendre_Pl
    integer  (i1), parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.d0/degree 
    real     (dp) :: tmp, unit_theta
    complex  (qp) :: sum 
    integer  (i4) :: i, j, k 

    unit_theta = 1_dp 
    if(op_degree == "Y") unit_theta = radian_to_degree

    open(file_psi1, file = "output/inner_u_0.d")
    sum = 0.d0 
    do i = 1, N
        sum = inner_u(0, i)
        write(file_psi1, form_psi) coord_r(i), dble(abs(sum)**2.d0)
    end do 
    close(file_psi1)

    open(file_psi2, file = "output/inner_psi.d")
    do i = 1, N, N/pr 
        do j = 0, ptheta
            sum = 0.d0 
            do k = 0, L 
                tmp = cos(coord_theta(j))
                sum = sum +inner_u(k, i)/coord_r(i)*gsl_sf_legendre_Pl(k, tmp)
            end do 
            write(file_psi2, form_psi) coord_r(i), coord_theta(j)*unit_theta, dble(abs(sum)**2.d0)
        end do 
        write(file_psi2, form_psi) 
    end do 
    if(allocated(inner_u) == .true.) deallocate(inner_u)
    close(file_psi2)
end subroutine PROC_inner_plot
! end inner plot -----------------------------------
end module inner
