module inner
    use kind_type 
    use global 
    use nr, only: plgndr_s
    implicit none
    complex(CDP), save, allocatable, private, protected :: inner_a(:)
contains 
    

subroutine inner_coeff(l)
    use math_const, only: i => math_i 
    integer(I4B), intent(in) :: l 
    real   (RSP) :: ka
    real   (RDP) :: tmp2 
    complex(CDP) :: tmp1 
    integer(I4B) :: j

    ka   = (2.d0*Mass*Kinet)**0.5d0*ra
    tmp1 = A(l)*(bessel_j(l, ka) -K(l)*bessel_y(l, ka))
    do j = 1, N 
        tmp2 = H(j, N)/(2.d0*Mass*(E(j) -Kinet))
        inner_a(j) = tmp1*tmp2/R(l)
    end do     
end subroutine inner_coeff










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_inner_achive(l)
    integer(I4B), intent(in) :: l 
    complex(CQP) :: sum 
    integer(I4B) :: i, j

    allocate(inner_a(1:N))
    call inner_coeff(l) 
    do i = 1, N 
        sum = 0.d0 
        do j = 1, N 
            sum = sum +inner_a(j)*H(j, i)
        end do 
        inner_u(l, i) = sum 
        write(*, *) sum 
    end do 
    deallocate(inner_a)
end subroutine PROC_inner_achive


subroutine PROC_inner_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    use math_const,  only: pi => math_pi, degree => math_degree 
    use hamiltonian, only: coord_r, coord_theta
    integer (I1B), parameter :: file_psi1 = 101, file_psi2 = 102, file_psi3 = 103
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real    (RDP), parameter :: radian_to_degree = 1.d0/degree 
    real    (RSP) :: tmp 
    complex (CQP) :: sum 
    integer (I4B) :: i, j, k 

    open(file_psi1, file = "inout/inner_u0.d")
    sum = 0.d0 
    do i = 1, N
        sum = inner_u(0, i)
        write(file_psi1, form_psi) coord_r(i), dble(abs(sum)**2.d0)
    end do 
    close(file_psi1)

    open(file_psi2, file = "inout/inner_psi.d")
    open(file_psi3, file = "inout/inner_psi-4pir.d")
    do i = 1, N, N/pr 
        do j = 0, ptheta
            sum = 0.d0 
            do k = 0, L 
                tmp = cos(coord_theta(j))
                sum = sum +inner_u(k, i)/coord_r(i)*plgndr_s(i, 0, tmp)
            end do 
!             write(file_psi2, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(sum)**2.d0)
            write(file_psi2, form_psi) coord_r(i), coord_theta(j), dble(abs(sum)**2.d0)
            sum = sum*4.d0*pi*coord_r(i)
!             write(file_psi3, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(sum)**2.d0)
            write(file_psi3, form_psi) coord_r(i), coord_theta(j), dble(abs(sum)**2.d0)
        end do 
        write(file_psi2, form_psi) 
        write(file_psi3, form_psi) 
    end do 
    close(file_psi2)
    close(file_psi3)
end subroutine PROC_inner_plot
end module inner
