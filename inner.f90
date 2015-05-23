module inner
    use kind_type 
    use global 
    implicit none
    complex(DP), save, allocatable, private, protected :: inner_a(:)
contains 
    

subroutine inner_coeff(l)
    use math_const, only: i => math_i 
    use nr, only: sphbes_s
    integer(I4), intent(in) :: l 
    real   (SP) :: ka, sb_j, sb_y, diff_j, diff_y 
    real   (DP) :: tmp2 
    complex(DP) :: tmp1 
    integer(I4) :: j

    ka = (2.d0*Mass*Kinet)**0.5d0*ra
    call sphbes_s(l, ka, sb_j, sb_y, diff_j, diff_y)
    
    tmp1 = A(l)*(sb_j -K(l)*sb_y)
    do j = 1, N 
        tmp2       = H(j, N)/(2.d0*Mass*(E(j) -Kinet))
        inner_a(j) = tmp1*tmp2/R(l)
    end do     
end subroutine inner_coeff









! ==================================================
! PROCESS
! ==================================================


subroutine PROC_inner_achive(l)
    integer(I4), intent(in) :: l 
    complex(QP) :: sum 
    integer(I4) :: i, j

    allocate(inner_a(1:N))
    call inner_coeff(l) 
    do i = 1, N 
        sum = 0.d0 
        do j = 1, N 
            sum = sum +inner_a(j)*H(j, i)
        end do 
        inner_u(l, i) = sum 
    end do 
    deallocate(inner_a)
end subroutine PROC_inner_achive


subroutine PROC_inner_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    use math_const,  only: pi => math_pi
    use hamiltonian, only: coord_r, coord_theta
    use nr, only: plgndr_s
    integer  (I1), parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real     (SP) :: tmp 
    complex  (QP) :: sum 
    integer  (I4) :: i, j, k 

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
                sum = sum +inner_u(k, i)/coord_r(i)*plgndr_s(k, 0, tmp)
            end do 
            write(file_psi2, form_psi) coord_r(i), coord_theta(j), dble(abs(sum)**2.d0)
        end do 
        write(file_psi2, form_psi) 
    end do 
    close(file_psi2)
end subroutine PROC_inner_plot
end module inner
