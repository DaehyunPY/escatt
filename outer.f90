module outer
    use kind_type
    use global 
    use nr, only: plgndr_s
    implicit none
    complex(CDP), save, allocatable, private, protected :: outer_f(:)
contains 
    

subroutine mat_f 
    use math_const, only: i => math_i, pi => math_pi 
    real   (RDP) :: k 
    integer(I4B) :: j 

    k = (2.d0*Mass*Kinet)**0.50
    do j = 0, L 
        outer_f(j) = (2.d0*dble(j) +1.d0)/(2.d0*i*k)*(S(j) -1.d0)
    end do 
end subroutine mat_f 


function outer_u(l, r)
    integer(I4B), intent(in) :: l 
    real   (RDP), intent(in) :: r 
    real   (RSP) :: kr
    complex(CDP) :: outer_u

    kr      = (2.d0*Mass*Kinet)**0.5d0*r
!     outer_u = A(l)*(bessel_j(l, kr) -K(l)*bessel_y(l, kr))*r 
    outer_u = A(l)*(bessel_j(l, kr) -K(l)*bessel_y(l, kr))
end function outer_u










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_CS_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use hamiltonian, only: coord_r, coord_theta
    integer (I1B), parameter :: file_dcs = 101
    character(30), parameter :: form_dcs = '(30ES25.10)'
    character(30), parameter :: form_out = '(1A15, 5X, 1ES25.10)'
    real    (RDP), parameter :: radian_to_degree = 1.d0/degree 
    real    (RSP) :: tmp 
    complex (CQP) :: sum1, sum2  
    real    (RDP) :: k 
    integer (I4B) :: i, j 

    k = (2.d0*Mass*Kinet)**0.50
    allocate(outer_f(0:L))
    
    call mat_f 

    sum1 = 0.d0 
    do i = 0, L 
        sum1 = sum1 +outer_f(i)
!         write(*, *) sum1 
    end do 
    sum1 = 4.d0*pi/k*sum1
    write(*, form_out) "total sigma: ", dble(aimag(sum1))

    open(file_dcs, file = "inout/diff_cs.d")
    do j = 0, ptheta 
        sum2 = 0.d0 
        do i = 0, L 
            tmp = cos(coord_theta(j))
            sum2 = sum2 +outer_f(i)*plgndr_s(i, 0, tmp)
        end do 
        write(file_dcs, form_dcs) coord_theta(j)*radian_to_degree, abs(sum2)**2.d0
!         write(file_dcs, form_dcs) coord_theta(j)*radian_to_degree, abs(sum2/sum1)**2.d0
    end do 
    close(file_dcs)
    deallocate(outer_f)
end subroutine PROC_CS_plot


subroutine PROC_outer_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use hamiltonian, only: coord_r, coord_theta
    integer (I1B), parameter :: file_psi1 = 101, file_psi2 = 102, file_psi3 = 103
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real    (RDP), parameter :: radian_to_degree = 1.d0/degree
    real    (RSP) :: tmp 
    complex (CQP) :: sum 
    integer (I4B) :: i, j, k 

    open(file_psi1, file = "inout/outer_u0.d")
    sum = 0.d0 
    do i = N +1, 2*N
        sum = outer_u(0, coord_r(i))
        write(file_psi1, form_psi) coord_r(i), dble(abs(sum)**2.d0)
    end do 
    close(file_psi1)

    open(file_psi2, file = "inout/outer_psi.d")
    open(file_psi3, file = "inout/outer_psi-4pir.d")
    do i = N +1, 2*N, N/pr 
        do j = 0, ptheta
            sum = 0.d0 
            do k = 0, L 
                tmp = cos(coord_theta(j))
                sum = sum +outer_u(k, coord_r(i))/coord_r(i)*plgndr_s(i, 0, tmp)
            end do 
!             write(file_psi2, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(sum))**2.d0 
            write(file_psi2, form_psi) coord_r(i), coord_theta(j), dble(abs(sum))**2.d0 
            if(j == 0) write(file_psi2, form_psi) coord_r(i), dble(abs(sum)**2.d0)
            sum = sum*4.d0*pi*coord_r(i)
!             write(file_psi3, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(sum))**2.d0
            write(file_psi3, form_psi) coord_r(i), coord_theta(j), dble(abs(sum))**2.d0
        end do 
        write(file_psi2, form_psi) 
        write(file_psi3, form_psi) 
    end do 
    close(file_psi2)
    close(file_psi3)
end subroutine PROC_outer_plot
end module outer
