module inner
    use global 
    implicit none
    complex(8), save, allocatable, protected :: a(:)
contains 
    

subroutine inner_coeff(l)
    integer,    intent(in) :: l 
    real(8),    parameter  :: pi = 2.0d0*acos(0.0d0)
    complex(8), parameter  :: i  = (0.d0, 1.d0)
    real(8)    :: ka, scka, tmp2
    complex(8) :: tmp1 
    integer    :: j

    ka   = (2.d0*Mass*Kinet)**0.5d0*ra
    scka = ka -dble(l)*pi/2.d0
    tmp1 = (exp(-i*scka) -S(l)*exp(i*scka))/ka

    do j = 1, N 
        tmp2 = H(j, N)/(2.d0*Mass*(E(j)**2.d0 -Kinet**2.d0))
        a(j) = tmp1*tmp2/R(l)
    end do     
end subroutine inner_coeff










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_inner_achive(l)
    integer, intent(in) :: l 
    complex(16) :: tmp 
    integer     :: i, j

    allocate(a(1:N))
    call inner_coeff(l) 
    do i = 1, N 
        tmp = 0.d0 
        do j = 1, N 
            tmp = tmp +a(j)*H(j, i)
        end do 
        inner_u(l, i) = tmp 
    end do 
    deallocate(a)
end subroutine PROC_inner_achive


subroutine PROC_inner_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer,       parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real(8), parameter :: pi = 2.0d0*acos(0.0d0) 
    real(8), parameter :: radian_to_degree = 180.d0/pi 
    complex(16) :: tmp 
    integer     :: i, j, k 

    open(file_psi1, file = "inout/inner_psi.d")
    open(file_psi2, file = "inout/inner_psi-4pir.d")
    do i = 1, N, N/pr 
        do j = 0, ptheta
            tmp = 0.d0 
            do k = 0, L 
                tmp = tmp +inner_u(k, i)/coord_r(i)*plgndr_s(i, 0, cos(coord_theta(j)))
            end do 
            write(file_psi1, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(tmp)**2.d0)
            tmp = tmp*4.d0*pi*coord_r(i)
            write(file_psi2, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(tmp)**2.d0)
        end do 
        write(file_psi1, form_psi) 
        write(file_psi2, form_psi) 
    end do 
    close(file_psi1)
    close(file_psi2)
end subroutine PROC_inner_plot
end module inner
