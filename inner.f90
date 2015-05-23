module inner
    use global 
    implicit none
    complex(8), save, allocatable :: a(:)
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
end module inner










module PROC_inner
    use inner 
    implicit none
contains


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


subroutine PROC_inner_plot(l) ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer,       intent(in) :: l 
    integer,       parameter  :: file_psi = 101
    character(30), parameter  :: form_psi = '(30ES25.10)'
    integer     :: i, j

    open(file_psi, file = "inout/inner_psi.d")
    write(file_psi, form_psi) 0.d0, 0.d0 
    do i = 1, N 
        write(file_psi, form_psi) coord_r(i), dble(abs(inner_u(l, i))**2.d0)
    end do 
    close(file_psi)
end subroutine PROC_inner_plot
end module PROC_inner
