module outer
    use global 
    implicit none
    real(8), save :: sigma
    complex(8), save, allocatable :: f(:)
contains 
    

subroutine mat_f 
    complex(8), parameter :: i  = (0.d0, 1.d0)
    real(8),    parameter :: pi = 2.0d0*acos(0.0d0)
    real(8) :: k 
    integer :: j 

    k = (2.d0*Mass*Kinet)**0.50
    do j = 0, L 
        f(j) = 1.d0/(2.d0*i*k)*(2.d0*dble(j) +1.d0)*(S(j) -1.d0)
    end do 
end subroutine mat_f 


function outer_u(l, r)
    integer,    intent(in) :: l 
    real(8),    intent(in) :: r 
    real(8),    parameter  :: pi = 2.0d0*acos(0.0d0)
    complex(8), parameter  :: i = (0.d0, 1.d0)
    real(8)    :: k, sckr 
    complex(8) :: outer_u

    k       = (2.d0*Mass*Kinet)**0.5d0
    sckr    = k*r -dble(l)*pi/2.d0
    outer_u = (exp(-i*sckr) -S(l)*exp(i*sckr))/k 
end function outer_u
end module outer










module PROC_outer
    use outer 
    implicit none
contains


subroutine PROC_CS_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    real(8), parameter :: pi = 2.0d0*acos(0.0d0) 
    complex(16) :: tmp 
    real(8)     :: k 
    integer     :: i, j 

    k = (2.d0*Mass*Kinet)**0.50
    allocate(f(0:L))
    
    call mat_f 

    tmp = 0.d0 
    do i = 0, L 
        tmp = tmp +f(i)
        write(*, *) tmp 
    end do 
    tmp = 4.d0*pi/k*tmp

    write(*, *) "total sigma: ", dble(aimag(tmp))

!     call differential_CS(l)
    deallocate(f)
end subroutine PROC_CS_plot


subroutine PROC_outer_plot(l) ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer, intent(in) :: l 
    integer,       parameter :: file_psi = 101
    character(30), parameter :: form_psi = '(30ES25.10)'
    integer :: i 

    open(file_psi, file = "inout/outer_psi.d")
    do i = N +1, 2*N 
        write(file_psi, form_psi) coord_r(i), abs(outer_u(l, coord_r(i)))**2.d0 
    end do 
    close(file_psi)
end subroutine PROC_outer_plot
end module PROC_outer
