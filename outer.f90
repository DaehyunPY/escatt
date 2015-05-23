module outer
    use global 
    implicit none
    integer, parameter :: M = 2*N 
    real(8),    save   :: sigma
    complex(8), save   :: f
contains 
    

subroutine total_CS(l)
    integer, intent(in) :: l 
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: tmp 

    tmp   = 1.d0/(1.d0 +1.d0/K(l)**2.d0)
    sigma = 4.d0*pi/(2.d0*Mass*Kinet)*(2.d0*dble(l) +1.d0)*tmp
    write(*, *) "sigma: ", sigma 
end subroutine total_CS


subroutine differential_CS(l)
    integer,    intent(in) :: l 
    complex(8), parameter  :: i = (0.d0, 1.d0)
    real(8) :: tmp 

    tmp = (2.d0*Mass*Kinet)**0.50
    f   = 1.d0/(2.d0*i*tmp)*(2.d0*dble(l) +1.d0)*(S(l) -1.d0)
    write(*, *) "f: ", f
end subroutine differential_CS


function outer_f(l, r)
    integer,    intent(in) :: l 
    real(8),    intent(in) :: r 
    real(8),    parameter  :: pi = 2.0d0*acos(0.0d0)
    complex(8), parameter  :: i = (0.d0, 1.d0)
    real(8)    :: k, sckr 
    complex(8) :: outer_f

    k       = (2.d0*Mass*Kinet)**0.5d0
    sckr    = k*r -dble(l)*pi/2.d0
    outer_f = (exp(-i*sckr) -S(l)*exp(i*sckr))/k 
end function outer_f
end module outer










module PROC_outer
    use outer 
    implicit none
contains


subroutine PROC_CS(l) ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer, intent(in) :: l 

    call total_CS(l)
    call differential_CS(l)
end subroutine PROC_CS


subroutine PROC_outer_plot(l) ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer, intent(in) :: l 
    integer,       parameter :: file_psi = 101
    character(30), parameter :: form_psi = '(30ES25.10)'
    integer :: i 

    open(file_psi, file = "inout/outer_psi.d")
    do i = N +1, M 
        write(file_psi, form_psi) cood_r(i), abs(outer_f(l, cood_r(i)))**2.d0 
    end do 
    close(file_psi)
end subroutine PROC_outer_plot
end module PROC_outer
