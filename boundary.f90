module boundary
    use global 
    implicit none
contains 


subroutine mat_R(l)
    character(30), parameter :: form_out = '(1A15, 1I15, 1ES15.3)'
    integer, intent(in) :: l 
    real(16) :: tmp 
    integer  :: i, j 

    tmp = 0.d0 
    do i = 1, N 
        tmp = tmp +H(i, N)**2.d0/(2.d0*Mass*(E(i) -Kinet))
    end do 
    R(l) = tmp/ra 
    write(*, form_out) "R: ", l, R(l)
end subroutine mat_R


subroutine mat_K(l)
    character(30), parameter :: form_out = '(1A15, 1I15, 1ES15.3)'
    integer, intent(in) :: l 
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: ka, scka 

    ka   = (2.d0*Mass*Kinet)**0.5d0*ra
    scka = ka -dble(l)*pi/2.d0
    K(l) = (-sin(scka) +R(l)*ka*cos(scka))/(cos(scka) +R(l)*ka*sin(scka))
    write(*, form_out) "K: ", l, K(l)
end subroutine mat_K


subroutine mat_S(l)
    character(60), parameter :: form_out = '(1A15, 1I15, 1ES15.3, 1ES15.3, "i")'
    integer,    intent(in) :: l 
    complex(8), parameter  :: i = (0.d0, 1.d0)

    S(l) = (1.d0 +i*K(l))/(1.d0 -i*K(l))
    write(*, form_out) "S: ", l, S(l)
end subroutine mat_S










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_boundary_mat(l) ! It must be called after PROC_input, PROC_H 
    integer, intent(in) :: l 

    call mat_R(l)
!     write(*, *) "Here. (1)"
    call mat_K(l)
!     write(*, *) "Here. (2)"
    call mat_S(l)  
!     write(*, *) "Here. (3)"
end subroutine PROC_boundary_mat
end module boundary
