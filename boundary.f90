module boundary
    use global 
    implicit none
contains 


subroutine mat_R(l)
    integer, intent(in) :: l 
    real(16) :: tmp 
    integer  :: i, j 

    tmp = 0.d0 
    do i = 1, N 
        tmp = tmp +H(i, N)**2.d0/(2.d0*Mass*(E(i)**2.d0 -Kinet**2.d0))
    end do 
    R(l) = tmp/ra 
    write(*, *) "R: ", l, R(l)
end subroutine mat_R


subroutine mat_K(l)
    integer, intent(in) :: l 
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: ka, scka 

    ka   = (2.d0*Mass*Kinet)**0.5d0*ra
    scka = ka -dble(l)*pi/2.d0
    K(l) = (-sin(scka) +R(l)*ka*cos(scka))/(cos(scka) +R(l)*ka*sin(scka))
    write(*, *) "K: ", l, K(l)
end subroutine mat_K


subroutine mat_S(l)
    integer,    intent(in) :: l 
    complex(8), parameter  :: i = (0.d0, 1.d0)

    S(l) = (1.d0 +i*K(l))/(1.d0 -i*K(l))
    write(*, *) "S: ", l, S(l)
end subroutine mat_S
end module boundary










module PROC_boundary
    use boundary
    implicit none
contains
    

subroutine PROC_boundary_mat(l) ! It must be called after PROC_input, PROC_H 
    integer, intent(in) :: l 

    call mat_R(l)
    call mat_K(l)
    call mat_S(l)    
end subroutine PROC_boundary_mat
end module PROC_boundary
