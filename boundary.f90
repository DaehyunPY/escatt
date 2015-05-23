module boundary
    use global 
    implicit none
contains 


subroutine mat_R
    real(16) :: tmp 
    integer  :: i, j 

    tmp = 0.d0 
    do i = 1, N 
        tmp = tmp +H(i, N)**2.d0/(2.d0*Mass*(E(i)**2.d0 -Kinet**2.d0))
    end do 
    R = tmp/ra 
    write(*, *) "R: ", R 
end subroutine mat_R


subroutine mat_K
    real(8) :: ka 

    ka = (2.d0*Mass*Kinet)**0.5d0*ra
    K = (-sin(ka) +R*ka*cos(ka))/(cos(ka) +R*ka*sin(ka))
    write(*, *) "K: ", K 
end subroutine mat_K


subroutine mat_S
    complex(8), parameter :: i = (0.d0, 1.d0)

    S = (1.d0 +i*K)/(1.d0 -i*K)
    write(*, *) "S: ", S
end subroutine mat_S
end module boundary










module PROC_boundary
    use boundary
    implicit none
contains
    

subroutine PROC_boundary_mat  ! It must be called after PROC_H 
    call mat_R
    call mat_K
    call mat_S    
end subroutine PROC_boundary_mat
end module PROC_boundary
