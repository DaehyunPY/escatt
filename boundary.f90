module boundary_process
    use global 
    implicit none
contains 
    

subroutine PROC_R
    real(16) :: tmp 
    integer  :: i, j 

    tmp = 0.d0 
    do i = 1, N 
        tmp = tmp +H(i, N)**2.d0/(2.d0*Mass*(E(i)**2.d0 -Kinet**2.d0))
    end do 
    R = tmp/r_a 
    write(*, *) "R: ", R 
end subroutine PROC_R


subroutine PROC_K
    real(8) :: tmp 

    tmp = (2.d0*Mass*Kinet)**0.5d0*r_a
    K = (-sin(tmp) +R*tmp*cos(tmp))/(cos(tmp) +R*tmp*sin(tmp))
    write(*, *) "K: ", K 
end subroutine PROC_K


subroutine PROC_S
    complex(8), parameter :: i = (0.d0, 1.d0)

    S = (1.d0 +i*K)/(1.d0 -i*K)
    write(*, *) "S: ", S
end subroutine PROC_S
end module boundary_process
