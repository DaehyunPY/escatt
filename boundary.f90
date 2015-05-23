module boundary_process
    use basis_global
    implicit none
    real(8), parameter :: Kinet = 1.d0 
    real(8)    :: R, K
    complex(8) :: S
contains 
    

subroutine PROC_R
    real(16) :: tmp 
    real(8)  :: conv1, conv2 
    integer  :: i, j 

    tmp = 0.d0 
    do i = 1, N 
        tmp = tmp +H(i, N)**2.d0/(2.d0*mass*(E(i)**2.d0 -Kinet**2.d0))
    end do 
    R = tmp/r_a 
    write(*, *) "R: ", R 

!     tmp = H(1, N)**2.d0/(2.d0*mass*(E(1)**2.d0 -Kinet**2.d0))
!     do i = 2, 1000
!         conv1 = H(i -1, N)**2.d0/(2.d0*mass*(E(i -1)**2.d0 -Kinet**2.d0))
!         conv2 = H(i, N)**2.d0/(2.d0*mass*(E(i)**2.d0 -Kinet**2.d0))
!         tmp = tmp +conv2
!         write(*, *) i, conv2/conv1 
!     end do 
!     R = tmp/r_a
!     write(*, *) "R: ", R 
end subroutine PROC_R


subroutine PROC_K
    real(8) :: p 

    p = (2.d0*mass*Kinet)**0.5d0*r_a
    K = (-sin(p) +R*p*cos(p))/(cos(p) +R*p*sin(p))
    write(*, *) "K: ", K 
end subroutine PROC_K


subroutine PROC_S
    complex(8), parameter :: i = (0.d0, 1.d0)

    S = (1.d0 +i*K)/(1.d0 -i*K)
    write(*, *) "S: ", S
end subroutine PROC_S
end module boundary_process
