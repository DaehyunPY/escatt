module cs_process
    use global 
    implicit none
    real(8) :: sigma
    complex(8) :: f
contains 
    

subroutine PROC_sigma
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)
    real(8) :: tmp 

    tmp = 1.d0/(1.d0 +1.d0/K**2.d0)
    sigma = 4.d0*pi/(2.d0*Mass*Kinet)*tmp 
    write(*, *) "sigma: ", sigma 
end subroutine PROC_sigma

subroutine PROC_dsigma
    complex(8), parameter :: i = (0.d0, 1.d0)
    real(8) :: tmp 

    tmp = (2.d0*Mass*Kinet)**0.50 
    f = 1.d0/(2.d0*i*tmp)*(S -1.d0)
    write(*, *) "f: ", f
end subroutine PROC_dsigma
end module cs_process
