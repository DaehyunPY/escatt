module process
    use input
    implicit none
contains


subroutine proc1
    integer i, j

    H(:, :) = 0.d0  
    do j = 1, CooN
    do i = 1, CooN
!         H(i, j) = nabla_r_close_open(dr, i, j)
        H(i, j) = nabla_r_close_close(dr, i, j)
    enddo
    enddo
    do i = 1, CooN
        H(i, i) = H(i, i) +Poten(CooR(i))
    enddo

    H(:, :) = -0.5d0*H(:, :)
end subroutine proc1


subroutine proc2_test
    integer,       parameter :: file_test = 100
    character(30), parameter :: form_test = '(1000ES)'
    integer i, j 
    open(file_test, file = "inout/test.d")
    do i = 1, CooN
        write(file_test, form_test) (H(i, j), j = 1, CooN)
    end do 
    close(file_test)
end subroutine proc2_test
subroutine proc2
    character(1) :: jobz = 'Vectors', uplo = 'Upper' 
    integer lda, n, lwork, info
    real(8), allocatable :: work(:)

    allocate(work(1:3*CooN -1))
      lda = CooN
        n = CooN
    lwork = 3*n -1
     info = 0

!   lwork = -1
!   call dsyev(jobz, uplo, n, H(1:n, 1:n), lda, E(1:n), work, lwork, info)
!   lwork = int(work(1))
    
    call dsyev(jobz, uplo, n, H(1:n, 1:n), lda, E(1:n), work, lwork, info)
    if(info /= 0) stop "SUBROUTINE diag_real: Something is wrong. (2)"  
end subroutine proc2


subroutine proc3 
    integer i, j 

    write(file_psi, form_psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 
    write(file_ana, form_ana) 0.d0, 0.d0, 0.d0, 0.d0  
    do i = 1, CooN 
        write(file_psi, form_psi) CooR(i), (abs(H(i, j))**2.d0, j = 1, 5)
        write(file_ana, form_ana) CooR(i), (f_ana(j, CooR(i))**2.d0, j = 1, 3)
        write(file_energy, form_energy) E(i)
    end do 
    
    write(*, '(10F10.3)') (E(i), i = 1, 10)
end subroutine proc3
end module process


program main
    use process
    implicit none

    open(file_psi, file = "inout/psi.d")
    open(file_energy, file = "inout/energy.d")
    open(file_ana, file = "inout/ana.d")

    call proc1
    call proc2_test 
    call proc2
    call proc3

    close(101)
    close(102)
    close(103)
end program main
