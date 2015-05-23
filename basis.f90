module basis_process
    use global 
    implicit none
contains


function nabla_x(ty_i, ty_j)
    integer, intent(in) :: ty_i, ty_j
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: nabla_x, i, j

    i = dble(ty_i)
    j = dble(ty_j)
    if( ty_i /= ty_j ) then 
        nabla_x = &
            2.d0/(i -j)**2.d0
    else
        nabla_x = &
            pi**2.d0 /3.d0
    end if
    nabla_x = nabla_x &
        *(-1.d0)/dr**2.d0 &
        *(-1.d0)**(i -j)
end function nabla_x


subroutine PROC_H
    real(8), allocatable :: tmp_H(:, :), tmp_E(:)
    real(8) :: sign
    integer :: i, j, k, num, ch 

    allocate(tmp_H(2*N, 2*N), tmp_E(2*N))
    tmp_H(:, :) = 0.d0  
    do j = 1, 2*N
        do i = 1, 2*N
            tmp_H(i, j) = -1.d0/(2.d0*Mass)*nabla_x(i, j)
        enddo
    enddo
    do i = 1, N
        tmp_H(i, i) = tmp_H(i, i) +Poten(cood_r(i))
        tmp_H(2*N +1 -i, 2*N +1 -i) = tmp_H(2*N +1 -i, 2*N +1 -i) +Poten(cood_r(i))
    enddo

    call PROC_diag(tmp_H, tmp_E)

    num = 0 
    do k = 1, N
        if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) > 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) < 0.d0) then 
                j = 2*k -1 
                num = num +1
                ch = -1 
        else if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) < 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) > 0.d0) then 
                j = 2*k 
                num = num +1
                ch = 0
        else if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) > 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) > 0.d0) then 
                j = 2*k +ch 
                num = num +1
                write(*, *) "PROCESS H: Warning. (1)", N, j 
        else 
            stop "PROCESS H: Somthing is wrong. (2)"
        end if 
        E(num) = tmp_E(j)
        sign = 1.d0 
        if(tmp_H(1, j) < 0.d0) then 
            sign = -1.d0 
        end if 
        do i = 1, N 
            H(num, i) = sign*(tmp_H(i, j) +tmp_H(2*N -i, j))*2.d0**(-0.5d0)
        end do 
    end do 
    if(num /= N) then 
        write(*, *) "PROCESS H: Something is wrong. (3)", N, num 
        stop 
    end if 
    deallocate(tmp_H, tmp_E)

    write(*,'(X, A, 10F10.3)') "Energy: ", (E(i), i = 1, 10)
end subroutine PROC_H


subroutine PROC_diag(H, E)
    real(8), intent(inout)  :: H(1:, 1:), E(1:)
    character(1), parameter :: jobz = 'Vectors', uplo = 'Upper' 
    integer :: n, lda, lwork, info, tmp 
    real(8), allocatable :: work(:)

        n = size(H(:, 1))
      tmp = size(E(:))
    if(n /= tmp) stop "PROCESS diag: Something is wrong. (1)"

      lda = size(H(1, :))
    lwork = 3*n -1
    allocate(work(1:3*n -1))
     info = 0

    !   lwork = -1
    !   call dsyev(jobz, uplo, n, H, lda, E, work, lwork, info)
    !   lwork = int(work(1))
    call DSYEV(jobz, uplo, n, H, lda, E, work, lwork, info)
    if(info /= 0) stop "PROCESS diag: Something is wrong. (2)"  
    deallocate(work)
end subroutine PROC_diag


subroutine PROC_psi
    real(16) :: sum 
    integer  :: i, j 

    do j = 1, N
        sum = 0.d0 
        do i = 1, N 
            sum = sum +H(j, i)*H(j, i)*dr
        end do 
        H(j, :) = H(j, :)/sum**0.5d0 
    end do  
end subroutine PROC_psi


subroutine PROC_plot
    integer,       parameter :: file_psi = 101
    character(30), parameter :: form_psi = '(30ES25.10)'
    integer :: i, j 

    open(file_psi, file = "inout/basis_psi.d")
    write(file_psi, form_psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    do i = 1, N 
        write(file_psi, form_psi) cood_r(i), (H(j, i), j = 1, 5)
    end do 
    close(file_psi)
end subroutine PROC_plot
end module basis_process
