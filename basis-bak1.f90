module basis
    use global 
    implicit none
    real(8), save, allocatable :: tmp_H(:, :), tmp_E(:)
contains


function nabla_x(int_i, int_j)
    integer, intent(in) :: int_i, int_j
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: nabla_x, i, j

    i = dble(int_i)
    j = dble(int_j)
    if(int_i /= int_j) then 
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


subroutine diag
    character(1), parameter :: jobz = 'Vectors', uplo = 'Upper' 
    integer :: n, lda, lwork, info, tmp 
    real(8), allocatable :: work(:)

!     n     = 2*N
!     lda   = n
    n     = size(tmp_H(:, 1))
    lda   = size(tmp_H(1, :))
    lwork = 3*n -1
    allocate(work(1:3*n -1))
    info  = 0

!     lwork = -1
!     call DSYEV(jobz, uplo, n, tmp_H, lda, tmp_E, work, lwork, info)
!     lwork = int(work(1))
    call DSYEV(jobz, uplo, n, tmp_H, lda, tmp_E, work, lwork, info)
    if(info /= 0) stop "subroutine diag: Error (2)"  
    deallocate(work)
end subroutine diag


subroutine stnad
    real(16) :: sum 
    integer  :: i, j 

    do j = 1, N
        sum = 0.d0 
        do i = 1, N 
            sum = sum +H(j, i)*H(j, i)*dr
        end do 
        H(j, :) = H(j, :)/sum**0.5d0 
    end do  
end subroutine stnad
end module basis










module PROC_basis
    use basis 
    implicit none
contains


subroutine PROC_H(l) ! It must be called after PROC_input
    integer, intent(in) :: l 
    character(30), parameter :: form_out = '(1X, 1A, 10F10.3)'
    real(8) :: sign, tmp 
    integer :: i, j, k, num, ch 

    allocate(tmp_H(2*N, 2*N), tmp_E(2*N))
    tmp_H(:, :) = 0.d0  
    do j = 1, 2*N
        do i = 1, 2*N
            tmp_H(i, j) = tmp_H(i, j) -1.d0/(2.d0*Mass)*nabla_x(i, j)
        enddo
    enddo
    do i = 1, N 
        tmp = 1.d0/(2.d0*Mass)*dble(l)*(dble(l) +1.d0)/cood_r(i)**2.d0 &
                +Poten(cood_r(i))
        tmp_H(i, i) = tmp_H(i, i) +tmp 
        tmp_H(2*N +1 -i, 2*N +1 -i) = tmp_H(2*N +1 -i, 2*N +1 -i) +tmp 
    enddo

    call diag

    num     = 0
    ch      = -1 
    H(:, :) = 0.d0
    E(:)    = 0.d0
    do k = 1, N
        if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) > 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) < 0.d0) then 
                j   = 2*k -1
                num = num +1
                ch  = -1
        else if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) < 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) > 0.d0) then 
                j   = 2*k
                num = num +1
                ch  = 0
        else if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) > 0.d0 & 
            .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) > 0.d0) then 
                j   = 2*k +ch
                num = num +1
                write(*, *) "subroutine PROC_H: Warning (1)", N, j 
        else 
                j   = 2*k +ch
                num = num +1
                write(*, *) "subroutine PROC_H: Warning (2)", N, j 
        end if 
        E(num) = tmp_E(j)
        sign   = 1.d0 
        if(tmp_H(1, j) < 0.d0) then 
            sign = -1.d0 
        end if 
        do i = 1, N 
            H(num, i) = sign*(tmp_H(i, j) +tmp_H(2*N -i, j))*2.d0**(-0.5d0)
        end do 
    end do 
    if(num /= N) then 
        write(*, *) "subroutine PROC_H: Error (3)", N, num 
        stop 
    end if 
    deallocate(tmp_H, tmp_E)

    call stnad
    write(*, form_out) "Energy: ", (E(i), i = 1, 10)
end subroutine PROC_H


subroutine PROC_basis_plot ! It must be called after PROC_input, PROC_H 
    integer,       parameter :: file_psi = 101,           file_ene = 102
    character(30), parameter :: form_psi = '(30ES25.10)', form_ene = '(1I5, 1ES25.10)'
    integer :: i, j 

    open(file_psi, file = "inout/basis_psi.d")
    open(file_ene, file = "inout/basis_energy.d")
    write(file_psi, form_psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    do i = 1, N 
        write(file_psi, form_psi) cood_r(i), (H(j, i), j = 1, 5)
        write(file_ene, form_ene) i, E(i)
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
end module PROC_basis
