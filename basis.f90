module basis
    use kind_type 
    use global 
    implicit none
    real(RDP), save, allocatable, private, protected :: tmp_H(:, :), tmp_E(:)
contains


subroutine diag
    character(1), parameter :: jobz = 'Vectors', uplo = 'Upper' 
    integer(I8B) :: n, lda, lwork
    integer(I1B) :: info 
    real   (RDP), allocatable :: work(:)

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
    if(info /= 0) stop "subroutine diag: Error. (2)"  
    deallocate(work)
end subroutine diag


subroutine correct(num, ty)
    character(30), parameter :: form_out = '(1A15, 1I15, 1ES15.3)'
    integer(I4B), intent(in)  :: num 
    integer(I1B), intent(out) :: ty   
    real   (RDP) :: tmp, sign, c1, c2     
    real   (RQP) :: sum 
    integer(I4B) :: i, j

    sum = 0.d0 
    do i = 1, N 
        tmp = (tmp_H(2*i -1, 2*num) +tmp_H(2*i, 2*num))*2.d0**(-0.5d0) &
                /dr**0.5d0 
        sum = sum +tmp*tmp*dr 
    end do 
    if(sum < 1.d-20) ty = 0 
    write(*, form_out) "correct: ", num, dble(sum) 

    if(tmp_H(1, 2*num -1)*tmp_H(1, 2*num) > 0.d0) then 
        sign = 1.d0 
    else 
        sign = -1.d0 
    end if 
    c1 = (1.d0 -sum)**0.5d0 
    c2 = sum**0.5d0 

    tmp_H(:, 2*num -1) = c1*tmp_H(:, 2*num -1) +sign*c2*tmp_H(:, 2*num)
!     tmp_E(2*num -1) = c1**2.d0*tmp_E(2*num -1) +c2**2.d0*tmp_E(2*num)
    tmp_E(2*num -1) = (1.d0 -sum)*tmp_E(2*num -1) +sum*tmp_E(2*num)
end subroutine correct


subroutine stnad
    real   (RQP) :: sum 
    integer(I4B) :: i, j 

    do j = 1, N
        sum = 0.d0 
        do i = 1, N 
            sum = sum +H(j, i)*H(j, i)*dr
        end do 
        write(*, *) j, dble(sum)
        H(j, :) = H(j, :)/sum**0.5d0 
    end do  
end subroutine stnad










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_H(l) 
    use hamiltonian, only: coord_r, nabla_x, Poten
    character(30), parameter  :: form_out = '(1A15, 10F10.3)'
    integer (I4B), intent(in) :: l 
    real    (RDP) :: sign, tmp 
    integer (I4B) :: i, j, u, v, num
    integer (I1B) :: ty 

    allocate(tmp_H(1:2*N, 1:2*N), tmp_E(1:2*N))
    tmp_H(:, :) = 0.d0 
    do i = 1, 2*N 
        do j = 1, 2*N 
            if(i <= N) then 
                u = 2*i -1 
            else 
                u = 2*(2*N -i +1)
            end if 
            if(j <= N) then 
                v = 2*j -1 
            else 
                v = 2*(2*N -j +1)
            end if 
            tmp_H(u, v) = tmp_H(u, v) -1.d0/(2.d0*Mass)*nabla_x(i, j)
        enddo
    enddo
    do i = 1, N 
        u = 2*i -1 
        v = 2*i 
        tmp = 1.d0/(2.d0*Mass)*dble(l)*(dble(l) +1.d0)/coord_r(i)**2.d0 &
                +Poten(coord_r(i))
        tmp_H(u, u) = tmp_H(u, u) +tmp 
        tmp_H(v, v) = tmp_H(v, v) +tmp 
    enddo

    call diag

    num     = 0 
    H(:, :) = 0.d0
    E(:)    = 0.d0 
    ty      = 1 
    do j = 1, 2*N, 2 
        num = num +1 
        if(ty == 1) call correct(num, ty)
        if(tmp_H(1, j) > 0.d0) then 
            sign = 1.d0 
        else 
            sign = -1.d0 
        end if 
        do i = 1, N 
            H(num, i) = sign*(tmp_H(2*i -1, j) +tmp_H(2*i, j))*2.d0**(-0.5d0) &
                        /dr**0.5d0 
        end do 
        E(num) = tmp_E(j)
    end do 
    if(num /= N) stop "subroutine PROC_H: Error. (1)"
    deallocate(tmp_H, tmp_E)

!     call stnad
    write(*, form_out) "Energy: ", (E(i), i = 1, 10)
end subroutine PROC_H


subroutine PROC_basis_plot 
    use hamiltonian, only: coord_r
    integer (I1B), parameter :: file_psi = 101,           file_ene = 102
    character(30), parameter :: form_psi = '(30ES25.10)', form_ene = '(1I5, 1ES25.10)'
    integer (I4B) :: i, j 

    open(file_psi, file = "inout/basis_psi.d")
    open(file_ene, file = "inout/basis_energy.d")
    write(file_psi, form_psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    do i = 1, N 
        write(file_psi, form_psi) coord_r(i), (H(j, i), j = 1, 5), H(N, i)
        write(file_ene, form_ene) i, E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
end module basis
