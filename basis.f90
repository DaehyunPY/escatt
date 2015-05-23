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
    nabla_x = nabla_x/dr**2.d0
    if(mod(int_i +int_j, 2) == 0) then 
        nabla_x = -nabla_x
    end if 
!        *(-1.d0)**(i -j +1.d0)
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
    if(info /= 0) stop "subroutine diag: Error. (2)"  
    deallocate(work)
end subroutine diag


subroutine correct(num, ty)
    integer, intent(in)  :: num 
    integer, intent(out) :: ty   
    real(8)  :: tmp, sign, c1, c2     
    real(16) :: sum 
    integer  :: i, j 

    sum = 0.d0 
    do i = 1, N 
        tmp = (tmp_H(2*i -1, 2*num -1) +tmp_H(2*i, 2*num -1))*2.d0**(-0.5d0) &
                /dr**0.5d0 
        sum = sum +tmp*tmp*dr 
    end do 
    write(*, *) "before corrected: ", num, dble(sum)

    sum = 0.d0 
    do i = 1, N 
        tmp = (tmp_H(2*i -1, 2*num) +tmp_H(2*i, 2*num))*2.d0**(-0.5d0) &
                /dr**0.5d0 
        sum = sum +tmp*tmp*dr 
    end do 
    if(sum < 1.d-20) ty = 0 
    write(*, *) "correct: ", num, dble(sum)

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

    sum = 0.d0 
    do i = 1, N 
        tmp = (tmp_H(2*i -1, 2*num -1) +tmp_H(2*i, 2*num -1))*2.d0**(-0.5d0) &
                /dr**0.5d0 
        sum = sum +tmp*tmp*dr 
    end do 
    write(*, *) "after corrected: ", num, dble(sum)
end subroutine correct


subroutine stnad
    real(16) :: sum 
    integer  :: i, j 

    do j = 1, N
        sum = 0.d0 
        do i = 1, N 
            sum = sum +H(j, i)*H(j, i)*dr
        end do 
        write(*, *) j, dble(sum)
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
    integer :: i, j, u, v, num, ty 

!     write(*, *) "Here. (1-1)"
    allocate(tmp_H(2*N, 2*N), tmp_E(2*N))
    tmp_H(:, :) = 0.d0 
    do i = 1, 2*N 
        do j = 1, 2*N 
            if(i <= N) then 
                u = 2*i -1 
!                 u = 2*(N -i +1) -1 
            else 
                u = 2*(2*N -i +1)
!                 u = 2*(i -N)
            end if 
            if(j <= N) then 
                v = 2*j -1 
!                 v = 2*(N -j +1) -1 
            else 
                v = 2*(2*N -j +1)
!                 v = 2*(j -N)
            end if 
            tmp_H(u, v) = tmp_H(u, v) -1.d0/(2.d0*Mass)*nabla_x(i, j)
        enddo
    enddo
!     write(*, *) "Here. (1-2)"
    do i = 1, N 
        u = 2*i -1 
!         u = 2*(N -i +1) -1 
        v = 2*i 
!         v = 2*(N -i +1) -1 
        tmp = 1.d0/(2.d0*Mass)*dble(l)*(dble(l) +1.d0)/cood_r(i)**2.d0 &
                +Poten(cood_r(i))
        tmp_H(u, u) = tmp_H(u, u) +tmp 
        tmp_H(v, v) = tmp_H(v, v) +tmp 
    enddo

!     write(*, *) "Here. (2)"
    call diag

!     write(*, *) "Here. (3-1)"
    num     = 0 
    H(:, :) = 0.d0
    E(:)    = 0.d0 
    ty      = 1 
    do j = 1, 2*N, 2 
        num = num +1 
!         call correct(num, ty)
        if(ty == 1) call correct(num, ty)
        if(tmp_H(1, j) > 0.d0) then 
            sign = 1.d0 
        else 
            sign = -1.d0 
        end if 
!         write(*, *) "Here. (3-2)"
        do i = 1, N 
            H(num, i) = sign*(tmp_H(2*i -1, j) +tmp_H(2*i, j))*2.d0**(-0.5d0) &
                        /dr**0.5d0 
        end do 
        E(num) = tmp_E(j)
    end do 
    if(num /= N) stop "subroutine PROC_H: Error. (1)"
!     write(*, *) "Here. (3-3)"
    deallocate(tmp_H, tmp_E)

!     write(*, *) "Here. (4)"
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
        write(file_psi, form_psi) cood_r(i), (H(j, i), j = 1, 5), H(N, i)
        write(file_ene, form_ene) i, E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
end module PROC_basis
