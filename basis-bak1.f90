module basis
    use kind_type 
    use global 
    implicit none
    real(dp), save, allocatable, private, protected :: tmp_H(:, :), tmp_E(:)
contains


! subroutine diag
!     character(1), parameter :: jobz = 'Vectors', uplo = 'Upper' ! 'Upper' or 'Lower'
!     integer (i8) :: n, lda, lwork
!     integer (i1) :: info 
!     real    (dp), allocatable :: work(:)

!     n     = size(tmp_H(:, 1))
!     lda   = size(tmp_H(1, :))
!     lwork = 3*n -1
!     allocate(work(1:lwork))
!     info  = 0

!     lwork = -1
!     call DSYEV(jobz, uplo, n, tmp_H, lda, tmp_E, work, lwork, info)
!     lwork = int(work(1))
!     deallocate(work)
!     allocate(work(1:lwork))

!     call DSYEV(jobz, uplo, n, tmp_H, lda, tmp_E, work, lwork, info)
!     if(info /= 0) stop "subroutine diag: Error. (2)"  
!     deallocate(work)
! end subroutine diag
subroutine diag
    character(1), parameter :: jobvl = 'No Vectors', jobvr = 'Vectors' ! 'Vectors' or 'No Vecters'
    integer (i8) :: n, lda, ldvl, ldvr, lwork 
    real    (dp), allocatable :: wr(:), wi(:), vl(:, :), vr(:, :), work(:)
    integer (i1) :: info 

    n     = size(tmp_H(:, 1))
    lda   = size(tmp_H(1, :))
    ldvl  = 1
    ldvr  = n 
    work  = 3*n -1
    info  = 0 
    allocate(wr(1:n), wi(1:n), vl(1, 1), vr(1:ldvr, 1:n), work(1:lwork))

    lwork = -1
    call DGEEV(jobvl, jobvr, n, tmp_H, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(1:lwork))

    write(*, *) "Here. (3)"
    call DGEEV(jobvl, jobvr, n, tmp_H, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    if(info /= 0) stop "subroutine diag: Error. (2)"  
    write(*, *) "Here. (4)"
    tmp_H(:, :) = vr(:, :)
    tmp_E(:)    = wr(:)
    write(*, *) "Here. (5)"
    deallocate(wr, wi, vl, vr, work)
    write(*, *) "Here. (6)"
end subroutine diag


! subroutine correct(num, ty)
!     use hamiltonian, only: dr 
!     character(30), parameter :: form_out = '(1A15, 1I15, 1ES15.3)'
!     integer  (i4), intent(in)  :: num 
!     integer  (i1), intent(out) :: ty   
!     real     (dp) :: tmp, sign, c1, c2
!     real     (qp) :: sum 
!     integer  (i4) :: i, j

!     sum = 0.d0 
!     do i = 1, N 
!         tmp = (tmp_H(2*i -1, 2*num) +tmp_H(2*i, 2*num))*2.d0**(-0.5d0) &
!                 /dr**0.5d0 
!         sum = sum +tmp*tmp*dr 
!     end do 
!     if(sum < 1.d-20) ty = 0 
!     write(file_log, form_out) "correct: ", num, dble(sum) 

!     if(tmp_H(1, 2*num -1)*tmp_H(1, 2*num) > 0.d0) then 
!         sign = 1.d0 
!     else 
!         sign = -1.d0 
!     end if 
!     c1 = (1.d0 -sum)**0.5d0 
!     c2 = sum**0.5d0 

!     tmp_H(:, 2*num -1) = c1*tmp_H(:, 2*num -1) +sign*c2*tmp_H(:, 2*num)
! !     tmp_E(2*num -1) = c1**2.d0*tmp_E(2*num -1) +c2**2.d0*tmp_E(2*num)
!     tmp_E(2*num -1) = (1.d0 -sum)*tmp_E(2*num -1) +sum*tmp_E(2*num)
! end subroutine correct
subroutine check_mat
    use hamiltonian, only: dr 
    integer  (i1), parameter :: file_check = 101
    character(30), parameter :: & 
        form_num   = '(20X, 5000(I15, 5X))', &
        form_check = '(20X, 5000(ES15.5, 5X))', &
        form_ch    = '(1I10)', &
        form_ch1   = '(1I15, 5X, ', &
        form_ch2   = '20X, ', &
        form_ch3   = '(ES15.5, 5X), 20X, ', &
        form_ch4   = '(ES15.5, 5X))'
    character(10) :: ch1, ch2 
    integer (i4) :: i, j 

    open(file_check, file = "output/check_mat.d")
!     write(file_check, form_num) 0, (j, j = 1, n)
!     do i = 1, n 
!         write(file_check, form_check) i, (tmp_H(i, j), j = 1, n)
!     end do 
!     write(form_check, *)
!     write(form_check, *)
!     write(form_check, *)
!     write(form_check, *)
!     write(form_check, *)
    write(file_check, form_num) (j, j = 1, n)
    do i = 1, n 
        if(i == 1) then 
            write(ch1, form_ch) 0 
            write(ch2, form_ch) n -1
            write(file_check, form_ch1//form_ch2//ch2//form_ch4) i, (tmp_H(i, j), j = 2, n)
            ch1 = ""
            ch2 = ""
        else if(i == n) then 
            write(ch1, form_ch) n -1 
            write(ch2, form_ch) 0 
            write(file_check, form_ch1//ch1//form_ch4) i, (tmp_H(i, j), j = 1, n -1) 
            ch1 = ""
            ch2 = ""
        else 
            write(ch1, form_ch) i -1 
            write(ch2, form_ch) n -i 
            write(file_check, form_ch1//ch1//form_ch3//ch2//form_ch4) i, (tmp_H(i, j), j = 1, i -1), (tmp_H(i, j), j = i +1, n)
            ch1 = ""
            ch2 = ""
        end if 
    end do 
    write(file_check, *)
    write(file_check, form_num) (j, j = 1, n)
    write(file_check, form_check) (tmp_H(j, j), j = 1, n)
    write(file_check, *)
    close(file_check)
end subroutine check_mat


subroutine stnad
    use hamiltonian, only: dr 
    real   (qp) :: sum 
    integer(i4) :: i, j 

    do j = 1, N
        sum = 0.d0 
        do i = 1, N 
            sum = sum +H(j, i)*H(j, i)*dr
        end do 
        write(file_log, *) j, dble(sum)
        H(j, :) = H(j, :)/sum**0.5d0 
    end do  
end subroutine stnad










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_H(l) 
    use hamiltonian, only: dr, coord_r, nabla_x, Poten
    character(30), parameter  :: form_out = '(1A15, 10F10.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: sign, tmp 
    integer  (i4) :: i, j, k

    allocate(tmp_H(1:N, 1:N), tmp_E(1:N))
    tmp_H(:, :) = 0.d0 
    do j = 1, 2*N -1 
        do i = 1, N             
            k   = j 
            tmp = -1.d0/(2.d0*Mass)*nabla_x(i, j) 
            if(j > N) k = 2*N -j 
            tmp_H(i, k) = tmp_H(i, k) +tmp 
        enddo
    enddo
    do i = 1, N 
        tmp         = 1.d0/(2.d0*Mass)*dble(l)*(dble(l) +1.d0)/coord_r(i)**2.d0 & ! azimuthal quantum term
                        +Poten(coord_r(i))                                        ! potential term
        tmp_H(i, i) = tmp_H(i, i) +tmp 
    enddo

    call check_mat
    call diag

    H(:, :) = 0_dp
    E(:)    = 0_dp 
    sign    = 1_dp 
    do j = 1, N 
        if(tmp_H(1, j) > 0.d0) then 
            sign = 1_dp 
        else 
            sign = -1_dp
        end if 
        H(j, :) = sign*tmp_H(:, j)/dr**0.5_dp
        E(j) = tmp_E(j)
    end do 
    deallocate(tmp_H, tmp_E)

!     call stnad ! for test 
    write(file_log, form_out) "Energy: ", (E(i), i = 1, 5)
end subroutine PROC_H


subroutine PROC_basis_plot(num)
    use hamiltonian, only: coord_r
    integer  (i1), parameter  :: file_psi = 101,           file_ene = 102
    character(30), parameter  :: form_psi = '(30ES25.10)', form_ene = '(1I5, 1ES25.10)'
    integer  (i4), intent(in) :: num 
    integer  (i4) :: i, j 
    character :: ch1*1, ch2*2, ch3*3, ch4*4

    if(L >= 0 .and. L < 10) then 
        write(ch1, '(I1.1)') num 
        open(file_psi, file = "output/basis_u_"//ch1//".d")
        open(file_ene, file = "output/basis_energy_"//ch1//".d")
    else if(L >= 10 .and. L < 100) then 
        write(ch2, '(I2.2)') num 
        open(file_psi, file = "output/basis_u_"//ch2//".d")
        open(file_ene, file = "output/basis_energy_"//ch2//".d")
    else if(L >= 100 .and. L < 1000) then 
        write(ch3, '(I3.3)') num 
        open(file_psi, file = "output/basis_u_"//ch3//".d")
        open(file_ene, file = "output/basis_energy_"//ch3//".d")
    else if(L >= 1000 .and. L < 10000) then 
        write(ch4, '(I4.4)') num 
        open(file_psi, file = "output/basis_u_"//ch4//".d")
        open(file_ene, file = "output/basis_energy_"//ch4//".d")
    end if 


    write(file_psi, form_psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 
    do i = 1, N 
!         write(file_psi, form_psi) coord_r(i), (H(j, i), j = 1, 5), H(n, i)
        write(file_psi, form_psi) coord_r(i), (H(j, i), j = n, n -4, -1), H(1, i)
        write(file_ene, form_ene) i, E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
end module basis
