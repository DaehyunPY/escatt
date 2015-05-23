module basis
    use kind_type 
    use global 
    implicit none
    real(dp), allocatable, private, protected :: mat_H(:, :)
contains


! ==================================================
! FUNCTION
! ==================================================
! kinetic term -------------------------------------
function mat_Kinet(i, j)
    use hamiltonian, only: dr_p_drho, Delta_rho 
    integer(i4), intent(in) :: i, j 
    real   (dp) :: mat_Kinet, tmp 
    tmp       = Delta_rho(i, j)/dr_p_drho**2_dp 
    mat_Kinet = -tmp/2_dp
end function mat_Kinet 
! potential term -----------------------------------
function mat_Poten(i)
    use hamiltonian, only: coord_weight, coord_r, Poten_r 
    integer(i4), intent(in) :: i
    real   (dp) :: mat_Poten
    mat_Poten = Poten_r(coord_r(i))
end function mat_Poten 
! angular term -------------------------------------
function mat_angular(l, i)
    use hamiltonian, only: coord_weight, coord_r, angular_r
    integer(i4), intent(in) :: l, i
    real   (dp) :: mat_angular
    mat_angular = angular_r(l, coord_r(i))
end function mat_angular


! ==================================================
! CALCULATION 
! ==================================================
! dsyev diagonalization ----------------------------
subroutine diag
    character(1), parameter   :: jobz = 'Vectors', uplo = 'Upper'
    real    (dp), allocatable :: work(:)
    integer (i8), save        :: lwork = -1, n, lda 
    integer (i1) :: info 
    
    if(lwork < 0) then 
        n    = size(mat_H(:, 1))
        lda  = size(mat_H(1, :))
        info = 0 
        if(allocated(work) == .false.) allocate(work(1))
        call DSYEV(jobz, uplo, n, mat_H, lda, E, work, lwork, info)
        lwork = int(work(1))
        if(allocated(work) == .true.) deallocate(work)
    end if 

    info = 0 
    if(allocated(work) == .false.) allocate(work(1:lwork))
    call DSYEV(jobz, uplo, n, mat_H, lda, E, work, lwork, info)
    if(info /= 0) stop "SUBROUTINE diag: Error. (2)"
    if(allocated(work) == .true.) deallocate(work)
end subroutine diag


! ==================================================
! TEST
! ==================================================
! check matrix -------------------------------------
! subroutine check_mat
!     integer  (i1), parameter :: file_check = 101
!     character(30), parameter :: & 
!         form_num   = '(20X, 5000(I15, 5X))', &
!         form_check = '(20X, 5000(ES15.5, 5X))', &
!         form_ch    = '(1I10)', &
!         form_ch1   = '(1I15, 5X, ', &
!         form_ch2   = '20X, ', &
!         form_ch3   = '(ES15.5, 5X), 20X, ', &
!         form_ch4   = '(ES15.5, 5X))'
!     character(10) :: ch1, ch2 
!     integer  (i4) :: n, i, j 

!     n = size(mat_H(:, 1))
!     open(file_check, file = "output/check_mat.d")
!     write(file_check, form_num) (j, j = 1, n)
!     do i = 1, n 
!         if(i == 1) then 
!             write(ch1, form_ch) 0 
!             write(ch2, form_ch) n -1
!             write(file_check, form_ch1//form_ch2//ch2//form_ch4) i, (mat_H(i, j), j = 2, n)
!             ch1 = ""
!             ch2 = ""
!         else if(i == n) then 
!             write(ch1, form_ch) n -1 
!             write(ch2, form_ch) 0 
!             write(file_check, form_ch1//ch1//form_ch4) i, (mat_H(i, j), j = 1, n -1) 
!             ch1 = ""
!             ch2 = ""
!         else 
!             write(ch1, form_ch) i -1 
!             write(ch2, form_ch) n -i 
!             write(file_check, form_ch1//ch1//form_ch3//ch2//form_ch4) i, (mat_H(i, j), j = 1, i -1), (mat_H(i, j), j = i +1, n)
!             ch1 = ""
!             ch2 = ""
!         end if 
!     end do 
!     write(file_check, *)
!     write(file_check, form_num) (j, j = 1, n)
!     write(file_check, form_check) (mat_H(j, j), j = 1, n)
!     write(file_check, *)
!     close(file_check)
! end subroutine check_mat










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    use hamiltonian, only: dr_p_drho, coord_weight
    character(30), parameter  :: form_out = '(1A15, 10F10.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: sign, tmp 
    integer  (i4) :: i, j

    if(allocated(mat_H) == .false.) allocate(mat_H(1:N, 1:N))
    mat_H(:, :) = 0_dp
    do i = 1, N 
        do j = 1, N
            mat_H(i, j) = mat_H(i, j) +mat_Kinet(i, j)
        enddo
        mat_H(i, i) = mat_H(i, i) +mat_Poten(i) +mat_angular(l, i)
    enddo
    call diag

    H(:, :) = 0_dp
    sign    = 1_dp 
    do j = 1, N 
        sign = 1_dp 
        if(mat_H(1, j) < 0_dp) sign = -1_dp 
        H(j, 0) = 0_dp 
        do i = 1, N 
            tmp = sign*mat_H(i, j)
            tmp = tmp
            H(j, i) = tmp/(coord_weight(i)*dr_p_drho)**0.5_dp
        end do 
    end do 
    if(allocated(mat_H) == .true.) deallocate(mat_H)
    write(file_log, form_out) "Energy: ", (E(i), i = 1, 5) 
end subroutine PROC_H
! end hamiltonian ----------------------------------
! basis plot ---------------------------------------
subroutine PROC_basis_plot(num)
    use hamiltonian, only: coord_rho, coord_r, coord_weight
    integer  (i1), parameter  :: file_psi = 101,                       file_ene = 102
    character(30), parameter  :: form_psi = '(1ES25.10, 1000ES25.10)', form_ene = '(1I5, 1ES25.10)'
    integer  (i4), intent(in) :: num 
    integer  (i4) :: i, j, n 
    character     :: ch1*1, ch2*2, ch3*3, ch4*4

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

    n = size(H(:, 1))
    write(file_psi, form_psi) (0_dp, j = 0, 20)
    do i = 1, n 
        write(file_psi, form_psi) coord_r(i), (H(j, i), j = 1, 10), (H(j, i), j = n -9, n)
!         write(file_psi, form_psi) (H(j, i), j = 1, 10), (H(j, i), j = n -9, n), coord_r(i)
        write(file_ene, form_ene) i, E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
! end basis plot -----------------------------------
end module basis
