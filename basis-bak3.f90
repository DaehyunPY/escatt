module basis
    use kind_type 
    use global 
    implicit none
    real(dp), save, allocatable, private, protected :: tmp_H(:, :), tmp_E(:)
contains


function mat_Kinet(i, j)
    use hamiltonian, only: nabla_rho 
    integer(i4), intent(in) :: i, j 
    real   (dp) :: mat_Kinet, tmp 

    tmp = 0_dp 
    tmp = tmp +nabla_rho(i, j) 
    tmp = tmp/ra**2_dp 
    mat_Kinet = -1_dp/(2_dp*Mass)*tmp 
end function mat_Kinet 


function mat_Poten(i)
    use hamiltonian, only: coord_r, Poten_r 
    integer(i4), intent(in) :: i
    real   (dp) :: mat_Poten

    mat_Poten = Poten_r(coord_r(i))
end function mat_Poten 


function mat_angular(l, i)
    use hamiltonian, only: coord_r, angular_r
    integer(i4), intent(in) :: l, i
    real   (dp) :: mat_angular

    mat_angular = angular_r(l, coord_r(i))
end function mat_angular


subroutine check_mat
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
    integer (i4) :: n, i, j 

    n = size(tmp_H(:, 1))
    open(file_check, file = "output/check_mat.d")
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










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    use hamiltonian, only: diag 
    character(30), parameter  :: form_out = '(1A15, 10F10.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: sign, tmp 
    integer  (i4) :: i, j 

    allocate(tmp_H(1:N, 1:N), tmp_E(1:N))
    tmp_H(:, :) = 0_dp
    do j = 1, N 
        do i = 1, j 
            tmp = 0_dp 
            tmp = tmp +mat_Kinet(i, j)
!             tmp = tmp +mat_Kinet(2_i4*N +1_i4 -i, j)
            tmp_H(i, j) = tmp_H(i, j) +tmp 
        enddo
        tmp_H(j, j) = tmp_H(j, j) +mat_Poten(j) +mat_angular(l, j)
    enddo
    call check_mat ! for test 
    call diag('Vectors', 'Upper', tmp_H, tmp_E)

    H(:, :) = 0_dp
    E(:)    = 0_dp 
    sign    = 1_dp 
    do j = 1, N
        if(tmp_H(1, j) > 0_dp) then 
            sign =  1_dp 
        else 
            sign = -1_dp 
        end if 
        H(j, :) = sign*tmp_H(:, j) 
        E(j) = tmp_E(j) 
    end do 
    deallocate(tmp_H, tmp_E) 
    write(file_log, form_out) "Energy: ", (E(i), i = 1, 5) 
end subroutine PROC_H
! end hamiltonian ----------------------------------
! basis plot ---------------------------------------
subroutine PROC_basis_plot(num)
    use hamiltonian, only: coord_rho, coord_r
    integer  (i1), parameter  :: file_psi = 101,           file_ene = 102
    character(30), parameter  :: form_psi = '(30ES25.10)', form_ene = '(1I5, 1ES25.10)'
    integer  (i4), intent(in) :: num 
    integer  (i4) :: i, j, n 
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

    n = size(H(:, 1))
    do i = 1, n 
        write(file_psi, form_psi) dble(i), (H(j, i), j = 1, 10), (H(j, i), j = n, n -9, -1)
!         write(file_psi, form_psi) coord_rho(i), (H(j, i), j = 1, 5), (H(j, i), j = n, n -4, -1)
!         write(file_psi, form_psi) coord_r(i), (H(j, i), j = 1, 5), (H(j, i), j = n, n -4, -1)
        write(file_ene, form_ene) i, E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
! end basis plot -----------------------------------
end module basis
