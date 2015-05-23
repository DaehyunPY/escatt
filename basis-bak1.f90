module basis_function
    implicit none
    real(8), parameter :: r_a = 20.d0 
    integer, parameter ::   N = 2000
    real(8), parameter ::  dr = r_a/dble(N)
    real(8) :: H(1:N, 1:N), E(1:N)
contains


function Poten(r)
    real(8), intent(in) :: r
    real(8) :: Poten
    Poten = 2.d0/r
end function Poten
function cood_r(i)
    integer, intent(in) :: i
    real(8) :: cood_r
    cood_r = dr*dble(i)
end function cood_r


function nabla_x(dr, ty_i, ty_j)
    real(8), intent(in) :: dr
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
function nabla_r(ty_i, ty_j)
    integer, intent(in) :: ty_i, ty_j
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    real(8) :: nabla_r, i, j
    i = dble(ty_i)
    j = dble(ty_j)
    if( ty_i /= ty_j ) then 
        nabla_r = &
            (pi**2*((i-j)**3*(j+i)**2*sin(pi*(j+i))-sin(pi*(i-j))*(i-j)**2*(j+i)**3) &
            -2*(i-j)**3*sin(pi*(j+i))+2*pi*(i-j)**3*(j+i)*cos(pi*(j+i)) &
            -2*pi*cos(pi*(i-j))*(i-j)*(j+i)**3+2*sin(pi*(i-j))*(j+i)**3) &
            /(pi*dr**2*(i-j)**3*(j+i)**3)
    else
        nabla_r = &
            (3*pi**2*(j+i)**2*sin(pi*(j+i))-6*sin(pi*(j+i))+6*pi*(j+i)*cos(pi*(j+i)) &
            -pi**3*(j+i)**3)/(pi*dr**2*(j+i)**3)/3.d+0
    end if
end function nabla_r
end module basis_function










module basis_process
    use basis_function
    implicit none
contains


subroutine PROC_set_nomal_H
    real(8), allocatable :: tmp_H(:, :), tmp_E(:)
    integer, parameter   :: file_Psi = 101, file_Energy = 102  
    real(8) :: sign 
    integer :: i, j

    allocate(tmp_H(1:2*N, 1:2*N), tmp_E(1:2*N))
    tmp_H(:, :) = 0.d0  
    do j = 1, 2*N
        do i = 1, 2*N
            tmp_H(i, j) = nabla_x(dr, i, j)
        enddo
    enddo
    do i = 1, 2*N
        tmp_H(i, i) = tmp_H(i, i) +Poten(cood_r(i))
    enddo
    tmp_H(:, :) = -0.5d0*tmp_H(:, :)
!     write(*, *) "Here. (1)"   

    call PROC_diag(tmp_H, tmp_E)
!     write(*, *) "Here. (2)"   

    do j = 1, N 
        E(j) = tmp_E(j)
        sign = 1.d0 
        if(tmp_H(1, j) < 0.d0) then 
            sign = -1.d0 
        end if 
        H(j, :) = sign*tmp_H(1:N, j)
    end do 
    deallocate(tmp_H, tmp_E)
!     write(*, *) "Here. (3)"   

    open(file_Psi, file = "inout/nomal_psi.d")
    open(file_Energy, file = "inout/nomal_energy.d")
    call PROC_plot(file_Psi, file_Energy)
    close(file_Psi)
    close(file_Energy)
!     write(*, *) "Here. (4)"   
end subroutine PROC_set_nomal_H


subroutine PROC_set_H
    real(8), allocatable :: tmp_H(:, :), tmp_E(:)
    integer, parameter   :: file_Psi = 101, file_Energy = 102 
!     real(8) :: sign 
    real(16) :: sign 
    integer :: i, j, k, num 

!     write(*, *) "Here. (1-1)"   
    allocate(tmp_H(2*N, 2*N), tmp_E(2*N))
!     write(*, *) "Here. (1-2)"   
    tmp_H(:, :) = 0.d0  
    do j = 1, 2*N
        do i = 1, 2*N
            tmp_H(i, j) = nabla_x(dr, i, j)
        enddo
    enddo
!     write(*, *) "Here. (1-3)"   
    do i = 1, N
        tmp_H(i, i) = tmp_H(i, i) +Poten(cood_r(i))
        tmp_H(2*N +1 -i, 2*N +1 -i) = tmp_H(2*N +1 -i, 2*N +1 -i) +Poten(cood_r(i))
    enddo
!     write(*, *) "Here. (1-4)"   
    tmp_H(:, :) = -0.5d0*tmp_H(:, :)
!     write(*, *) "Here. (1)"   

    call PROC_diag(tmp_H, tmp_E)
!     write(*, *) "Here. (2)"   

    num = 0 
    do j = 1, 2*N
        sign = tmp_H(1, j)*tmp_H(2*N -1, j)
        if(sign > 0.d0) then 
            num = num +1 
            E(num) = tmp_E(j)
            sign = 1.d0 
            if(tmp_H(1, j) < 0.d0) then 
                sign = -1.d0 
            end if 
            do i = 1, N 
                H(num, i) = sign*(tmp_H(i, j) +tmp_H(2*N -i, j))*2.d0**(-0.5d0)
            end do 
        end if 
    end do 
!     do k = 1, N
!         if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) > 0.d0 & 
!             .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) < 0.d0) then 
!                 j = 2*k -1 
!                 num = num +1
!         else if(tmp_H(1, 2*k -1)*tmp_H(2*N -1, 2*k -1) < 0.d0 & 
!             .and. tmp_H(1, 2*k)*tmp_H(2*N -1, 2*k) > 0.d0) then 
!                 j = 2*k 
!                 num = num +1
!         else
!                 j = 2*k -1
!                 num = num +1
!         end if 
!         E(num) = tmp_E(j)
!         sign = 1.d0 
!         if(tmp_H(1, j) < 0.d0) then 
!             sign = -1.d0 
!         end if 
!         do i = 1, N 
!             H(num, i) = sign*(tmp_H(i, j) +tmp_H(2*N -i, j))*2.d0**(-0.5d0)
!         end do 
!     end do 
!     write(*, *) "Here. (3-1)", num, N   
    if(num /= N) stop "PROCESS set H: Something is wrong." 
    deallocate(tmp_H, tmp_E)
!     write(*, *) "Here. (3)"   

    open(file_Psi, file = "inout/basis_psi.d")
    open(file_Energy, file = "inout/basis_energy.d")
    call PROC_plot(file_Psi, file_Energy)
    close(file_Psi)
    close(file_Energy)
!     write(*, *) "Here. (4)"   
end subroutine PROC_set_H


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


subroutine PROC_plot(file_Psi, file_Energy)
    integer,       intent(in) :: file_Psi, file_Energy
    character(30), parameter  :: form_Psi = '(30ES25.10)', form_Energy = '(ES25.10)', form_out = '(10F10.3)'
    integer :: i, j 
    write(file_Psi, form_Psi) 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    do i = 1, N 
        write(file_Psi, form_Psi) cood_r(i), (H(j, i), j = 1, 5)
        write(file_Energy, form_Energy) E(i)
    end do 
    write(*, form_out) (E(i), i = 1, 10)
end subroutine PROC_plot
end module basis_process
