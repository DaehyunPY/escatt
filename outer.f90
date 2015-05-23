module outer
    use global 
    implicit none
    complex(8), save, allocatable, protected :: f(:)
contains 
    

subroutine mat_f 
    complex(8), parameter :: i  = (0.d0, 1.d0)
    real(8),    parameter :: pi = 2.0d0*acos(0.0d0)
    real(8) :: k 
    integer :: j 

    k = (2.d0*Mass*Kinet)**0.50
    do j = 0, L 
        f(j) = 1.d0/(2.d0*i*k)*(2.d0*dble(j) +1.d0)*(S(j) -1.d0)
    end do 
end subroutine mat_f 


function outer_u(l, r)
    integer,    intent(in) :: l 
    real(8),    intent(in) :: r 
    real(8),    parameter  :: pi = 2.0d0*acos(0.0d0)
    complex(8), parameter  :: i  = (0.d0, 1.d0)
    real(8)    :: k, sckr 
    complex(8) :: outer_u

    k       = (2.d0*Mass*Kinet)**0.5d0
    sckr    = k*r -dble(l)*pi/2.d0
    outer_u = (exp(-i*sckr) -S(l)*exp(i*sckr))/k 
end function outer_u










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_CS_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer,       parameter :: file_dcs = 101
    character(30), parameter :: form_dcs = '(30ES25.10)'
    character(30), parameter :: form_out = '(1A15, 5X, 1ES25.10)'
    real(8), parameter :: pi = 2.0d0*acos(0.0d0) 
    real(8), parameter :: radian_to_degree = 180.d0/pi 
    complex(16) :: tmp1, tmp2  
    real(8)     :: k 
    integer     :: i, j 

    k = (2.d0*Mass*Kinet)**0.50
    allocate(f(0:L))
    
    call mat_f 

    tmp1 = 0.d0 
    do i = 0, L 
        tmp1 = tmp1 +f(i)
!         write(*, *) tmp1 
    end do 
    tmp1 = 4.d0*pi/k*tmp1
    write(*, form_out) "total sigma: ", dble(aimag(tmp1))

    open(file_dcs, file = "inout/dcs.d")
    tmp2 = 0.d0 
    do j = 0, ptheta 
        do i = 0, L 
            tmp2 = tmp2 +f(i)*plgndr_s(i, 0, cos(coord_theta(j)))
        end do 
    write(file_dcs, form_dcs) coord_theta(j)*radian_to_degree, abs(tmp2)**2.d0
!     write(file_dcs, form_dcs) coord_theta(j)*radian_to_degree, abs(tmp2/tmp1)**2.d0
    end do 
    close(file_dcs)
    deallocate(f)
end subroutine PROC_CS_plot


subroutine PROC_outer_plot ! It must be called after PROC_input, PROC_H, PROC_boundary
    integer,       parameter :: file_psi1 = 101, file_psi2 = 102, file_psi3 = 103
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real(8), parameter :: pi = 2.0d0*acos(0.0d0) 
    real(8), parameter :: radian_to_degree = 180.d0/pi 
    complex(16) :: tmp 
    integer :: i, j, k 

    open(file_psi1, file = "inout/outer_u0.d")
    tmp = 0.d0 
    do i = N +1, 2*N
!         tmp = outer_u(0, coord_r(i))/coord_r(i)
        tmp = outer_u(0, coord_r(i))
        write(file_psi1, form_psi) coord_r(i), dble(abs(tmp)**2.d0)
!         write(file_psi1, form_psi) coord_r(i), dble(real(tmp))
!         write(file_psi1, form_psi) coord_r(i), dble(aimag(tmp))
    end do 
    close(file_psi1)

    open(file_psi2, file = "inout/outer_psi.d")
    open(file_psi3, file = "inout/outer_psi-4pir.d")
    do i = N +1, 2*N, N/pr 
        do j = 0, ptheta
            tmp = 0.d0 
            do k = 0, L 
                tmp = tmp +outer_u(k, coord_r(i))/coord_r(i)*plgndr_s(i, 0, cos(coord_theta(j)))
            end do 
!             write(file_psi2, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(tmp))**2.d0 
            write(file_psi2, form_psi) coord_r(i), coord_theta(j), dble(abs(tmp))**2.d0 
            if(j == 0) write(file_psi2, form_psi) coord_r(i), dble(abs(tmp)**2.d0)
            tmp = tmp*4.d0*pi*coord_r(i)
!             write(file_psi3, form_psi) coord_r(i), coord_theta(j)*radian_to_degree, dble(abs(tmp))**2.d0
            write(file_psi3, form_psi) coord_r(i), coord_theta(j), dble(abs(tmp))**2.d0
        end do 
        write(file_psi2, form_psi) 
        write(file_psi3, form_psi) 
    end do 
    close(file_psi2)
    close(file_psi3)
end subroutine PROC_outer_plot
end module outer
