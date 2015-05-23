module hamiltonian 
    use kind_type
    use global 
    implicit none
contains


function coord_r(i)
    integer(I4B), intent(in) :: i
    real(RDP) :: coord_r

    coord_r = dr*dble(i)
end function coord_r


function coord_theta(i)
    integer(I4B), intent(in) :: i
    real(RDP) :: coord_theta

    coord_theta = dtheta*dble(i)
end function coord_theta


function nabla_x(int_i, int_j)
    use math_const, only: pi => math_pi 
    integer(I4B), intent(in) :: int_i, int_j
    real   (RDP) :: nabla_x, i, j

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


function Poten(r)
    real(RDP), intent(in) :: r
    real(RDP) :: Poten

    Poten = Charge*r**2.d0*exp(-r)
end function Poten










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_input
    use math_const, only: pi => math_pi
    integer (I1B), parameter :: file_input       = 101
    character(30), parameter :: form_particle    = '(////, 2(40X, 1F10.4, /), /)'
    character(30), parameter :: form_potential   = '(////, 1(40X, 1F10.4, /), /)'
    character(60), parameter :: form_calculation = '(////, 1(40X, 1F10.4, /), 2(40X, 1I10, /), /)'
    character(30), parameter :: form_plot        = '(////, 2(40X, 1I10, /), /)'
!     real(RDP),     parameter :: pi = 2.0d0*acos(0.0d0) 

    open(file_input, file = "inout/input.d")
    read(file_input, form_particle)    Mass, Kinet
    read(file_input, form_potential)   Charge
    read(file_input, form_calculation) ra, N, L 
    read(file_input, form_plot)        pr, ptheta 
    close(file_input) 

    dr     = ra/dble(N) 
    dtheta = pi/dble(ptheta)
!     dtheta = (2.d0*pi)/dble(ptheta)

    allocate(H(1:N, 1:N), E(1:N))
    allocate(R(0:L), K(0:L), S(0:L))
    allocate(inner_u(0:L, 1:N))
    H(:, :)       = 0.d0
    E(:)          = 0.d0
    R(:)          = 0.d0
    K(:)          = 0.d0
    S(:)          = 0.d0
    inner_u(:, :) = 0.d0
end subroutine PROC_input


subroutine PROC_Poten_plot
    integer (I1B), parameter :: file_poten = 101
    character(30), parameter :: form_poten = '(30ES25.10)'
    integer (I4B) :: i 

    open(file_poten, file = "inout/poten.d")
    do i = 1, N, N/pr 
        write(file_poten, form_poten) coord_r(i), Poten(coord_r(i))
    end do 
    close(file_poten)
end subroutine PROC_Poten_plot


subroutine PROC_out 
    deallocate(H, E)
    deallocate(R, K, S)
    deallocate(inner_u)
end subroutine PROC_out 
end module hamiltonian 