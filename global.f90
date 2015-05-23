module global 
    implicit none
    real(8), save :: ra, dr, Mass, Kinet, Charge 
    integer, save :: N, L
    real(8),    allocatable, save :: H(:, :), E(:), R(:), K(:)
    complex(8), allocatable, save :: S(:), inner_u(:, :)
contains


function Poten(r)
    real(8), intent(in) :: r
    real(8) :: Poten

    Poten = -Charge/r
end function Poten


function coord_r(i)
    integer, intent(in) :: i
    real(8) :: coord_r

    coord_r = dr*dble(i)
end function coord_r
end module global 










module PROC_global
    use global
    implicit none
contains


subroutine PROC_input
    integer,       parameter :: file_input       = 101
    character(50), parameter :: form_calculation = '(////, 1(40X, 1F, /), 2(40X, 1I, /), /)'
    character(30), parameter :: form_particle    = '(////, 2(40X, 1F, /), /)'
    character(30), parameter :: form_potential   = '(////, 1(40X, 1F, /), /)'

    open(file_input, file = "inout/input.d")
    read(file_input, form_calculation) ra, N, L 
    read(file_input, form_particle)    Mass, Kinet
    read(file_input, form_potential)   Charge
    close(file_input) 

    dr = ra/dble(N) 

    allocate(H(1:N, 1:N), E(1:N))
    allocate(R(0:L), K(0:L), S(0:L))
    allocate(inner_u(0:L, 1:N))
    inner_u(:, :) = 0.d0 
end subroutine PROC_input


subroutine PROC_out 
    deallocate(H, E)
    deallocate(R, K, S)
    deallocate(inner_u)
end subroutine PROC_out 
end module PROC_global
