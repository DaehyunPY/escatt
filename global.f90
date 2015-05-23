module global 
    implicit none
    real(8), parameter :: ra    = 20.d0
    integer, parameter :: N     = 2000
    integer, parameter :: L     = 0
    real(8), parameter :: dr    = ra/dble(N)
    real(8), parameter :: Mass  = 1.d0
    real(8), parameter :: Kinet = 1.d0
    real(8), save :: H(1:N, 1:N), E(1:N)
    real(8), save :: R(0:L), K(1:L)
    complex(8), save :: S(0:L)
contains


function Poten(r)
    real(8), intent(in) :: r
    real(8) :: Poten

    Poten = -1.d0/r
end function Poten


function cood_r(i)
    integer, intent(in) :: i
    real(8) :: cood_r

    cood_r = dr*dble(i)
end function cood_r
end module global 










! module PROC_global
!     implicit none
! contains


!     subroutine PROC_input
!         argument type, intent(inout) :: 
        
!     end subroutine PROC_input
! end module PROC_global
