module global 
    implicit none
    real(8), parameter ::   r_a = 20.d0
    integer, parameter ::     N = 1000
    real(8), parameter ::    dr = r_a/dble(N)
    real(8), parameter ::  Mass = 1.d0
    real(8), parameter :: Kinet = 1.d0 
    real(8), save :: H(1:N, 1:N), E(1:N)
    real(8), save :: R, K
    complex(8), save :: S
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
