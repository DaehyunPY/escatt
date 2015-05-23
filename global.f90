module global 
    use kind_type 
    use nr, only: bessj0, bessj1, bessj, bessy0, bessy1, bessy
    implicit none
    real   (RDP), save :: Mass, Kinet, Charge 
    real   (RDP), save :: ra, dr, dtheta
    integer(I4B), save :: N, L, pr, ptheta 
    real   (RDP), save, allocatable :: H(:, :), E(:), R(:), K(:)
    complex(CDP), save, allocatable :: S(:), A(:), inner_u(:, :)
contains 

    function bessel_j(l, x)
        integer(I4B), intent(in) :: l 
        real   (RSP), intent(in) :: x 
        integer(I4B) :: n 
        real   (RDP) :: bessel_j, sign 

        n    = abs(l) 
        sign = 1.d0 
        if(l < 0 .and. mod(-l, 2) == 1) sign = -1.d0 
        if(n == 0) then 
            bessel_j = bessj0(x) 
        else if(n == 1) then 
            bessel_j = bessj1(x)
        else 
            bessel_j = bessj(n, x)
        end if 
        bessel_j = sign*bessel_j 
    end function

    function bessel_y(l, x)
        integer(I4B), intent(in) :: l 
        real   (RSP), intent(in) :: x 
        integer(I4B) :: n 
        real   (RDP) :: bessel_y, sign 

        n    = abs(l)
        sign = 1.d0 
        if(l < 0 .and. mod(-l, 2) == 1) sign = -1.d0 
        if(n == 0) then 
            bessel_y = bessy0(x) 
        else if(n == 1) then 
            bessel_y = bessy1(x)
        else 
            bessel_y = bessy(n, x)
        end if 
        bessel_y = sign*bessel_y 
    end function
end module ! global 
