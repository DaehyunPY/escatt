module global 
    use kind_type 
    implicit none
    real   (RDP), save :: Mass, Kinet, Charge 
    real   (RDP), save :: ra, dr, dtheta
    integer(I4B), save :: N, L, pr, ptheta 
    real   (RDP), save, allocatable :: H(:, :), E(:), R(:), K(:)
    complex(CDP), save, allocatable :: S(:), inner_u(:, :)

    interface
        function plgndr_s(l,m,x)
            use kind_type
            integer(I4B), intent(in) :: l,m
            real(RDP), intent(in) :: x
            real(RDP) :: plgndr_s
        end function ! plgndr_s
    end interface
end module ! global 
