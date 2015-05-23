module global 
    use kind_type 
    implicit none
    real   (RDP), save :: Mass, Kinet, Charge 
    real   (RDP), save :: ra, dr, dtheta
    integer(I4B), save :: N, L, pr, ptheta 
    real   (RDP), save, allocatable :: H(:, :), E(:), R(:), K(:)
    complex(CDP), save, allocatable :: S(:), A(:), inner_u(:, :)
end module ! global 
