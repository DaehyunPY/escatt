module global 
    use kind_type 
    implicit none
    integer  (I1), parameter :: file_input = 5, file_log = 6 
    real     (DP), save :: Mass, Kinet, Charge 
    real     (DP), save :: ra, dr, dtheta
    integer  (I4), save :: N, L, pr, ptheta 
    real     (DP), save, allocatable :: H(:, :), E(:), R(:), K(:)
    complex  (DP), save, allocatable :: S(:), A(:), inner_u(:, :)
end module ! global 
