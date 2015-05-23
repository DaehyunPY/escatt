module global 
    use kind_type 
   implicit none
    integer  (I1), parameter :: file_input = 5, file_log = 6 
    real     (DP), save :: Mass, Kinet, Charge, ra 
    integer  (I4), save :: N, M, L, pr, ptheta 
    real     (DP), save, allocatable :: H(:, :), E(:), R(:), K(:), CS(:) 
    complex  (DP), save, allocatable :: S(:), A(:), inner_u(:, :)
    character(1),  save :: &
        op_poten, op_basis, op_dcs, op_inner, op_outer, &
        op_tcs, op_phase, op_lt 
end module ! global 
