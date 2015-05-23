module boundary
    use kind_type
    use global 
    implicit none
contains 


subroutine mat_R(l)
    character(30), parameter  :: form_out = '(1A15, 1I15, 1ES15.3)'
    integer (I4B), intent(in) :: l 
    real    (RQP) :: sum 
    integer (I4B) :: i, j 

    sum = 0.d0 
    do i = 1, N 
        sum = sum +H(i, N)**2.d0/(2.d0*Mass*(E(i) -Kinet)) 
    end do 
    R(l) = sum/ra 
    write(*, form_out) "R: ", l, R(l)
end subroutine mat_R


! subroutine mat_K(l)
!     character(30), parameter  :: form_out = '(1A15, 1I15, 1ES15.3)'
!     integer (I4B), intent(in) :: l 
!     real    (RDP) :: ka, agamma, diff_j, diff_y, tmp1, tmp2 

!     ka     = (2.d0*Mass*Kinet)**0.5d0*ra
!     agamma = 1.d0/R(l) -1.d0
!     diff_j = (bessel_jn(l -1, ka) -bessel_jn(l +1, ka))/2.d0
!     diff_y = (bessel_yn(l -1, ka) -bessel_yn(l +1, ka))/2.d0 
!     tmp1   = ka*diff_j -agamma*bessel_jn(l, ka)
!     tmp2   = ka*diff_y -agamma*bessel_yn(l, ka)
!     K(l)   = tmp1/tmp2
!     write(*, form_out) "K: ", l, K(l)
! end subroutine mat_K
subroutine mat_K(l)
    use nr, only: sphbes_s
    character(30), parameter  :: form_out = '(1A15, 1I15, 1ES15.3)'
    integer (I4B), intent(in) :: l 
    real    (RSP) :: ka, sb_j, sb_y, diff_j, diff_y 
    real    (RDP) :: agamma, tmp1, tmp2 

    ka = (2.d0*Mass*Kinet)**0.5d0*ra
    call sphbes_s(l, ka, sb_j, sb_y, diff_j, diff_y)

    agamma = 1.d0/R(l) -1.d0
    tmp1   = ka*diff_j -agamma*sb_j 
    tmp2   = ka*diff_y -agamma*sb_y 
    K(l)   = tmp1/tmp2
    write(*, form_out) "K: ", l, K(l)
end subroutine mat_K


subroutine mat_S(l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1I15, 1ES15.3, 1ES15.3, "i")'
    integer (I4B), intent(in) :: l 

    S(l) = (1.d0 +i*K(l))/(1.d0 -i*K(l))
    write(*, form_out) "S: ", l, S(l)
end subroutine mat_S


subroutine mat_A(l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1I15, 1ES15.3, 1ES15.3, "i")'
    integer (I4B), intent(in) :: l 
    complex (CDP) :: sign 
    real    (RDP) :: tmp1, tmp2 

    if(mod(l, 4) == 0) then 
        sign = 1.d0 
    else if(mod(l, 4) == 1) then 
        sign = i 
    else if(mod(l, 4) == 2) then 
        sign = -1.d0 
    else if(mod(l, 4) == 3) then 
        sign = -i 
    end if 
    tmp1 = (1.d0 +real(S(l)))/2.d0 
    tmp2 = aimag(S(l))/2.d0 
    A(l) = sign*(2.d0*dble(l) +1.d0)*exp(tmp1 +i*tmp2)
    write(*, form_out) "A: ", l, A(l)
end subroutine mat_A












! ==================================================
! PROCESS
! ==================================================


subroutine PROC_boundary_mat(l) ! It must be called after PROC_input, PROC_H 
    integer(I4B), intent(in) :: l 

    call mat_R(l)
    call mat_K(l)
    call mat_S(l)  
    call mat_A(l)
end subroutine PROC_boundary_mat
end module boundary
