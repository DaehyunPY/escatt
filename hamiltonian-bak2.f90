module hamiltonian 
    use kind_type
    use global 
    implicit none
    real(dp), save, allocatable, protected :: coord_rho(:) ! coordination 
    real(dp), save, protected :: dtheta, dE 
    real(dp), save, private, protected :: & ! potential 
        Z, alphab, cutoff, &   ! general 
        z0, z1, z2, &          ! II-1 GELTMAN
        pD, pH, pdelta, &      ! II-2 GREEN & YUKAWA
        pA, pB, alpha, beta, & ! II-3 KYLSTRA & BUCKINGHAM
        alpha1, alpha2, alpha3 ! III  CHEN & NOUS 
    integer(i1), save, private, protected :: ty 
contains


! ==================================================
! COORDNATION
! ==================================================
! r ------------------------------------------------
function coord_r(i) 
    integer(i4), intent(in) :: i
    real   (dp) :: coord_r
    coord_r = ra*(-abs(coord_rho(i)) +1_dp)
end function coord_r
! theta --------------------------------------------
function coord_theta(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_theta
    coord_theta = dtheta*dble(i)
end function coord_theta
! E ------------------------------------------------
function coord_E(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_E 
    coord_E = dE*dble(i)
end function coord_E


! ==================================================
! OPERATOR 
! ==================================================
! nabla --------------------------------------------
function nabla_rho(int_j, int_k)
    use math_const, only: pi => math_pi 
    integer(i4), intent(in) :: int_j, int_k
    real   (dp) :: nabla_rho, j, k, qj, qk, M 
    j  = dble(int_j)
    k  = dble(int_k)
    qj = coord_rho(int_j)
    qk = coord_rho(int_k)
    M  = dble(2*N -1)
    if(int_j /= int_k) then 
        nabla_rho = &
            (2_dp*qj/(qj -qk)/(1_dp -qj**2_dp) -2_dp/(qj -qk)**2_dp) &
            *((1_dp -qk**2_dp)/(1_dp -qj**2_dp))**0.5_dp 
    else if(int_j == int_k) then 
        nabla_rho = & 
            (-(M -1_dp)*(M +2_dp) +(M**2_dp +M +6_dp)*qj**2_dp) &
            /3_dp/(1_dp -qj**2_dp)**2_dp  
    end if
end function nabla_rho
! potential ----------------------------------------
function Poten_r(r)
    real(dp), intent(in) :: r
    real(dp) :: Poten_r, stat, pol 
    if(ty == 1) then 
        Poten_r = -Z/r 
    else if(ty == 2) then 
        stat  = -exp(-2_dp*z0*r)*(z1 +z2/r)
        pol   = -alphab/(2_dp*r**4_dp)*(1_dp -exp(-r))**6_dp 
        Poten_r = stat +pol
    else if(ty == 3) then 
        stat  = -(Z/pH)*(exp(-r/pD)/r)*(1_dp +(pH -1_dp)*exp(-r*pH/pD))
        pol   = -alphab/(2_dp*(r**2_dp +pdelta**2_dp)**2_dp)
        Poten_r = stat +pol
    else if(ty == 4) then 
        stat  = -Z*( &
                        +pA**2_dp/(4_dp*alpha**3_dp)*exp(-2_dp*alpha*r)*(alpha +1_dp/r) &
                        +4_dp*pA*pB/(alpha +beta)**3_dp*exp(-(alpha +beta)*r)*((alpha +beta)/2_dp +1_dp/r) & 
                        +pB**2_dp/(4_dp*beta**3_dp)*exp(-2_dp*beta*r)*(beta +1_dp/r))
        pol   = -alphab/(2_dp*(cutoff**2_dp +r**2_dp)**2_dp)
        Poten_r = stat +pol
    else if(ty == 5) then 
        stat  = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol   = -alphab/(2_dp*r**4_dp) * (1_dp -exp(-(r/cutoff)**3_dp))**2_dp 
        Poten_r = stat +pol 
    end if 
end function Poten_r
! angular ------------------------------------------
function angular_r(l, r)
    real   (dp), intent(in) :: r 
    integer(i4), intent(in) :: l 
    real   (dp) :: angular_r 
    angular_r = 1_dp/(2_dp*Mass)*dble(l)*(dble(l) +1_dp)/r**2_dp
end function angular_r


! ==================================================
! CALCULATION 
! ==================================================
! diagonalization ----------------------------------
subroutine diag(jobz, uplo, EF, EV)
    real    (dp), intent(inout) :: EF(1:, 1:)
    real    (dp), intent(out)   :: EV(1:)
    character(1), intent(in)    :: jobz, uplo 
        ! jobz: No vectors, Vectors   uplo: Upper, Lower 
    integer (i8) :: n, lda, lwork
    integer (i1) :: info 
    real    (dp), allocatable :: work(:)

    n     = size(EF(:, 1))
    lda   = size(EF(1, :))
    lwork = 3*n -1
    allocate(work(1:lwork))
    info  = 0

    lwork = -1
    call DSYEV(jobz, uplo, n, EF, lda, EV, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(1:lwork))

    call DSYEV(jobz, uplo, n, EF, lda, EV, work, lwork, info)
    if(info /= 0) stop "SUBROUTINE diag: Error."
    deallocate(work)
end subroutine diag









! ==================================================
! PROCESS
! ==================================================
! input --------------------------------------------
subroutine PROC_input
    use math_const, only: pi => math_pi
    use unit_const, only: other_e_eV, au_hartree
    character(30), parameter :: &
        form_part  = '(4/, 2(45X, 1F15.8, /), /)'
    character(60), parameter :: & 
        form_poten = '(4/, 1(45X, 1I15, /), 3(45X, 1F15.8, /))'
    character(30), parameter :: &
        form_p1    = '(18/)', &
        form_p2    = '(  /, 3(45X, 1F15.8, /), 14/)', &
        form_p3    = '( 5/, 3(45X, 1F15.8, /), 10/)', &
        form_p4    = '( 9/, 4(45X, 1F15.8, /),  5/)', &
        form_p5    = '(14/, 3(45X, 1F15.8, /),   /)'
    character(60), parameter :: &
        form_cal   = '(4/, 1(45X, 1F15.8, /), 5(45X, 1I15, /), /)'
    character(90), parameter :: &
        form_opt   = '(6/, 2(45X, 6X, 1A1, /), /, 5(45X, 6X, 1A1, /), 3/, 3(45X, 6X, 1A1, /))'
    real     (dp), parameter :: eV_to_au = other_e_ev/au_hartree
    real     (dp) :: unit_e

    open(file_input, file = "input.d")
    read(file_input, form_part)  Mass, Scatt
    read(file_input, form_poten) ty, Z, alphab, cutoff
    if(ty == 1) then 
        read(file_input, form_p1) 
    else if(ty == 2) then 
        read(file_input, form_p2) z0, z1, z2
    else if(ty == 3) then 
        read(file_input, form_p3) pD, pH, pdelta
        if(pD < 0_dp) then 
            pD = pH/Z**0.4_dp 
        end if 
        if(pH < 0_dp)  then 
            pH = pD*Z**0.4_dp 
        end if 
        if(pD < 0_dp .and. pH < 0_dp) stop "SUBROUTINE PROC_input: Check potential coefficinet."
        if(pdelta < 0_dp) then 
            pdelta = (alphab/(2_dp*Z**(1_dp/3_dp)))**(0.25_dp)
        end if 
    else if(ty == 4) then 
        read(file_input, form_p4) pA, pB, alpha, beta
    else if(ty == 5) then 
        read(file_input, form_p5) alpha1, alpha2, alpha3
    else 
        stop "SUBROUTINE PROC_input: Check potential type."
    end if 
    read (file_input, form_cal) ra, N, M, L, pr, ptheta
    read (file_input, form_opt) op_ev, op_degree, op_poten, op_basis, op_dcs, op_inner, op_outer, op_tcs, op_phase, op_lt 
    close(file_input) 
    open (file_log, file = "output/log.d")

    unit_e = 1_dp 
    if(op_ev == "Y") unit_e = eV_to_au
    
    if(pr > N) pr = N 
    Scatt  = Scatt*unit_e
    dtheta = pi/dble(ptheta)
    dE     = Scatt/dble(M) 

    allocate(coord_rho(1:2*N -1))
!     allocate(H(1:N, 1:N), E(1:N))
    allocate(H(1:2*N -1, 1:2*N -1), E(1:2*N -1)) ! for test 
    allocate(R(0:L), K(0:L), S(0:L), A(0:L))
    H(:, :)       = 0_dp
    E(:)          = 0_dp
    R(:)          = 0_dp
    K(:)          = 0_dp
    S(:)          = 0_dp

    if(op_tcs == "Y" .or. op_phase == "Y" .or. op_lt == "Y") then 
        op_basis = "N"
        op_dcs   = "N"
        op_inner = "N" 
        op_outer = "N"
    else if(op_tcs == "N" .and. op_phase == "N" .and. op_lt == "N") then 
        M  = 1
        dE = Scatt
    else 
        stop "SUBROUTINE PROC_input: Check option type."
    end if 
end subroutine PROC_input
! end input ----------------------------------------
! information --------------------------------------
subroutine PROC_inform
    use unit_const, only: other_e_eV, au_hartree
    character(30), parameter :: form_out = '(1A15, 1ES15.3)'
    real     (dp), parameter :: au_to_eV = au_hartree/other_e_ev
    
    write(file_log, *) "================================================================="
    write(file_log, *) "PARTICLE: ELECTRON"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " MASS                                   [au] ", Mass 
    write(file_log, *) " KINETIC ENERGY                         [au] ", Scatt
    write(file_log, *) "                                        [eV] ", Scatt*au_to_eV
    write(file_log, *) " - "
    write(file_log, *) " - "

    if(ty == 1) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: COULOMB"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 2) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: GELTMAN"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " COEFFICIENT z0                         [au] ", z0
        write(file_log, *) " COEFFICIENT z1                         [au] ", z1
        write(file_log, *) " COEFFICIENT z2                         [au] ", z2
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 3) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: GREEN & YUKAWA"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " COEFFICIENT D                          [au] ", pD 
        write(file_log, *) " COEFFICIENT H                          [au] ", pH 
        write(file_log, *) " COEFFICIENT delta                      [au] ", pdelta
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 4) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: KYLSTRA & BUCKINGHAM"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " CUTOFF                                 [au] ", cutoff
        write(file_log, *) " COEFFICIENT A                          [au] ", pA 
        write(file_log, *) " COEFFICIENT B                          [au] ", pB 
        write(file_log, *) " COEFFICIENT alpha                      [au] ", alpha 
        write(file_log, *) " COEFFICIENT beta                       [au] ", beta 
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 5) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: CHEN & NOUS"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " CUTOFF                                 [au] ", cutoff
        write(file_log, *) " COEFFICIENT alpha1                     [au] ", alpha1
        write(file_log, *) " COEFFICIENT alpha2                     [au] ", alpha2
        write(file_log, *) " COEFFICIENT alpha3                     [au] ", alpha3 
        write(file_log, *) " - "
        write(file_log, *) " - "
    end if 

    write(file_log, *) "================================================================="
    write(file_log, *) "UNIT"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " DEFAULT                                               au        "
    if(op_ev     == "Y") write(file_log, *) " ENERGY                                                eV        "
    if(op_degree == "Y") write(file_log, *) " ANGULAR                                           degree        "
    write(file_log, *) " - "
    write(file_log, *) " - "

    write(file_log, *) "================================================================="
    write(file_log, *) "CALCULATION"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " BOUNDARY SIZE                          [au] ", ra 
    write(file_log, *) " GRID NUMBER OF r COORDINATES            [1] ", N 
    write(file_log, *) " MAXIUM OF QUANTUM NUMBER L              [1] ", L 
    write(file_log, *) " - "
    write(file_log, *) " - "
end subroutine PROC_inform
! end information ----------------------------------
! coordination -------------------------------------
subroutine PORC_coord
    real   (dp), allocatable :: X(:, :)
    real   (dp) :: tmp 
    integer(i4) :: i
    
    allocate(X(1:2*N -1, 1:2*N -1))
    X = 0_dp 
    do i = 1, 2*N -2
        tmp = dble(i)
        X(i, i +1) = tmp/((2_dp*tmp -1_dp)*(2_dp*tmp +1_dp))**0.5_dp 
    end do 
    call diag('No vectors', 'Upper', X, coord_rho)
    deallocate(X)
end subroutine PORC_coord
! end coordination ---------------------------------
! potential plot -----------------------------------
subroutine PROC_Poten_plot
    integer  (i1), parameter :: file_poten = 101
    character(30), parameter :: form_poten = '(30ES25.10)'
    integer  (i4) :: i 

    open(file_poten, file = "output/poten.d")
    do i = 1, N, N/pr 
        write(file_poten, form_poten) coord_r(i), Poten_r(coord_r(i))
    end do 
    close(file_poten)
end subroutine PROC_Poten_plot
! end potential plot -------------------------------
! out ----------------------------------------------
subroutine PROC_out 
    deallocate(coord_rho)
    deallocate(H, E)
    deallocate(R, K, S, A)
    close(file_log)
end subroutine PROC_out 
! end out ------------------------------------------
end module hamiltonian 
