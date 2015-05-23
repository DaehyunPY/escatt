module hamiltonian 
    use kind_type
    use global 
    implicit none
    real(dp), save, protected :: dr, dtheta, dE 
    real(dp), save, private, protected :: & 
        Z, alphab, cutoff, &
        z0, z1, z2, &          ! II-1 potential GELTMAN
        pD, pH, pdelta, &      ! II-2 potential GREEN & YUKAWA
        pA, pB, alpha, beta, & ! II-3 potential KYLSTRA & BUCKINGHAM
        alpha1, alpha2, alpha3 ! III  potential CHEN & NOUS 
    integer(i1),  save, private, protected :: ty 
contains


function coord_r(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_r
    coord_r = dr*dble(i)
end function coord_r
function weight_r(i)
    integer(i4), intent(in) :: i 
    real   (dp) :: weight_r, q, alpha 
    q        = coord_r(i)
    alpha    = 0_dp
    weight_r = q**(alpha)*exp(-q)
end function weight_r
function coord_theta(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_theta
    coord_theta = dtheta*dble(i)
end function coord_theta
function coord_E(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_E 
    coord_E = dE*dble(i)
end function coord_E


! function nabla(int_j, int_k) ! dvr 
!     use math_const, only: pi => math_pi 
!     integer(i4), intent(in) :: int_j, int_k
!     real   (dp) :: nabla, j, k, qj, qk  
!     j = dble(int_j)
!     k = dble(int_k)
!     qj = coord_r(j)
!     qk = coord_r(k)
!     if(int_j /= int_k) then 
!         nabla = &
!             (2_dp*qj/(qj -qk)/(1_dp -qj**2_dp) -2_dp/(qj -qk)**2_dp) &
!             *((1_dp -qk**2_dp)/(1_dp -qj**2_dp))**0.5_dp 
!     else if(int_j == int_k) then 
!         nabla = & 
!             (-(N -1_dp)*(N +2_dp) +(N**2_dp +N +6_dp)*qj**2_dp) &
!             /3_dp/(1_dp -qj**2_dp)**2_dp  
!     end if
! end function nabla
! function nabla(int_i, int_j) ! nomal r 
!     use math_const, only: pi => math_pi 
!     integer(i4), intent(in) :: int_i, int_j
!     real   (dp) :: nabla, i, j
!     i = dble(int_i)
!     j = dble(int_j)
!     if(int_i /= int_j) then 
!         nabla = &
!             2_dp/(i -j)**2_dp -2_dp/(i -j)**2_dp 
!     else
!         nabla = &
!             pi**2_dp/3_dp -0.5_dp/i**2_dp 
!     end if
!     nabla = nabla/dr**2_dp
!     if(mod(int_i +int_j, 2) == 0) then 
!         nabla = -nabla
!     end if 
! end function nabla
function nabla(int_i, int_j) ! nomal x
    use math_const, only: pi => math_pi 
    integer(i4), intent(in) :: int_i, int_j
    real   (dp) :: nabla, i, j
    i = dble(int_i)
    j = dble(int_j)
    if(int_i /= int_j) then 
        nabla = &
            2_dp/(i -j)**2_dp
    else
        nabla = &
            pi**2_dp /3_dp
    end if
    nabla = nabla/dr**2_dp
    if(mod(int_i +int_j, 2) == 0) then 
        nabla = -nabla
    end if 
end function nabla
function Poten(r) ! potential 
    real(dp), intent(in) :: r
    real(dp) :: Poten, stat, pol 
    if(ty == 1) then 
        Poten = -Z/r 
    else if(ty == 2) then 
        stat  = -exp(-2_dp*z0*r)*(z1 +z2/r)
        pol   = -alphab/(2_dp*r**4_dp)*(1_dp -exp(-r))**6_dp 
        Poten = stat +pol
    else if(ty == 3) then 
        stat  = -(Z/pH)*(exp(-r/pD)/r)*(1_dp +(pH -1_dp)*exp(-r*pH/pD))
        pol   = -alphab/(2_dp*(r**2_dp +pdelta**2_dp)**2_dp)
        Poten = stat +pol
    else if(ty == 4) then 
        stat  = -Z*( &
                        +pA**2_dp/(4_dp*alpha**3_dp)*exp(-2_dp*alpha*r)*(alpha +1_dp/r) &
                        +4_dp*pA*pB/(alpha +beta)**3_dp*exp(-(alpha +beta)*r)*((alpha +beta)/2_dp +1_dp/r) & 
                        +pB**2_dp/(4_dp*beta**3_dp)*exp(-2_dp*beta*r)*(beta +1_dp/r))
        pol   = -alphab/(2_dp*(cutoff**2_dp +r**2_dp)**2_dp)
        Poten = stat +pol
    else if(ty == 5) then 
        stat  = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol   = -alphab/(2_dp*r**4_dp) * (1_dp -exp(-(r/cutoff)**3_dp))**2_dp 
        Poten = stat +pol 
    end if 
end function Poten
function angular(l, r) ! angular 
    real   (dp), intent(in) :: r 
    integer(i4), intent(in) :: l 
    real   (dp) :: angular 
    angular = 1_dp/(2_dp*Mass)*dble(l)*(dble(l) +1_dp)/r**2_dp
end function angular










! ==================================================
! PROCESS
! ==================================================


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
    dr     = ra/dble(N) 
    dtheta = pi/dble(ptheta)
    dE     = Scatt/dble(M) 

    allocate(H(1:N, 1:N), E(1:N))
    allocate(R(0:L), K(0:L), S(0:L), A(0:L))
    allocate(inner_u(0:L, 1:N)) 
    H(:, :)       = 0_dp
    E(:)          = 0_dp
    R(:)          = 0_dp
    K(:)          = 0_dp
    S(:)          = 0_dp
    inner_u(:, :) = 0_dp

    if(op_tcs == "Y" .or. op_phase == "Y" .or. op_lt == "Y") then 
        op_basis = "N"
        op_dcs   = "N"
        op_inner = "N" 
        op_outer = "N"
        allocate(CS(1:M))
        CS(:) = 0_dp 
    else if(op_tcs == "N" .and. op_phase == "N" .and. op_lt == "N") then 
        M  = 1
        dE = Scatt
    else 
        stop "SUBROUTINE PROC_input: Check option type."
    end if 
end subroutine PROC_input


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


subroutine PROC_Poten_plot
    integer  (i1), parameter :: file_poten = 101
    character(30), parameter :: form_poten = '(30ES25.10)'
    integer  (i4) :: i 

    open(file_poten, file = "output/poten.d")
    do i = 1, N, N/pr 
        write(file_poten, form_poten) coord_r(i), Poten(coord_r(i))
    end do 
    close(file_poten)
end subroutine PROC_Poten_plot


subroutine PROC_out 
    deallocate(H, E)
    deallocate(R, K, S, A)
    deallocate(inner_u)
    close(file_log)
end subroutine PROC_out 
end module hamiltonian 
