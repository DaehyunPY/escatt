module hamiltonian 
    use kind_type
    use global 
    implicit none
    real(RDP), save, private, protected :: & 
        Z, alphab, &
        z0, z1, z2, &               ! first potential
        pH, pD, pdelta, &           ! second potential
        pA, pB, alpha, beta, cutoff ! third potential
    integer(I1B), save, private, protected :: ty 
contains


function coord_r(i)
    integer(I4B), intent(in) :: i
    real   (RDP) :: coord_r
    coord_r = dr*dble(i)
end function coord_r


function coord_theta(i)
    integer(I4B), intent(in) :: i
    real   (RDP) :: coord_theta
    coord_theta = dtheta*dble(i)
end function coord_theta


function nabla_x(int_i, int_j)
    use math_const, only: pi => math_pi 
    integer(I4B), intent(in) :: int_i, int_j
    real   (RDP) :: nabla_x, i, j

    i = dble(int_i)
    j = dble(int_j)
    if(int_i /= int_j) then 
        nabla_x = &
            2.d0/(i -j)**2.d0
    else
        nabla_x = &
            pi**2.d0 /3.d0
    end if
    nabla_x = nabla_x/dr**2.d0
    if(mod(int_i +int_j, 2) == 0) then 
        nabla_x = -nabla_x
    end if 
!        *(-1.d0)**(i -j +1.d0)
end function nabla_x


function Poten(r) 
    real(RDP), intent(in) :: r
    real(RDP) :: Poten, stat, pol 

    if(ty == 1) then 
        stat  = -exp(-2.d0*z0*r)*(z1 +z2/r)
        pol   = -alphab/(2.d0*r**4.d0)*(1.d0 -exp(-r))**6.d0 
        Poten = stat +pol

    else if(ty == 2) then 
        stat  = -(Z/pH)*(exp(-r/pD)/r)*(1.d0 +(pH -1.d0)*exp(-r*pH/pD))
        pol   = -alphab/(2.d0*(r**2.d0 +pdelta**2.d0)**2.d0)
        Poten = stat +pol

    else if(ty == 3) then 
        stat  = -Z*( &
                        +pA**2.d0/(4.d0*alpha**3.d0)*exp(-2.d0*alpha*r)*(alpha +1.d0/r) &
                        +4.d0*pA*pB/(alpha +beta)**3.d0*exp(-(alpha +beta)*r)*((alpha +beta)/2.d0 +1.d0/r) & 
                        +pB**2.d0/(4.d0*beta**3.d0)*exp(-2.d0*beta*r)*(beta +1.d0/r))
        pol   = -alphab/(2.d0*(cutoff**2.d0 +r**2.d0)**2.d0)
        Poten = stat +pol

    end if 
end function Poten










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_input
    use math_const, only: pi => math_pi
    use unit_const, only: other_e_eV, au_hartree
!     integer (I1B), parameter :: file_input       = 101
    character(30), parameter :: form_particle    = '(////, 2(40X, 1F10.4, /), /)'
    character(60), parameter :: form_potential   = '(////, 1(40X, 1I10, /), 12(40X, 1F10.4, /), /)'
    character(60), parameter :: form_calculation = '(////, 1(40X, 1F10.4, /), 2(40X, 1I10, /), /)'
    character(30), parameter :: form_plot        = '(////, 2(40X, 1I10, /), /)'
    real    (RDP), parameter :: eV_to_au         = other_e_ev/au_hartree

!     open(file_input, file = "inout/input.d")
    read(*, form_particle)    Mass, Kinet
    read(*, form_potential)   ty, Z, alphab, z0, z1, z2, pH, pdelta, pA, pB, alpha, beta, cutoff
    read(*, form_calculation) ra, N, L 
    read(*, form_plot)        pr, ptheta 
!     close(file_input) 

    if(pr > N) pr = N 
    Kinet  = Kinet*eV_to_au
    dr     = ra/dble(N) 
    dtheta = pi/dble(ptheta)
    pD     = pH/Z**0.4d0 

    allocate(H(1:N, 1:N), E(1:N))
    allocate(R(0:L), K(0:L), S(0:L), A(0:L))
    allocate(inner_u(0:L, 1:N)) 
    H(:, :)       = 0.d0
    E(:)          = 0.d0
    R(:)          = 0.d0
    K(:)          = 0.d0
    S(:)          = 0.d0
    inner_u(:, :) = 0.d0
end subroutine PROC_input


subroutine PROC_Poten_plot
    integer (I1B), parameter :: file_poten = 101
    character(30), parameter :: form_poten = '(30ES25.10)'
    integer (I4B) :: i 

    open(file_poten, file = "inout/poten.d")
    do i = 1, N, N/pr 
        write(file_poten, form_poten) coord_r(i), Poten(coord_r(i))
    end do 
    close(file_poten)
    write(*, *) "Potential Type ", ty 
end subroutine PROC_Poten_plot


subroutine PROC_out 
    deallocate(H, E)
    deallocate(R, K, S, A)
    deallocate(inner_u)
end subroutine PROC_out 
end module hamiltonian 
