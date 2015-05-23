module local 
    implicit none
    real(8), save :: Z, alpha1, alpha2, alpha3, alphab, cutoff
    integer, save :: tp 
end module local 


module global 
    implicit none
    real(8),    save, protected   :: Mass, Kinet, ra
    real(8),    save, protected   :: dr, dtheta 
    integer,    save, protected   :: N, L, pr, ptheta
    real(8),    save, allocatable :: H(:, :), E(:), R(:), K(:)
    complex(8), save, allocatable :: S(:), inner_u(:, :)
contains


function Poten(r) 
    use local 
    real(8), intent(in) :: r
    real(8) :: Poten, stat, pol 

    stat  = -(Z/r)*exp(-alpha1*r)-alpha2*exp(-alpha3*r) 
    pol   = 0.d0 
    if(tp == 1) then 
    else if(tp == 2) then 
        pol   = -(alphab/(2.d0*r**4.d0))*(1.d0 -exp((r/cutoff)**3.d0))**2.d0 
    else 
        stop "FUNCTION Poten: Something is wrong."
    end if 
    Poten = stat +pol
end function Poten


function coord_r(i)
    use local 
    integer, intent(in) :: i
    real(8) :: coord_r

    coord_r = dr*dble(i)
end function coord_r


function coord_theta(i)
    use local 
    integer, intent(in) :: i
    real(8) :: coord_theta

    coord_theta = dtheta*dble(i)
end function coord_theta


FUNCTION plgndr_s(l,m,x)
    USE nrtype; USE nrutil, ONLY : arth,assert
    INTEGER(I4B), INTENT(IN) :: l,m
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: plgndr_s
    INTEGER(I4B) :: ll
    REAL(DP) :: pll,pmm,pmmp1,somx2
    call assert(m >= 0, m <= l, abs(x) <= 1.0, 'plgndr_s args')
    pmm=1.0
    if (m > 0) then
        somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
        pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
        if (mod(m,2) == 1) pmm=-pmm
    end if
    if (l == m) then
        plgndr_s=pmm
    else
        pmmp1=x*(2*m+1)*pmm
        if (l == m+1) then
            plgndr_s=pmmp1
        else
            do ll=m+2,l
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                pmm=pmmp1
                pmmp1=pll
            end do
            plgndr_s=pll
        end if
    end if
END FUNCTION plgndr_s










! ==================================================
! PROCESS
! ==================================================


subroutine PROC_input
    use const 
    use local 
    integer,       parameter :: file_input       = 101
    character(30), parameter :: form_particle    = '(////, 2(40X, 1F10.4, /), /)'
    character(60), parameter :: form_potential   = '(////, 1(40X, 1I10, /), 6(40X, 1F10.4, /), /)'
    character(60), parameter :: form_calculation = '(////, 1(40X, 1F10.4, /), 2(40X, 1I10, /), /)'
    character(30), parameter :: form_plot        = '(////, 2(40X, 1I10, /), /)'
    real(8),       parameter :: pi       = math_pi
    real(8),       parameter :: eV_to_au = other_e_eV/au_hartree

    open(file_input, file = "inout/input.d")
    read(file_input, form_particle)    Mass, Kinet
    read(file_input, form_potential)   tp, Z, alpha1, alpha2, alpha3, alphab, cutoff
    read(file_input, form_calculation) ra, N, L 
    read(file_input, form_plot)        pr, ptheta 
    close(file_input) 

    Kinet  = Kinet*eV_to_au
    dr     = ra/dble(N) 
    dtheta = pi/dble(ptheta)

    allocate(H(1:N, 1:N), E(1:N))
    allocate(R(0:L), K(0:L), S(0:L))
    allocate(inner_u(0:L, 1:N))
    H(:, :) = 0.d0 
    E(:) = 0.d0 
    R(:) = 0.d0 
    K(:) = 0.d0 
    S(:) = 0.d0 
    inner_u(:, :) = 0.d0 
end subroutine PROC_input


subroutine PROC_Poten_plot 
    integer,       parameter :: file_poten = 101
    character(30), parameter :: form_poten = '(30ES25.10)'
    integer :: i 

    open(file_poten, file = "inout/poten.d")
    do i = 1, N, N/pr 
        write(file_poten, form_poten) coord_r(i), Poten(coord_r(i))
    end do 
    close(file_poten)
end subroutine PROC_Poten_plot


subroutine PROC_out 
    deallocate(H, E)
    deallocate(R, K, S)
    deallocate(inner_u)
end subroutine PROC_out 
end module global 
