module input
    implicit none
!     real(8), parameter :: Ra = 500.d0 
    real(8), parameter :: Ra = 20.d0 
    integer, parameter :: CooN = 5000
    real(8), parameter :: dr = Ra/dble(CooN)
    real(8) H(1:CooN, 1:CooN), E(1:CooN)
    integer,       parameter :: file_psi = 101,     file_energy = 102,    file_ana = 103
    character(30), parameter :: form_psi = '(6ES)', form_energy = '(ES)', form_ana = '(6ES)'
contains


function Poten(r)
    real(8), intent(in) :: r
    real(8) Poten
    Poten = 2.d0/r
end function Poten
function CooR(i)
    integer, intent(in) :: i
    real(8) CooR
    CooR = dr*dble(i)
end function CooR


function nabla_r_close_open(dx, ii, jj)
    real(8), intent( in) :: dx
    integer, intent( in) :: ii, jj
    real(8) :: nabla_r_close_open, i, j
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)
    i = dble(ii)
    j = dble(jj)
    if( ii /= jj ) then 
        nabla_r_close_open = &
            (pi**2*((i-j)**3*(j+i)**2*sin(pi*(j+i))-sin(pi*(i-j))*(i-j)**2*(j+i)**3) &
            -2*(i-j)**3*sin(pi*(j+i))+pi*(2*(i-j)**3*(j+i)*cos(pi*(j+i)) &
            -2*cos(pi*(i-j))*(i-j)*(j+i)**3)+2*sin(pi*(i-j))*(j+i)**3) &
            /(pi*dx**2*(i-j)**3*(j+i)**3)
    else
        nabla_r_close_open = &
            -(-3*pi**2*(j+i)**2*sin(pi*(j+i))+6*sin(pi*(j+i))-6*pi*(j+i)*cos(pi*(j+i)) &
            +pi**3*(j+i)**3)/(pi*dx**2*(j+i)**3)/3.d+0
    end if
end function nabla_r_close_open
function nabla_r_close_close(dx, ii, jj)
    real(8), intent(in) :: dx
    integer, intent(in) :: ii, jj
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)
    real(8) nabla_r_close_close, i, j
    i = dble(ii)
    j = dble(jj)
    if( ii /= jj ) then 
        nabla_r_close_close = &
!           -(-1.d0)**(i-j)/dx**2.d0 *(2.d0/(i-j)**2.d0 -2.d0/(i+j)**2.d0)
            (pi**2*((i-j)**3*(j+i)**2*sin(pi*(j+i))-sin(pi*(i-j))*(i-j)**2*(j+i)**3) &
            -2*(i-j)**3*sin(pi*(j+i))+2*pi*(i-j)**3*(j+i)*cos(pi*(j+i)) &
            -2*pi*cos(pi*(i-j))*(i-j)*(j+i)**3+2*sin(pi*(i-j))*(j+i)**3) &
            /(pi*dx**2*(i-j)**3*(j+i)**3)
    else
        nabla_r_close_close = &
!           -(-1.d0)**(i-j)/dx**2.d0 *(pi**2.d0/3.d0 -1/2.d0/i**2.d0) 
            (3*pi**2*(j+i)**2*sin(pi*(j+i))-6*sin(pi*(j+i))+6*pi*(j+i)*cos(pi*(j+i)) &
            -pi**3*(j+i)**3)/(pi*dx**2*(j+i)**3)/3.d+0
    end if
end function nabla_r_close_close


function f_ana(ty, r)
    integer, intent(in) :: ty
    real(8), intent(in) :: r
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)
    real(8), parameter :: a = 1.d0 
    real(8) f_ana

    if(ty == 1) then 
        f_ana = 2.d0/a**1.5d0 *exp(-r/a)

    else if(ty == 2) then 
        f_ana = 1.d0/2.d0**1.5d0 *(2.d0 -r/a) *exp(-r/2.d0/a)

    else if(ty == 3) then 
        f_ana = 2.d0/3.d0**4.5d0 *(27.d0 -18.d0*r/a +2.d0*(r/a)**2.d0) *exp(-r/3.d0/a)

    end if 
    f_ana = f_ana*r
end function f_ana 
end module input
