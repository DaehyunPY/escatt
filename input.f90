module input
    implicit none
    real(8), parameter :: Ra = 100.d0 
    integer, parameter :: CooN = 1000
    real(8), parameter :: dr = Ra/dble(CooN)

    real(8), parameter :: Rb = 20.d0 
    integer, parameter :: CooK = int(Rb/dr)

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


function nabla_limit(dx, kk, ll, ii, jj)
    real(8), intent(in) :: dx
    real(8), parameter  :: pi = 2.0d0*acos(0.0d0)
    integer, intent(in) :: ii, jj, kk, ll
    real(8)  :: nabla_limit, i, j, k, l, N
    real(16) :: tmp 
    k = dble(kk)
    l = dble(ll)
    i = dble(ii)
    j = dble(jj)
    N = 10.d0**5 
    if( ii /= jj ) then 
        tmp = &
!             (((4*pi**2*l**2-4*pi**2*l+pi**2)*cos((j+i)*(2*pi*l+5*pi)/k/2.d+0) &
!             +(-20*pi**2*l**2+12*pi**2*l+3*pi**2)*cos((j+i)*(2*pi*l+3*pi)/k/2.d+0) &
!             +(40*pi**2*l**2-8*pi**2*l-14*pi**2)*cos((j+i) &
!             *(2*pi*l+pi)/k/2.d+0)+(-40*pi**2*l**2-8*pi**2*l+14*pi**2)*cos((j+i) &
!             *(2*pi*l-pi)/k/2.d+0)+(20*pi**2*l**2+12*pi**2*l-3*pi**2) &
!             *cos((j+i)*(2*pi*l-3*pi)/k/2.d+0)+(-4*pi**2*l**2-4*pi**2*l-pi**2) &
!             *cos((j+i)*(2*pi*l-5*pi)/k/2.d+0))/(8*cos(3*pi*(j+i)/k) &
!             *k**2-48*cos(2*pi*(j+i)/k)*k**2+120*cos(pi*(j+i)/k)*k**2-80*k**2) &
!             -((4*pi**2*l**2-4*pi**2*l+pi**2)*cos((i-j)*(2*pi*l+5*pi)/k/2.d+0) &
!             +(-20*pi**2*l**2+12*pi**2*l+3*pi**2)*cos((i-j)*(2*pi*l+3*pi)/k/2.d+0) &
!             +(40*pi**2*l**2-8*pi**2*l-14*pi**2)*cos((i-j)*(2*pi*l+pi)/k/2.d+0) &
!             +(-40*pi**2*l**2-8*pi**2*l+14*pi**2)*cos((i-j)*(2*pi*l-pi)/k/2.d+0) &
!             +(20*pi**2*l**2+12*pi**2*l-3*pi**2)*cos((i-j)*(2*pi*l-3*pi)/k/2.d+0) &
!             +(-4*pi**2*l**2-4*pi**2*l-pi**2)*cos((i-j)*(2*pi*l-5*pi)/k/2.d+0)) &
!             /(8*cos(3*pi*(i-j)/k)*k**2-48*cos(2*pi*(i-j)/k)*k**2+120 &
!             *cos(pi*(i-j)/k)*k**2-80*k**2))/(dx**2*l)
            (((4*pi**2*N**2-4*pi**2*N+pi**2)*cos((2*pi*(j+i)*N+5*pi*(j+i))/k/2.d+0) &
            +(-20*pi**2*N**2+12*pi**2*N+3*pi**2)*cos((2*pi*(j+i)*N &
            +3*pi*(j+i))/k/2.d+0)+(40*pi**2*N**2-8*pi**2*N-14*pi**2) &
            *cos((2*pi*(j+i)*N+pi*(j+i))/k/2.d+0)+(-40*pi**2*N**2-8*pi**2*N+14*pi**2) &
            *cos((2*pi*(j+i)*N-pi*(j+i))/k/2.d+0)+(20*pi**2*N**2 &
            +12*pi**2*N-3*pi**2)*cos((2*pi*(j+i)*N-3*pi*(j+i))/k/2.d+0) &
            +(-4*pi**2*N**2-4*pi**2*N-pi**2)*cos((2*pi*(j+i)*N-5*pi*(j+i))/k/2.d+0)) &
            /(8*cos(3*pi*(j+i)/k)*k**2-48*cos(2*pi*(j+i)/k)*k**2 &
            +120*cos(pi*(j+i)/k)*k**2-80*k**2)-((4*pi**2*N**2-4*pi**2*N+pi**2) &
            *cos((2*pi*(i-j)*N+5*pi*(i-j))/k/2.d+0) &
            +(-20*pi**2*N**2+12*pi**2*N+3*pi**2)*cos((2*pi*(i-j)*N+3*pi*(i-j))/k/2.d+0) &
            +(40*pi**2*N**2-8*pi**2*N-14*pi**2)*cos((2*pi*(i-j)*N+pi*(i-j))/k/2.d+0) &
            +(-40*pi**2*N**2-8*pi**2*N+14*pi**2)*cos((2*pi*(i-j)*N-pi*(i-j))/k/2.d+0) &
            +(20*pi**2*N**2+12*pi**2*N-3*pi**2)*cos((2*pi*(i-j)*N-3*pi*(i-j))/k/2.d+0)+( &
            -4*pi**2*N**2-4*pi**2*N-pi**2)*cos((2*pi*(i-j)*N-5*pi*(i-j))/k/2.d+0)) &
            /(8*cos(3*pi*(i-j)/k)*k**2-48*cos(2*pi*(i-j)/k)*k**2+120*cos(pi*(i-j)/k)*k**2-80*k**2))/(dx**2*l)
    else
        tmp = &
!             (((4*pi**2*l**2-4*pi**2*l+pi**2)*cos((j+i)*(2*pi*l+5*pi)/k/2.d+0) &
!             +(-20*pi**2*l**2+12*pi**2*l+3*pi**2)*cos((j+i)*(2*pi*l+3*pi)/k/2.d+0) &
!             +(40*pi**2*l**2-8*pi**2*l-14*pi**2)*cos((j+i)*(2*pi*l+pi)/k/2.d+0) &
!             +(-40*pi**2*l**2-8*pi**2*l+14*pi**2)*cos((j+i)*(2*pi*l-pi)/k/2.d+0) &
!             +(20*pi**2*l**2+12*pi**2*l-3*pi**2)*cos((j+i)*(2*pi*l-3*pi)/k/2.d+0) &
!             +(-4*pi**2*l**2-4*pi**2*l-pi**2)*cos((j+i)*(2*pi*l-5*pi)/k/2.d+0)) &
!             /(8*cos(3*pi*(j+i)/k)*k**2-48*cos(2*pi*(j+i)/k)*k**2 &
!             +120*cos(pi*(j+i)/k)*k**2-80*k**2)+pi**2*(-(2*l**3+3*l**2+l)/6.d+0 &
!             +(l**2+l)/2.d+0-l/4.d+0)/k**2)/(dx**2*l)
            (((4*pi**2*N**2-4*pi**2*N+pi**2)*cos((2*pi*(j+i)*N+5*pi*(j+i))/k/2.d+0) &
            +(-20*pi**2*N**2+12*pi**2*N+3*pi**2)*cos((2*pi*(j+i)*N+3*pi*(j+i))/k/2.d+0) &
            +(40*pi**2*N**2-8*pi**2*N-14*pi**2) &
            *cos((2*pi*(j+i)*N+pi*(j+i))/k/2.d+0)+(-40*pi**2*N**2-8*pi**2*N+14*pi**2) &
            *cos((2*pi*(j+i)*N-pi*(j+i))/k/2.d+0)+(20*pi**2*N**2+12*pi**2*N-3*pi**2) &
            *cos((2*pi*(j+i)*N-3*pi*(j+i))/k/2.d+0) &
            +(-4*pi**2*N**2-4*pi**2*N-pi**2)*cos((2*pi*(j+i)*N-5*pi*(j+i))/k/2.d+0)) &
            /(8*cos(3*pi*(j+i)/k)*k**2-48*cos(2*pi*(j+i)/k)*k**2 &
            +120*cos(pi*(j+i)/k)*k**2-80*k**2)+pi**2*(-(2*N**3+3*N**2+N)/6.d+0+(N**2+N)/2.d+0-N/4.d+0)/k**2)/(dx**2*l)
    end if
    nabla_limit = tmp 
end function nabla_limit


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
