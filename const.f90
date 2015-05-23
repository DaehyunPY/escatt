module kind_type
    use iso_fortran_env 
    implicit none
    integer, parameter :: &
        I1B = INT8, &  ! -128                       ~ 127
        I2B = INT16, & ! -32,768                    ~ 32,767
        I4B = INT32, & ! -2,147,483,648             ~ 2,147,483,647
        I8B = INT64, & ! -9,223,372,036,854,775,808 ~ 9,223,372,036,854,775,807
        RSP = REAL32, &
        RDP = REAL64, &
        RQP = REAL128, &
        CSP = REAL32, &
        CDP = REAL64, &
        CQP = REAL128 
end module ! kind_type


module math_const
    use kind_type, only: RDP, CDP 
    implicit none 
! real 
    real(RDP), parameter :: &
        math_pi     = 2.d0*acos(0.d0), &
        math_e      = exp(1.d0), &
        math_degree = math_pi/180.d0
! complex 
    complex(CDP), parameter :: &
        math_i = (0.d0, 1.d0)
end module ! math_const 


module unit_const
    use kind_type,  only: RDP 
    use math_const, only: math_pi
    implicit none 
!si
    real(RDP), parameter :: & 
        si_c       = 2.99792458d8, &                 !  (m s-)si
        si_Na      = 6.02214129d+23, &               !  (mol-)si
        si_kb      = 1.38065042d-23, &               !  (J K-)si
        si_epsilon = 1.d+7/(4*math_pi*si_c**2.d0), & !  (F m-)si
        si_mu      = 1.d-7*(4*math_pi), &            ! (N A-2)si
        si_mass    = 9.10938262d-31, &               !    (kg)si
        si_charge  = 1.60217653d-19, &               !     (C)si
        si_hbar    = 1.05457173d-34                  !   (J s)si
!other
    real(RDP), parameter :: &
        other_w_debye = 1.d-21/si_c, &                       ! (C m)si other-(debye)
        other_w_pcm   = (2.d0*math_pi*si_c)/1.d-2, &         !  (Hz)si other-(cm-)
        other_e_pcm   = (2.d0*math_pi*si_c)*si_hbar/1.d-2, & !   (J)si other-(cm-)
        other_e_eV    = si_charge                            !   (J)si other-(eV)
!au
    real(RDP), parameter :: &
        au_hartree     = si_mass*si_charge**4.d0/(4.d0*math_pi*si_epsilon*si_hbar)**2.d0, & ! (J)si au-(hartree)
        au_bohr        = 4.d0*math_pi*si_epsilon*si_hbar**2.d0/(si_mass*si_charge**2.d0), & ! (m)si au-(bohr)
        au_mass        = si_mass, &                        !   (kg)si au-
        au_time        = si_hbar/au_hartree, &             !    (s)si au-
        au_temperature = au_hartree/si_kb, &               !    (K)si au-
        au_E_charge    = si_charge, &                      !    (C)si au-
        au_E_field     = au_hartree/au_bohr/au_E_charge, & ! (V m-)si au-
        au_B_charge    = au_hartree/au_E_charge*au_time, & !   (Wb)si au-
        au_B_field     = au_hartree/au_bohr/au_B_charge    ! (A m-)si au-
end module ! unit_const
