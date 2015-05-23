program main
    use global 
    use basis
    use boundary
    use inner
    use outer
    implicit none
    character(30), parameter :: form_out = '(1A5, 1I5)'
    integer :: i 

    call PROC_input  
    call PROC_Poten_plot 
!     i = L  
    do i = 0, L 
        write(*, form_out) "ROUND", i 
        call PROC_H(i) 
!         call PROC_basis_plot
        call PROC_boundary_mat(i) 
        call PROC_inner_achive(i)
        write(*, *) 
    end do 
    call PROC_CS_plot 
    call PROC_inner_plot
    call PROC_outer_plot 
    call PROC_out 

end program main
