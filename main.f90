program main
    use kind_type
    use global
    use hamiltonian, only: PROC_input, PROC_Poten_plot, PROC_out 
    use basis,       only: PROC_H, PROC_basis_plot
    use boundary,    only: PROC_boundary_mat
    use inner,       only: PROC_inner_achive, PROC_inner_plot
    use outer,       only: PROC_CS_plot, PROC_outer_plot
    implicit none
    integer(I4B) :: i
    real   (RSP) :: t1, t2 

    write(*, *) "INPUT"
    call cpu_time(t1)
        call PROC_input  
!         call PROC_Poten_plot 
    call cpu_time(t2)
    write(*, *) "INPUT RUNNING TIME", t2 -t1 
    write(*, *) 

!     i = L  
    do i = 0, L 
        write(*, *) "ROUND", i 
        call cpu_time(t1)
            call PROC_H(i) 
!             if(i == L) call PROC_basis_plot
            call PROC_boundary_mat(i) 
            call PROC_inner_achive(i)
        call cpu_time(t2) 
        write(*, *) "ROUND RUNNING TIME", t2 -t1
        write(*, *) 
    end do 

    write(*, *) " PLOT"
    call cpu_time(t1)
        call PROC_CS_plot 
!         call PROC_inner_plot
!         call PROC_outer_plot 
    call cpu_time(t2)
    write(*, *) " PLOT RUNNING TIME", t2 -t1 
    write(*, *)

    write(*, *) "PROGRAM OVER"
    call PROC_out 
end program main
