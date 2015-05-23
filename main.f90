program main
    use kind_type
    use global
    use hamiltonian, only: PROC_input, PROC_inform, PROC_Poten_plot, PROC_out 
    use basis,       only: PROC_H, PROC_basis_plot
    use boundary,    only: PROC_boundary_mat
    use inner,       only: PROC_inner_achive, PROC_inner_plot
    use outer,       only: PROC_CS_plot, PROC_outer_plot
    implicit none
    integer(I4) :: i
    real   (SP) :: t1, t2 

    write(*, *) "Reading input file..."
    call cpu_time(t1)
    call PROC_input  
    write(*, *) "Reading over."

    write(file_log, *) "INPUT"
        call PROC_inform
        call PROC_Poten_plot 
    call cpu_time(t2)
    write(file_log, *) "INPUT RUNNING TIME", t2 -t1 
    write(file_log, *) 
    write(file_log, *) 

!     i = L  
    write(*, *) "Calculating..."
    do i = 0, L 
        write(file_log, *) "ROUND", i 
        call cpu_time(t1)
            call PROC_H(i) 
            if(i == L) then 
                write(file_log, *) "BASIS FUNCTION PLOTED"
                write(file_log, *) "~~~~~~~~~~~~~~~~~~~~~"
                call PROC_basis_plot
            end if 
            call PROC_boundary_mat(i) 
            call PROC_inner_achive(i)
        call cpu_time(t2) 
        write(file_log, *) "ROUND RUNNING TIME", t2 -t1
        write(file_log, *) 
    end do 
    write(file_log, *) 
    write(*, *) "Calculating over."

    write(*, *) "Ploting..."
    write(file_log, *) " PLOT"
    call cpu_time(t1)
        call PROC_CS_plot 
        call PROC_inner_plot
        call PROC_outer_plot 
    call cpu_time(t2)
    write(file_log, *) " PLOT RUNNING TIME", t2 -t1 
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Calculating over."

    write(*, *) "Program over."
    write(file_log, *) "PROGRAM OVER"
    write(file_log, *)
    call PROC_out 
    close(file_log)
end program main
