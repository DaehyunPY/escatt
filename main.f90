program main
    use kind_type
    use global
    use hamiltonian, only: coord_E, PROC_input, PROC_out 
    use basis,       only: PROC_H
    use boundary,    only: PROC_boundary_mat
    use outer,       only: PROC_CS_achive, PROC_E_CS_plot
    implicit none
    integer(I4B) :: i, j 
    real   (RSP) :: t1, t2 

    write(*, *) "INPUT"
    call cpu_time(t1)
        call PROC_input  
    call cpu_time(t2)
    write(*, *) "INPUT RUNNING TIME", t2 -t1 
    write(*, *) 

    do j = 1, pE 
        Kinet = coord_E(j)
        write(*, *) "ROUND", j, Kinet
        call cpu_time(t1)
!             do i = 0, L 
                i = 0
                call PROC_H(i) 
                call PROC_boundary_mat(i) 
!             end do ! i 
            call PROC_CS_achive(j)
        call cpu_time(t2) 
        write(*, *) "ROUND RUNNING TIME", t2 -t1
        write(*, *) 
    end do ! j 

    write(*, *) " PLOT"
    call cpu_time(t1)
        call PROC_E_CS_plot
    call cpu_time(t2)
    write(*, *) " PLOT RUNNING TIME", t2 -t1 
    write(*, *)

    write(*, *) "PROGRAM OVER"
    call PROC_out 
end program main
