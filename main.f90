program main
    use basis_process
    use boundary_process
    implicit none
    real(4) :: t1, t2 

    call PROC_H
    call PROC_Psi 
    call PROC_plot 
    call PROC_R
    call PROC_K 
    call PROC_S 
end program main
