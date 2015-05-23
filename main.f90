program main
!     use PROC_global 
    use PROC_basis
    use PROC_boundary
    use PROC_inner
    use PROC_outer
    implicit none
    integer :: i 

!     call PROC_input  
    
    i = l 
!     do i = 0, L 
        call PROC_H(i) 
!         call PROC_basis_plot
!         call PROC_boundary_mat(i) 
!         call PROC_inner_plot(i)
!         call PROC_CS(i) 
!         call PROC_outer_plot(i) 
!     end do 
end program main
