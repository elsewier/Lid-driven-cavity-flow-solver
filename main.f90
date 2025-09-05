! Main program for the 2D lid-driven caivty solver.
! 09/04/2025
! Ali Yesildag


program L_shaped_cavity_solver

        use mod_params
        use mod_bspline
        use mod_solver
        implicit none

        ! TODO: we can write a variable definition module
        integer :: num_colloc_pts, N
        real(dp), allocatable :: knots_x(:), knots_y(:)
        real(dp), allocatable :: colloc_pts(:,:)
        integer, allocatable  :: boundary_types(:)
        real(dp), allocatable :: A(:,:), b(:)
        real(dp), allocatable :: coeffs(:)

        write(*,*) "======================================="
        write(*,*) "L-shaped Cavity flow solver"
        write(*,*) "======================================="

        num_colloc_pts  = N_x * N_y 
        N               = 2 * N_x * N_y ! dimension of matrix 

        ! memory allocation 
        ! TODO: we can write a memory allocator module
        allocate(knots_x(N_x + P_x + 1), knots_y(N_y + P_y + 1))
        allocate(colloc_pts(num_colloc_pts,2))
        allocate(boundary_types(num_colloc_pts))
        allocate(A(N, N), b(N))
        allocate(coeffs(N))

        ! generate knot vectors
        write(*,*) "Generating knot vectors"
        call generate_knot_vector(P_x, N_x, 2.5d0, knots_x) ! TODO: remove the stretch magic number
        call generate_knot_vector(P_y, N_y, 2.5d0, knots_y) ! TODO: remove the stretch magic number

        ! define collocation grid
        write(*,*) "Setting up collocation grid"
        call setup_collocation_grid(colloc_pts, boundary_types)

        ! assemble linear system
        call assemble_system(A, b, colloc_pts, boundary_types, knots_x, knots_y)

        ! now we can solve the system
        ! TODO:  
        

        ! clean up TODO: write a module for this maybe 
        deallocate(knots_x, knots_y, colloc_pts, boundary_types, A, b, coeffs)
        write(*,*) "Finished!!!!"

        


end program L_shaped_cavity_solver
