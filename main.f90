! Main program for the 2D lid-driven cavity solver.
! Solves the 2D unsteady Navier-stokes equations in vorticity-streamfunction form 
! for an L-shaped cavity using a multi-block B-spline collocation method
! 09/04/2025
! Ali Yesildag


program L_shaped_cavity_solver

  use mod_params
  use mod_bspline
  use mod_grid
  use mod_solver
  implicit none

  ! TODO: we can write a variable definition module
  
  ! Simulation parameters 
  real(dp)    :: t_end = 0.01
  real(dp)    :: dt = 0.001d0 
  integer     :: n_steps 
  integer     :: output_interval = 100
  real(dp)    :: current_time 

  ! multiblock grid data structure 
  type(BLOCK_TYPE), allocatable :: blocks(:)

  ! global solver arrays 
  integer :: N , total_unknowns
  ! time independent matrices 
  real(dp), allocatable :: A_psi(:,:), M_psi(:,:), M_omega(:,:)
  ! time dependent coefficient vectors 
  real(dp), allocatable :: d_psi_old(:), d_omega_old(:), d_psi_new(:), d_omega_new(:)
  ! RHS vectors 
  real(dp), allocatable :: b_psi_total(:), b_psi_bc(:), b_omega_rhs(:)

  ! loop counters 
  integer :: t, iblock, k 


  write(*,*) "======================================="
  write(*,*) "L-shaped Cavity flow solver"
  write(*,*) "======================================="

  write(*,*) "1. Initializing grid and solver..."
  allocate(blocks(NUM_BLOCKS))
  do iblock = 1, NUM_BLOCKS 
    call initialize_block(blocks(iblock), iblock)
  enddo
    call generate_global_grid(blocks(iblock))

  ! initialize solvers global numbering system 
  call initialize_solver(blocks)

  N                     = (NX_B1 * NY_B1) + (NX_B2 * NY_B2)
  total_unknowns        = 2 * N


  ! TODO: we can write a memory allocator module
  allocate(A_psi(N, N), M_psi(N, N), M_omega(N, N))
  allocate(d_psi_old(N), d_omega_old(N), d_psi_new(N), d_omega_new(N))
  allocate(b_psi_total(N), b_psi_bc(N), b_omega_rhs(N))

  ! initial conditions 
  d_psi_old = 0.0d0; d_omega_old = 0.0d0; d_psi_new = 0.0d0; d_omega_new = 0.0d0

  write(*,*) "2. Assembling time-independent matrices..."
  call assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)

  write(*,*) "3. Starting time-stepping loop..."
  n_steps = ceiling(t_end / dt)

  do t = 1, n_steps 
    current_time = t * dt 

    ! Apply time-dependent bcs to psi 
    call apply_psi_bcs(blocks, A_psi, b_psi_bc, current_time)
    ! Calculate the RHS of the vorticity equation 
    ! RHS = -(u.grad)omega + (1/Re) * nabla^2(omega)
    ! call calculate_vorticity_rhs(blocks, d_psi_old, d_omega_old, b_omega_rhs)

    ! update interior vorticity with forward stepping 
    ! omega_new = omega_old + dt * RHS 
    d_omega_new = d_omega_old + dt * b_omega_rhs 

    ! enforce wall vorticity bcs at new time 
    ! call apply_omega_bcs(blocks, d_psi_old, d_omega_new)

    ! solve for the new streamfunction 
    ! A_psi * d_psi_new = -M_omega * d_omega_new + b_psi_bc 
    ! b_psi_total = matmul(-M_omega, d_omega_new) + b_psi_bc 
    !
    ! CALL LAPACK SOLVER 
    ! call DGESV(N, 1, A_psi, N, ipiv, b_psi_total, N, info )
    ! d_psi_new = b_psi_total 


    ! update for next iteration 
    ! d_psi_old = d_psi_new 
    ! d_omega_old = d_omega_new 
    !
    ! output 
    if (mod(t, output_interval) == 0) then 

    endif 
  enddo

  write(*,*) "Deallocating memory..."
  ! clean up TODO: write a module for this maybe 
  deallocate(A_psi, M_omega, d_psi_old, d_omega_old, d_psi_new, d_omega_new)
  deallocate(b_psi_total, b_psi_bc, b_omega_rhs)
  do iblock = 1, NUM_BLOCKS 
    deallocate(blocks(iblock)%colloc_pts, blocks(iblock)%boundary_types)
    deallocate(blocks(iblock)%knots_x, blocks(iblock)%knots_y)
  enddo 
  deallocate(blocks)
  write(*,*) "Finished!!!!"

        


end program L_shaped_cavity_solver
