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

  ! laplacian matrix [K] from the interior part of A_psi 
  real(dp), allocatable   :: K_omega(:,:)
  ! L-S RK3 Coefficients from Spalart, Moser & Rogers (1991)
  ! real(dp), dimension(3)  :: rk_alpha = (/ 29.0d0 / 60.0d0, -3.0d0 / 40.0d0, 1.0d0 / 6.0d0 /)
  ! real(dp), dimension(3)  :: rk_beta  = (/ 37.0d0 / 160.0d0, 5.0d0 / 24.0d0, 1.0d0 / 6.0d0 /)
  ! real(dp), dimension(3)  :: rk_gamma = (/ 8.0d0  / 15.0d0, 5.0d0 / 12.0d0, 3.0d0 / 4.0d0 /)
  ! real(dp), dimension(3)  :: rk_zeta  = (/ 0.0d0, -17.0d0 / 60.0d0, -5.0d0 / 12.0d0 /)
  ! intermediate storage vectors for RK3
  real(dp), allocatable   :: d_omega_s1(:), d_omega_s2(:) ! step 1 and 2 results 
  real(dp), allocatable   :: N_omega_n(:), N_omega_s1(:), N_omega_s2(:) ! advection terms 
  real(dp), allocatable   :: A_step(:,:), b_step(:)
  


  write(*,*) "======================================="
  write(*,*) "L-shaped Cavity flow solver"
  write(*,*) "======================================="

  write(*,*) "1. Initializing grid and solver..."
  allocate(blocks(NUM_BLOCKS))
  do iblock = 1, NUM_BLOCKS 
    call initialize_block(blocks(iblock), iblock)
    call generate_grid(blocks(iblock))
  enddo

  ! initialize solvers global numbering system 
  call initialize_solver(blocks)

  N                     = (NX_B1 * NY_B1) + (NX_B2 * NY_B2)
  total_unknowns        = 2 * N


  ! TODO: we can write a memory allocator module
  allocate(A_psi(N, N), M_psi(N, N), M_omega(N, N))
  allocate(d_psi_old(N), d_omega_old(N), d_psi_new(N), d_omega_new(N))
  allocate(b_psi_total(N), b_psi_bc(N), b_omega_rhs(N))
  allocate(K_omega(N, N))
  allocate(d_omega_s1(N), d_omega_s2(N))
  allocate(N_omega_n(N), N_omega_s1(N), N_omega_s2(N))
  allocate(A_step(N, N), b_step(N))
  allocate(piv(N)) ! lapack 
  allocate(A_psi_copy(N,N)) ! copy for psi during RK3


  ! initial conditions 
  d_psi_old = 0.0d0; d_omega_old = 0.0d0; d_psi_new = 0.0d0; d_omega_new = 0.0d0

  write(*,*) "2. Assembling time-independent matrices..."
  call assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
  call extract_laplacian_operator(blocks, A_psi, K_omega)

  write(*,*) "3. Starting time-stepping loop..."
  n_steps = ceiling(t_end / dt)

! The main time loop
  do t = 1, n_steps

    ! ====================================================================
    ! STEP 0: UPDATE STATE FOR THE CURRENT STEP
    ! ====================================================================
    d_psi_old = d_psi_new
    d_omega_old = d_omega_new
    current_time = (t - 1) * dt

    ! ====================================================================
    ! L-S RK3 STAGE 1
    ! ====================================================================
    ! Calculate N(omega_n) using the state at the beginning of the step
    call calculate_advection_term(blocks, d_psi_old, d_omega_old, N_omega_n)
    write(*,*) 'advection term calculated for step 1'

    ! Assemble the system for step 1
    call assemble_rk3_step(blocks, 1, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)

    ! solve the system A_step * x = b_step 
    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    d_omega_s1 = b_step ! dgesv puts the solution into b_step
    write(*,*) 'step1 complete'

    ! ====================================================================
    ! L-S RK3 STAGE 2
    ! ====================================================================
    ! Calculate N(omega_s1)
    call calculate_advection_term(blocks, d_psi_old, d_omega_s1, N_omega_s1)
    write(*,*) 'advection term calculated for step 2'

    call assemble_rk3_step(blocks, 2, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)

    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    d_omega_s2 = b_step ! dgesv puts the solution into b_step
    write(*,*) 'step2 complete'

    ! ====================================================================
    ! L-S RK3 STAGE 3
    ! ====================================================================
    ! Calculate N(omega_s2)
    call calculate_advection_term(blocks, d_psi_old, d_omega_s2, N_omega_s2)
    write(*,*) 'advection term calculated for step 3'

    call assemble_rk3_step(blocks, 3, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)

    ! for d_omega_new (the final vorticity at the end of the time step)
    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    d_omega_new = b_step
    write(*,*) 'step3 complete'

    ! ====================================================================
    ! FINAL PSI UPDATE
    ! ====================================================================
    ! Apply boundary conditions at the END of the time step (t_n+1)
    call apply_psi_bcs(blocks, A_psi, b_psi_bc, t * dt) ! 
    b_psi_total = matmul(-M_omega, d_omega_new) + b_psi_bc

    A_psi_copy = A_psi 

    call DGESV(N, 1, A_psi_copy, N, piv, b_psi_total, N, linfo)
    d_psi_new = b_psi_total

    ! output
    if (mod(t, output_interval) == 0) then
      write(*,'(A, I6, A, I6, A, F8.4)') "Timestep:", t, "/", n_steps, "  Time:", t * dt
    endif
  enddo

  write(*,*) "Deallocating memory..."
  ! clean up TODO: write a module for this maybe 
  deallocate(A_psi, M_omega, d_psi_old, d_omega_old, d_psi_new, d_omega_new)
  deallocate(b_psi_total, b_psi_bc, b_omega_rhs)
  deallocate(K_omega, A_step, b_step, d_omega_s1, d_omega_s2, N_omega_n, N_omega_s1, N_omega_s2)
  deallocate(piv, A_psi_copy)
  do iblock = 1, NUM_BLOCKS 
    deallocate(blocks(iblock)%colloc_pts, blocks(iblock)%boundary_types)
    deallocate(blocks(iblock)%knots_x, blocks(iblock)%knots_y)
  enddo 
  deallocate(blocks)
  write(*,*) "Finished!!!!"

        


end program L_shaped_cavity_solver
