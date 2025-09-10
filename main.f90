! =============================================================================
! FILE: main.f90 (Corrected Version)
! =============================================================================
program L_shaped_cavity_solver

  use mod_params
  use mod_bspline
  use mod_grid
  use mod_solver
  implicit none
  
  real(dp)    :: t_end = 0.1
  real(dp)    :: dt = 0.001d0 
  integer     :: n_steps, output_interval = 1
  real(dp)    :: current_time 
  type(BLOCK_TYPE), allocatable :: blocks(:)

  integer :: N , total_unknowns
  real(dp), allocatable :: A_psi(:,:), M_psi(:,:), M_omega(:,:)
  real(dp), allocatable :: d_psi_old(:), d_omega_old(:), d_psi_new(:), d_omega_new(:)
  real(dp), allocatable :: b_psi_total(:), b_psi_bc(:), b_omega_rhs(:)
  integer :: t, iblock
  real(dp), allocatable   :: K_omega(:,:)
  real(dp), allocatable   :: d_omega_s1(:), d_omega_s2(:)
  real(dp), allocatable   :: N_omega_n(:), N_omega_s1(:), N_omega_s2(:)
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

  call initialize_solver(blocks)
  N = (NX_B1 * NY_B1) + (NX_B2 * NY_B2)
  total_unknowns = 2 * N

  allocate(A_psi(N, N), M_psi(N, N), M_omega(N, N))
  allocate(d_psi_old(N), d_omega_old(N), d_psi_new(N), d_omega_new(N))
  allocate(b_psi_total(N), b_psi_bc(N), b_omega_rhs(N))
  allocate(K_omega(N, N))
  allocate(d_omega_s1(N), d_omega_s2(N))
  allocate(N_omega_n(N), N_omega_s1(N), N_omega_s2(N))
  allocate(A_step(N, N), b_step(N))
  allocate(piv(N))
  allocate(A_psi_copy(N,N))

  d_psi_old = 0.0d0; d_omega_old = 0.0d0; d_psi_new = 0.0d0; d_omega_new = 0.0d0

  write(*,*) "2. Assembling time-independent matrices..."
  call assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
  call extract_laplacian_operator(blocks, A_psi, K_omega)

  write(*,*) "3. Starting time-stepping loop..."
  n_steps = ceiling(t_end / dt)

  do t = 1, n_steps
    d_psi_old = d_psi_new
    d_omega_old = d_omega_new
    current_time = (t - 1) * dt

    ! ====================================================================
    ! L-S RK3 STAGE 1
    ! ====================================================================
    call calculate_advection_term(blocks, d_psi_old, d_omega_old, N_omega_n)
    call assemble_rk3_step(blocks, 1, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)
    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    ! *** FIX: ADD ERROR CHECKING ***
    if (linfo /= 0) then
        write(*,*) "FATAL: DGESV failed in RK3 Step 1 at t=", t, " with code=", linfo
        stop
    endif
    d_omega_s1 = b_step

    ! ====================================================================
    ! L-S RK3 STAGE 2
    ! ====================================================================
    call calculate_advection_term(blocks, d_psi_old, d_omega_s1, N_omega_s1)
    call assemble_rk3_step(blocks, 2, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)
    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    ! *** FIX: ADD ERROR CHECKING ***
    if (linfo /= 0) then
        write(*,*) "FATAL: DGESV failed in RK3 Step 2 at t=", t, " with code=", linfo
        stop
    endif
    d_omega_s2 = b_step

    ! ====================================================================
    ! L-S RK3 STAGE 3
    ! ====================================================================
    call calculate_advection_term(blocks, d_psi_old, d_omega_s2, N_omega_s2)
    call assemble_rk3_step(blocks, 3, A_step, b_step, dt, K_omega, M_omega, &
                           d_omega_old, d_psi_old, d_omega_s1, d_omega_s2, &
                           N_omega_n, N_omega_s1, N_omega_s2)
    call DGESV(N, 1, A_step, N, piv, b_step, N, linfo)
    ! *** FIX: ADD ERROR CHECKING ***
    if (linfo /= 0) then
        write(*,*) "FATAL: DGESV failed in RK3 Step 3 at t=", t, " with code=", linfo
        stop
    endif
    d_omega_new = b_step

    ! ====================================================================
    ! FINAL PSI UPDATE
    ! ====================================================================
    ! *** FIX: Use a fresh copy of A_psi each time to avoid corrupting the original. ***
    A_psi_copy = A_psi ! Use a fresh copy of the clean, original Laplacian matrix
    call apply_psi_bcs(blocks, A_psi_copy, b_psi_bc, t * dt) ! Modify the copy
    b_psi_total = matmul(-M_omega, d_omega_new) + b_psi_bc
    
    call DGESV(N, 1, A_psi_copy, N, piv, b_psi_total, N, linfo)
    ! *** FIX: ADD ERROR CHECKING ***
    if (linfo /= 0) then
        write(*,*) "FATAL: DGESV failed to solve for PSI at t=", t, " with code=", linfo
        stop
    endif
    d_psi_new = b_psi_total

    if (mod(t, output_interval) == 0) then
      write(*,'(A, I6, A, I6, A, F8.4)') "Timestep:", t, "/", n_steps, "  Time:", t * dt
      call output_results(blocks, d_psi_new, d_omega_new, t, "results/solution_")
    endif
  enddo

  write(*,*) "Deallocating memory..."
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