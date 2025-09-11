program check_poisson
    use mod_params
    use mod_grid
    use mod_bspline
    use mod_solver
    implicit none

    type(block_type), allocatable :: blocks(:)
    integer :: iblock, k, i_basis, j_basis
    integer :: N, info, global_point_offset, row_index, col_idx
    real(dp), parameter :: pi = acos(-1.0d0)
    real(dp) :: x, y

    ! linear system arrays
    real(dp), allocatable :: A_psi(:,:), b_rhs(:), d_psi_numerical(:)
    integer, allocatable  :: ipiv(:)

    ! exact solution arrays
    real(dp), allocatable :: M_psi(:,:), psi_exact_vec(:), d_psi_exact(:)

    ! local variables
    real(dp)  :: N_val, M_val, dM_dy, d2N_dx2, d2M_dy2

    ! Error calculation
    real(dp) :: l2_err_num, l2_err_den, l2_relative_error

    write(*,*) "================================================"
    write(*,*) "   Poisson solver verification (MMS)  "
    write(*,*) "================================================"

    ! 1. SETUP GRID AND SOLVER
    write(*,*) "1. Initializing grid and solver..."
    allocate(blocks(NUM_BLOCKS))
    do iblock = 1, NUM_BLOCKS
        call initialize_block(blocks(iblock), iblock)
        call generate_grid(blocks(iblock))
    enddo
    call initialize_solver(blocks)
    N = (NX_B1 * NY_B1) + (NX_B2 * NY_B2)
    allocate(A_psi(N, N), b_rhs(N), d_psi_numerical(N), ipiv(N))
    A_psi = 0.0d0; b_rhs = 0.0d0;

    write(*,*) "2. Assembling MMS Poisson system..."
    global_point_offset = 0
    do iblock = 1, NUM_BLOCKS
      do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
        row_index = global_point_offset + k
        x = blocks(iblock)%colloc_pts(k, 1)
        y = blocks(iblock)%colloc_pts(k, 2)

        A_psi(row_index, :) = 0.0d0

        select case (blocks(iblock)%boundary_types(k))

        case (BTYPE_INTERIOR)
          do j_basis = 1, blocks(iblock)%N_y; do i_basis = 1, blocks(iblock)%N_x
            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
            N_val   = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
            M_val   = bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
            d2N_dx2 = bspline_deriv2_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
            d2M_dy2 = bspline_deriv2_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
            A_psi(row_index, col_idx) = d2N_dx2 * M_val + N_val * d2M_dy2
          enddo; enddo
          b_rhs(row_index) = f_exact(x, y)

        case (BTYPE_WALL)
          do j_basis = 1, blocks(iblock)%N_y; do i_basis = 1, blocks(iblock)%N_x
            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
            N_val   = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
            M_val   = bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
            A_psi(row_index, col_idx) = N_val * M_val
          enddo; enddo
          b_rhs(row_index) = psi_exact(x, y)

        case (BTYPE_MOVING_LID)
          do j_basis = 1, blocks(iblock)%N_y; do i_basis = 1, blocks(iblock)%N_x
            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
            N_val   = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
            dM_dy   = bspline_deriv1_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
            A_psi(row_index, col_idx) = N_val * dM_dy
          enddo; enddo
          b_rhs(row_index) = dpsi_dy_exact(x, y)

        case (BTYPE_INTERFACE)
          b_rhs(row_index) = 0.0d0
          if (iblock == 1) then
            ! C0 continuity: psi_b1 - psi_b2 = 0
            do j_basis = 1, blocks(1)%N_y; do i_basis = 1, blocks(1)%N_x
              col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
              N_val   = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
              M_val   = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
              A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + (N_val * M_val)
            enddo; enddo
            do j_basis = 1, blocks(2)%N_y; do i_basis = 1, blocks(2)%N_x
              col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
              N_val   = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
              M_val   = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
              A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - (N_val * M_val)
            enddo; enddo
          else if (iblock == 2) then
            ! C1 continuity: d(psi_b1)/dy - d(psi_b2)/dy = 0
            do j_basis = 1, blocks(1)%N_y; do i_basis = 1, blocks(1)%N_x
              col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
              N_val   = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
              dM_dy   = bspline_deriv1_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
              A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + (N_val * dM_dy)
            enddo; enddo
            do j_basis = 1, blocks(2)%N_y; do i_basis = 1, blocks(2)%N_x
              col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
              N_val   = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
              dM_dy   = bspline_deriv1_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
              A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - (N_val * dM_dy)
            enddo; enddo
          endif
        end select
      enddo
      global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    enddo

    write(*,*) "3. Solving for numerical B-spline coefficients..."
    d_psi_numerical = b_rhs
    call dgesv(N, 1, A_psi, N, ipiv, d_psi_numerical, N, info)
    if (info /= 0) then
        write(*,*) "DGESV failed to solve MMS system. INFO=", info
        stop
    endif

    write(*,*) "4. Finding exact B-spline coefficients..."
    allocate(M_psi(N, N), psi_exact_vec(N), d_psi_exact(N))
    M_psi = 0.0d0

    global_point_offset = 0
    do iblock = 1, NUM_BLOCKS; do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
        row_index = global_point_offset + k
        x = blocks(iblock)%colloc_pts(k, 1)
        y = blocks(iblock)%colloc_pts(k, 2)
        do j_basis = 1, blocks(iblock)%N_y; do i_basis = 1, blocks(iblock)%N_x
            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
            N_val   = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
            M_val   = bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
            M_psi(row_index, col_idx) = N_val * M_val
        enddo; enddo
    enddo; global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y); enddo

    global_point_offset = 0
    do iblock = 1, NUM_BLOCKS; do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
        row_index = global_point_offset + k
        x = blocks(iblock)%colloc_pts(k, 1)
        y = blocks(iblock)%colloc_pts(k, 2)
        psi_exact_vec(row_index) = psi_exact(x, y)
    enddo; global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y); enddo

    d_psi_exact = psi_exact_vec
    call dgesv(N, 1, M_psi, N, ipiv, d_psi_exact, N, info)
    if (info /= 0) then
        write(*,*) "DGESV failed to find exact coeffs. INFO=", info
        stop
    endif

    l2_err_num = sum((d_psi_numerical - d_psi_exact)**2)
    l2_err_den = sum(d_psi_exact**2)
    if (l2_err_den == 0.0_dp) then
        l2_relative_error = 0.0_dp
    else
        l2_relative_error = sqrt(l2_err_num / l2_err_den)
    endif

    write(*,*) "5. Results:"
    write(*, '(A, E12.5)') " L2 relative error in coefficients:", l2_relative_error

    deallocate(blocks, A_psi, b_rhs, d_psi_numerical, ipiv)
    deallocate(M_psi, psi_exact_vec, d_psi_exact)

contains

  function psi_exact(x, y) result(val)
    implicit none
    real(dp), intent(in)  :: x, y
    real(dp)              :: val
    val = sin(pi * x) * cos(pi * y)
  end function psi_exact

  function f_exact(x, y) result(val)
    implicit none
    real(dp), intent(in)  :: x, y
    real(dp)              :: val
    val = -2.0d0 * pi**2 * sin(pi * x) * cos(pi * y)
  end function f_exact

  function dpsi_dy_exact(x, y) result(val)
    implicit none
    real(dp), intent(in)  :: x, y
    real(dp)              :: val
    val = -pi * sin(pi * x) * sin(pi * y)
  end function dpsi_dy_exact

end program check_poisson
