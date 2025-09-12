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

    ! Linear system arrays
    real(dp), allocatable :: A_psi(:, :), b_rhs(:), d_psi_numerical(:)
    real(dp), allocatable :: A_psi_copy(:, :) ! For residual check
    integer, allocatable :: ipiv(:)

    ! Exact solution arrays
    real(dp), allocatable :: M_psi(:, :), psi_exact_vec(:), d_psi_exact(:)
    real(dp), allocatable :: M_psi_copy(:, :) ! For residual check

    ! Local variables
    real(dp) :: x, y, N_val, M_val, dN_dx, dM_dy, d2N_dx2, d2M_dy2

    ! Pointwise error calculation variables
    real(dp) :: psi_numerical_at_point, psi_exact_at_point, error_at_point
    real(dp) :: l2_err_num_phys, l2_err_den_phys, l2_relative_error_phys

    ! Residual check variables
    real(dp), allocatable :: res_vec(:)
    real(dp) :: res_norm

    write(*, *) "================================================"
    write(*, *) "   Poisson solver verification (MMS) - 3 Blocks "
    write(*, *) "================================================"

    ! 1. SETUP GRID AND SOLVER
    write(*, *) "1. Initializing grid and solver..."
    allocate(blocks(NUM_BLOCKS))
    do iblock = 1, NUM_BLOCKS
        call initialize_block(blocks(iblock), iblock)
        call generate_grid(blocks(iblock))
    end do
    call initialize_solver(blocks)

    N = 0
    do iblock = 1, NUM_BLOCKS
        N = N + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    end do
    
    allocate(A_psi(N, N), b_rhs(N), d_psi_numerical(N), ipiv(N))

    ! 2. ASSEMBLE MMS POISSON SYSTEM
    write(*, *) "2. Assembling MMS Poisson System..."
    global_point_offset = 0
    A_psi = 0.0d0
    b_rhs = 0.0d0
    do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            row_index = global_point_offset + k
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            A_psi(row_index, :) = 0.0d0

            select case (blocks(iblock)%boundary_types(k))
                case (BTYPE_INTERIOR)
                    do j_basis = 1, blocks(iblock)%N_y
                        do i_basis = 1, blocks(iblock)%N_x
                            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            d2N_dx2 = bspline_deriv2_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
                            d2M_dy2 = bspline_deriv2_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                            N_val = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax)
                            M_val = bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                            A_psi(row_index, col_idx) = d2N_dx2 * M_val + N_val * d2M_dy2
                        end do
                    end do
                    b_rhs(row_index) = f_exact(x, y)
                case (BTYPE_WALL)
                    do j_basis = 1, blocks(iblock)%N_y
                        do i_basis = 1, blocks(iblock)%N_x
                            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            A_psi(row_index, col_idx) = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax) * &
                                                        bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                        end do
                    end do
                    b_rhs(row_index) = psi_exact(x, y)
                case (BTYPE_MOVING_LID)
                    do j_basis = 1, blocks(iblock)%N_y
                        do i_basis = 1, blocks(iblock)%N_x
                            col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            A_psi(row_index, col_idx) = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax) * &
                                                        bspline_deriv1_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                        end do
                    end do
                    b_rhs(row_index) = dpsi_dy_exact(x, y)
                
                case (BTYPE_INTERFACE)
                    b_rhs(row_index) = 0.0d0
                    select case (iblock)
                    case (1) ! -- Block 1 is MASTER (lower ID), enforces C0 --
                        if (abs(x - CORNER_X) < 1.0e-9) then ! Interface with Block 3
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_idx = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - N_val * M_val
                                end do
                            end do
                        else ! Interface with Block 2
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - N_val * M_val
                                end do
                            end do
                        endif
                    case (2) ! -- Block 2 can be MASTER or SLAVE --
                        if (abs(y - CORNER_Y) < 1.0e-9) then ! Interface with B1. B2 is SLAVE, enforces C1.
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    dM_dy = bspline_deriv1_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + N_val * dM_dy
                                end do
                            end do
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    dM_dy = bspline_deriv1_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - N_val * dM_dy
                                end do
                            end do
                        else ! Interface with B3. B2 is MASTER, enforces C0.
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_idx = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - N_val * M_val
                                end do
                            end do
                        endif
                    case (3) ! -- Block 3 is always SLAVE, enforces C1 --
                        if (y < CORNER_Y) then ! Interface with B1
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_idx = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + dN_dx * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_idx = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - dN_dx * M_val
                                end do
                            end do
                        else ! Interface with B2
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_idx = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) + dN_dx * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_idx = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_idx) = A_psi(row_index, col_idx) - dN_dx * M_val
                                end do
                            end do
                        endif
                    end select
            end select
        end do
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    end do
    
    ! 3. SOLVE
    write(*, *) "3. Solving for numerical B-spline coefficients..."
    d_psi_numerical = b_rhs
    allocate(A_psi_copy(N, N)); A_psi_copy = A_psi
    call dgesv(N, 1, A_psi, N, ipiv, d_psi_numerical, N, info)
    if (info /= 0) then
        write(*, *) "DGESV failed. INFO=", info
        stop
    end if
    deallocate(A_psi)

    ! --- Residual Check ---
    allocate(res_vec(N))
    res_vec = matmul(A_psi_copy, d_psi_numerical) - b_rhs
    if (sum(b_rhs**2) > 1.0e-30_dp) then
        res_norm = sqrt(sum(res_vec**2)) / sqrt(sum(b_rhs**2))
    else
        res_norm = sqrt(sum(res_vec**2))
    end if
    write(*, '(A, E12.5)') "   Relative Residual ||A*d_num - b||/||b||: ", res_norm
    deallocate(res_vec, A_psi_copy)

    ! 4. FIND EXACT COEFFICIENTS
    write(*, *) "4. Finding exact B-spline coefficients..."
    allocate(M_psi(N, N), psi_exact_vec(N), d_psi_exact(N))
    M_psi = 0.0d0
    global_point_offset = 0
    do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            row_index = global_point_offset + k
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            do j_basis = 1, blocks(iblock)%N_y
                do i_basis = 1, blocks(iblock)%N_x
                    col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                    M_psi(row_index, col_idx) = bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax) * &
                                                bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                end do
            end do
        end do
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    end do
    global_point_offset = 0
    do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            row_index = global_point_offset + k
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            psi_exact_vec(row_index) = psi_exact(x, y)
        end do
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    end do
    d_psi_exact = psi_exact_vec
    allocate(M_psi_copy(N, N)); M_psi_copy = M_psi
    call dgesv(N, 1, M_psi, N, ipiv, d_psi_exact, N, info)
    if (info /= 0) then
        write(*, *) "DGESV failed for exact coeffs. INFO=", info
        stop
    end if
    
    ! --- Residual Check ---
    allocate(res_vec(N))
    res_vec = matmul(M_psi_copy, d_psi_exact) - psi_exact_vec
    if (sum(psi_exact_vec**2) > 1.0e-30_dp) then
        res_norm = sqrt(sum(res_vec**2)) / sqrt(sum(psi_exact_vec**2))
    else
        res_norm = sqrt(sum(res_vec**2))
    end if
    write(*, '(A, E12.5)') "   Relative Residual ||M*d_exc - p||/||p||: ", res_norm
    deallocate(res_vec, M_psi_copy)

    ! 5. CALCULATE ERROR
    write(*, *) "5. Calculating pointwise error and writing to file..."
    open(unit=30, file='poisson_error.csv', status='replace')
    write(30, '(A)') '"x","y","psi_numerical","psi_exact","error"'
    global_point_offset = 0
    l2_err_num_phys = 0.0_dp
    l2_err_den_phys = 0.0_dp
    do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            row_index = global_point_offset + k
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            psi_numerical_at_point = 0.0_dp
            do j_basis = 1, blocks(iblock)%N_y
                do i_basis = 1, blocks(iblock)%N_x
                    col_idx = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                    psi_numerical_at_point = psi_numerical_at_point + d_psi_numerical(col_idx) * &
                                             bspline_basis_physical(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax) * &
                                             bspline_basis_physical(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                end do
            end do
            psi_exact_at_point = psi_exact(x, y)
            error_at_point = psi_numerical_at_point - psi_exact_at_point
            l2_err_num_phys = l2_err_num_phys + error_at_point**2
            l2_err_den_phys = l2_err_den_phys + psi_exact_at_point**2
            write(30, '(F8.4, A, F8.4, A, E15.7, A, E15.7, A, E15.7)') x, ',', y, ',', psi_numerical_at_point, ',', psi_exact_at_point, ',', error_at_point
        end do
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    end do
    close(30)
    write(*, *) "   ... Done. Data saved to poisson_error.csv"

    ! 6. RESULTS
    if (l2_err_den_phys < 1.0e-30_dp) then
        l2_relative_error_phys = sqrt(l2_err_num_phys)
    else
        l2_relative_error_phys = sqrt(l2_err_num_phys / l2_err_den_phys)
    end if
    write(*, *) "6. Results:"
    write(*, '(A, E12.5)') " L2 relative error in solution values:", l2_relative_error_phys

    ! 7. CLEAN UP
    deallocate(blocks, b_rhs, d_psi_numerical, ipiv)
    deallocate(M_psi, psi_exact_vec, d_psi_exact)
    write(*,*) "================================================"

contains

    function psi_exact(x, y) result(val)
        implicit none
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = sin(pi * x) * cos(pi * y)
    end function psi_exact

    function f_exact(x, y) result(val)
        implicit none
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = -2.0d0 * pi**2 * sin(pi * x) * cos(pi * y)
    end function f_exact

    function dpsi_dy_exact(x, y) result(val)
        implicit none
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = -pi * sin(pi * x) * sin(pi * y)
    end function dpsi_dy_exact

end program check_poisson
