! =============================================================================
! PROGRAM: check_derivatives
!
! PURPOSE:
!   Verifies the correctness of the B-spline derivative calculations by
!   comparing them against the exact derivatives of a known analytical function.
!
! METHOD:
!   1. An analytical function f(x,y) with known derivatives is chosen.
!   2. The B-spline coefficients that represent this function on the grid are
!      found by solving the linear system M*d = f, where M is the mass matrix.
!   3. The first and second derivatives of the B-spline representation are
!      calculated at each collocation point.
!   4. These approximated derivatives are compared to the exact analytical
!      derivatives, and the L2 norm of the error is computed.
! =============================================================================
PROGRAM check_derivatives
    USE mod_params
    USE mod_grid
    USE mod_bspline
    USE mod_solver
    IMPLICIT NONE

    TYPE(block_type), ALLOCATABLE :: blocks(:)
    INTEGER :: iblock, i, j, k
    INTEGER :: N, info
    INTEGER, ALLOCATABLE :: ipiv(:)
    REAL(dp), ALLOCATABLE :: M(:,:), f_vec(:), d_f(:)
    REAL(dp) :: x, y

    ! Analytical function and its derivatives
    REAL(dp) :: dfdx_exact, dfdy_exact, laplacian_exact

    ! B-spline approximated derivatives
    REAL(dp) :: dfdx_approx, dfdy_approx, laplacian_approx

    ! Error calculation
    REAL(dp) :: l2_err_dfdx, l2_err_dfdy, l2_err_laplacian
    INTEGER :: global_point_offset, row_index

    ! --- INTERFACE BLOCK REMOVED ---
    ! The linker will find the correct interface for the external DGESV subroutine.


    WRITE(*,*) "================================================"
    WRITE(*,*) "   B-Spline Derivative Checker"
    WRITE(*,*) "================================================"

    ! 1. SETUP GRID AND SOLVER
    WRITE(*,*) "1. Initializing grid and solver..."
    ALLOCATE(blocks(NUM_BLOCKS))
    DO iblock = 1, NUM_BLOCKS
        CALL initialize_block(blocks(iblock), iblock)
        CALL generate_grid(blocks(iblock))
    END DO
    CALL initialize_solver(blocks)
    N = (NX_B1 * NY_B1) + (NX_B2 * NY_B2)
    ALLOCATE(M(N, N), f_vec(N), d_f(N), ipiv(N))
    M = 0.0_dp

    WRITE(*,*) "2. Assembling Mass Matrix and Function Vector..."
    global_point_offset = 0
    DO iblock = 1, NUM_BLOCKS
        DO k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            row_index = global_point_offset + k
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            
            ! Evaluate analytical function at the collocation point
            f_vec(row_index) = test_func(x, y)

            ! Assemble the corresponding row in the mass matrix
            DO j = 1, blocks(iblock)%N_y
                DO i = 1, blocks(iblock)%N_x
                    M(row_index, get_global_psi_index(iblock, i, j)) = &
                        bspline_basis_norm(i, blocks(iblock)%P_x, blocks(iblock)%knots_x, x, blocks(iblock)%xmin, blocks(iblock)%xmax) * &
                        bspline_basis_norm(j, blocks(iblock)%P_y, blocks(iblock)%knots_y, y, blocks(iblock)%ymin, blocks(iblock)%ymax)
                END DO
            END DO
        END DO
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    END DO

    WRITE(*,*) "3. Solving for B-spline coefficients..."
    d_f = f_vec
    CALL DGESV(N, 1, M, N, ipiv, d_f, N, info)
    IF (info /= 0) THEN
        WRITE(*,*) "ERROR: DGESV failed to solve the system. INFO =", info
        STOP
    END IF

    WRITE(*,*) "4. Calculating derivatives and errors at each location..."
    l2_err_dfdx = 0.0_dp
    l2_err_dfdy = 0.0_dp
    l2_err_laplacian = 0.0_dp
    global_point_offset = 0

    OPEN(unit=20, file='derivative_check.csv', status='replace')
    WRITE(20, '(A)') '"Block","X","Y","Exact_dFdx","Approx_dFdx","Exact_dFdy","Approx_dFdy","Exact_Lapl","Approx_Lapl"'

    DO iblock = 1, NUM_BLOCKS
        DO k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)

            ! Calculate exact derivatives
            CALL get_exact_derivatives(x, y, dfdx_exact, dfdy_exact, laplacian_exact)
            
            ! Calculate approximated derivatives from B-spline representation
            CALL get_approx_derivatives(blocks(iblock), d_f, x, y, dfdx_approx, dfdy_approx, laplacian_approx)
            
            ! Write data to file
            WRITE(20,'(I2,A,F8.4,A,F8.4,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,F12.6)') &
                iblock, ',', x, ',', y, ',', dfdx_exact, ',', dfdx_approx, ',', &
                dfdy_exact, ',', dfdy_approx, ',', laplacian_exact, ',', laplacian_approx

            ! Accumulate squared errors for L2 norm
            l2_err_dfdx = l2_err_dfdx + (dfdx_approx - dfdx_exact)**2
            l2_err_dfdy = l2_err_dfdy + (dfdy_approx - dfdy_exact)**2
            l2_err_laplacian = l2_err_laplacian + (laplacian_approx - laplacian_exact)**2
        END DO
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
    END DO
    CLOSE(20)

    ! Finalize L2 error norm
    l2_err_dfdx = SQRT(l2_err_dfdx / dble(N))
    l2_err_dfdy = SQRT(l2_err_dfdy / dble(N))
    l2_err_laplacian = SQRT(l2_err_laplacian / dble(N))

    WRITE(*,*) "5. Results:"
    WRITE(*,'(A, E12.5)') "   RMS Error (dF/dx):     ", l2_err_dfdx
    WRITE(*,'(A, E12.5)') "   RMS Error (dF/dy):     ", l2_err_dfdy
    WRITE(*,'(A, E12.5)') "   RMS Error (Laplacian): ", l2_err_laplacian
    WRITE(*,*) "Detailed results saved to 'derivative_check.csv'."
    
    ! CLEAN UP
    DEALLOCATE(blocks, M, f_vec, d_f, ipiv)
    WRITE(*,*) "================================================"

CONTAINS

    !> Analytical test function f(x, y)
    FUNCTION test_func(x, y) RESULT(f_val)
        USE mod_params, ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: x, y
        REAL(dp) :: f_val
        REAL(dp), PARAMETER :: PI = ACOS(-1.0_dp)
        f_val = SIN(PI * x) * COS(PI * y)
    END FUNCTION test_func

    !> Calculates the exact derivatives of the test function
    SUBROUTINE get_exact_derivatives(x, y, dfdx, dfdy, laplacian)
        USE mod_params, ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: x, y
        REAL(dp), INTENT(OUT) :: dfdx, dfdy, laplacian
        REAL(dp), PARAMETER :: PI = ACOS(-1.0_dp)
        
        dfdx = PI * COS(PI * x) * COS(PI * y)
        dfdy = -PI * SIN(PI * x) * SIN(PI * y)
        laplacian = -2.0_dp * PI**2 * SIN(PI * x) * COS(PI * y)
    END SUBROUTINE get_exact_derivatives

    !> Calculates derivatives from the B-spline representation
    SUBROUTINE get_approx_derivatives(block, coeffs, x, y, dfdx, dfdy, laplacian)
        USE mod_params, ONLY: dp, block_type
        USE mod_bspline
        USE mod_solver, ONLY: get_global_psi_index
        IMPLICIT NONE
        TYPE(block_type), INTENT(IN) :: block
        REAL(dp), INTENT(IN) :: coeffs(:)
        REAL(dp), INTENT(IN) :: x, y
        REAL(dp), INTENT(OUT) :: dfdx, dfdy, laplacian

        INTEGER :: i, j, global_idx
        REAL(dp) :: N, M, dN_dx, dM_dy, d2N_dx2, d2M_dy2
        
        dfdx = 0.0_dp; dfdy = 0.0_dp; laplacian = 0.0_dp

        DO j = 1, block%N_y
            DO i = 1, block%N_x
                global_idx = get_global_psi_index(block%id, i, j)

                ! Get basis functions and their physical derivatives
                N       = bspline_basis_norm(i, block%P_x, block%knots_x, x, block%xmin, block%xmax)
                M       = bspline_basis_norm(j, block%P_y, block%knots_y, y, block%ymin, block%ymax)
                dN_dx   = bspline_deriv1_norm(i, block%P_x, block%knots_x, x, block%xmin, block%xmax)
                dM_dy   = bspline_deriv1_norm(j, block%P_y, block%knots_y, y, block%ymin, block%ymax)
                d2N_dx2 = bspline_deriv2_norm(i, block%P_x, block%knots_x, x, block%xmin, block%xmax)
                d2M_dy2 = bspline_deriv2_norm(j, block%P_y, block%knots_y, y, block%ymin, block%ymax)
                
                ! Sum contributions for each derivative
                dfdx = dfdx + coeffs(global_idx) * dN_dx * M
                dfdy = dfdy + coeffs(global_idx) * N * dM_dy
                laplacian = laplacian + coeffs(global_idx) * (d2N_dx2 * M + N * d2M_dy2)
            END DO
        END DO
    END SUBROUTINE get_approx_derivatives

    !> Wrapper for bspline_basis to handle scaling from physical to normalized coordinates
    FUNCTION bspline_basis_norm(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        USE mod_bspline, ONLY: bspline_basis
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u
        u = (x - xmin) / (xmax - xmin)
        val = bspline_basis(i, p, knots, u)
    END FUNCTION bspline_basis_norm

    !> Wrapper for bspline_deriv1 to handle scaling via the chain rule
    FUNCTION bspline_deriv1_norm(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        USE mod_bspline, ONLY: bspline_deriv1
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u, du_dx
        u = (x - xmin) / (xmax - xmin)
        du_dx = 1.0_dp / (xmax - xmin)
        val = bspline_deriv1(i, p, knots, u) * du_dx
    END FUNCTION bspline_deriv1_norm

    !> Wrapper for bspline_deriv2 to handle scaling via the chain rule
    FUNCTION bspline_deriv2_norm(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        USE mod_bspline, ONLY: bspline_deriv2
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u, du_dx
        u = (x - xmin) / (xmax - xmin)
        du_dx = 1.0_dp / (xmax - xmin)
        val = bspline_deriv2(i, p, knots, u) * (du_dx**2)
    END FUNCTION bspline_deriv2_norm

END PROGRAM check_derivatives