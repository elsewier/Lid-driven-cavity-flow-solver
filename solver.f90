module mod_solver
    use mod_params
    use mod_bspline
    use mod_grid
    implicit none

    private

    integer, allocatable, save :: num_coeffs(:)
    integer, allocatable, save :: psi_offset(:)
    integer, allocatable, save :: omega_offset(:)

    public :: initialize_solver, get_global_psi_index, get_global_omega_index, assemble_constant_matrices, apply_psi_bcs, calculate_advection_term

contains

    ! =============================================================================
    ! GENERALIZED: Calculates and stores coefficient counts and offsets for N blocks.
    ! =============================================================================
    subroutine initialize_solver(blocks)
        type(block_type), intent(in) :: blocks(:)
        integer :: iblock, total_psi_coeffs

        if (allocated(num_coeffs)) deallocate(num_coeffs)
        if (allocated(psi_offset)) deallocate(psi_offset)
        if (allocated(omega_offset)) deallocate(omega_offset)

        allocate(num_coeffs(NUM_BLOCKS))
        allocate(psi_offset(NUM_BLOCKS))
        allocate(omega_offset(NUM_BLOCKS))

        total_psi_coeffs = 0
        do iblock = 1, NUM_BLOCKS
            num_coeffs(iblock) = blocks(iblock)%N_x * blocks(iblock)%N_y
            psi_offset(iblock) = total_psi_coeffs
            total_psi_coeffs = total_psi_coeffs + num_coeffs(iblock)
        end do

        do iblock = 1, NUM_BLOCKS
            omega_offset(iblock) = total_psi_coeffs + psi_offset(iblock)
        end do
    end subroutine initialize_solver

    ! =============================================================================
    ! GENERALIZED: Converts a local (block_id, i, j) to a global psi index.
    ! =============================================================================
    function get_global_psi_index(blocks, block_id, i, j) result(global_index)
        type(block_type), intent(in) :: blocks(:)
        integer, intent(in)          :: block_id, i, j
        integer                      :: global_index, local_index

        local_index = (j - 1) * blocks(block_id)%N_x + i
        global_index = psi_offset(block_id) + local_index
    end function get_global_psi_index

    ! =============================================================================
    ! GENERALIZED: Converts a local (block_id, i, j) to a global omega index.
    ! =============================================================================
    function get_global_omega_index(blocks, block_id, i, j) result(global_index)
        type(block_type), intent(in) :: blocks(:)
        integer, intent(in)          :: block_id, i, j
        integer                      :: global_index, local_index

        local_index = (j - 1) * blocks(block_id)%N_x + i
        global_index = omega_offset(block_id) + local_index
    end function get_global_omega_index

    ! =============================================================================
    ! Assembles the time-independent matrices A_psi, M_psi, and M_omega.
    ! =============================================================================
    subroutine assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
        type(block_type), intent(in)  :: blocks(:)
        real(dp), intent(out)         :: A_psi(:,:), M_psi(:,:), M_omega(:,:)
        integer                       :: iblock, k, row_index, i_basis, j_basis
        integer                       :: col_psi, global_point_offset
        real(dp)                      :: x, y, N_val, M_val, d2N_dx2, d2M_dy2
        type(block_type)              :: current_block

        A_psi = 0.0d0
        M_psi = 0.0d0
        M_omega = 0.0d0
        global_point_offset = 0

        do iblock = 1, NUM_BLOCKS
            current_block = blocks(iblock)
            do k = 1, current_block%N_x * current_block%N_y
                row_index = global_point_offset + k
                x = current_block%colloc_pts(k, 1)
                y = current_block%colloc_pts(k, 2)

                do j_basis = 1, current_block%N_y
                    do i_basis = 1, current_block%N_x
                        col_psi = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                        N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                        M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                        M_psi(row_index, col_psi) = N_val * M_val
                        M_omega(row_index, col_psi) = N_val * M_val
                    end do
                end do

                if (current_block%boundary_types(k) == BTYPE_INTERIOR) then
                    do j_basis = 1, current_block%N_y
                        do i_basis = 1, current_block%N_x
                            col_psi = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            d2N_dx2 = bspline_deriv2_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                            d2M_dy2 = bspline_deriv2_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                            N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                            M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                            A_psi(row_index, col_psi) = d2N_dx2 * M_val + N_val * d2M_dy2
                        end do
                    end do
                end if
            end do
            global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
        end do
    end subroutine assemble_constant_matrices

    ! =============================================================================
    ! Modifies A_psi and builds b_psi_bc for the 3-block geometry.
    ! =============================================================================
    subroutine apply_psi_bcs(blocks, A_psi, b_psi_bc, current_time)
        type(block_type), intent(in)  :: blocks(:)
        real(dp), intent(inout)       :: A_psi(:,:)
        real(dp), intent(out)         :: b_psi_bc(:)
        real(dp), intent(in)          :: current_time
        integer                       :: iblock, k, row_index, i_basis, j_basis
        integer                       :: col_psi, global_point_offset
        real(dp)                      :: x, y, u_lid, N_val, M_val, dN_dx, dM_dy
        type(block_type)              :: current_block

        b_psi_bc = 0.0d0
        global_point_offset = 0

        do iblock = 1, NUM_BLOCKS
            current_block = blocks(iblock)
            do k = 1, current_block%N_x * current_block%N_y
                row_index = global_point_offset + k
                
                if (current_block%boundary_types(k) == BTYPE_INTERIOR) cycle
                
                A_psi(row_index, :) = 0.0d0
                x = current_block%colloc_pts(k, 1)
                y = current_block%colloc_pts(k, 2)
                
                select case (current_block%boundary_types(k))
                case (BTYPE_WALL)
                    do j_basis = 1, current_block%N_y
                        do i_basis = 1, current_block%N_x
                            col_psi = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            A_psi(row_index, col_psi) = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax) &
                                                      * bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                        end do
                    end do
                    b_psi_bc(row_index) = 0.0d0
                case (BTYPE_MOVING_LID)
                    u_lid = 1.0d0
                    do j_basis = 1, current_block%N_y
                        do i_basis = 1, current_block%N_x
                            col_psi = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                            A_psi(row_index, col_psi) = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax) &
                                                      * bspline_deriv1_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                        end do
                    end do
                    b_psi_bc(row_index) = u_lid
                case (BTYPE_INTERFACE)
                    b_psi_bc(row_index) = 0.0d0
                    select case (iblock)
                    case (1) ! -- Block 1 is MASTER (lower ID), enforces C0 --
                        if (abs(x - CORNER_X) < 1.0e-9) then ! Interface with Block 3
                            ! Enforce C0: psi_1 - psi_3 = 0
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_psi = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_psi = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - N_val * M_val
                                end do
                            end do
                        else ! Interface with Block 2
                            ! Enforce C0: psi_1 - psi_2 = 0
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_psi = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_psi = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - N_val * M_val
                                end do
                            end do
                        endif
                    case (2) ! -- Block 2 is MASTER to B3, SLAVE to B1 --
                        if (abs(y - CORNER_Y) < 1.0e-9) then ! Interface with B1. B2 is SLAVE -> C1
                            ! Enforce C1: d(psi_1)/dy - d(psi_2)/dy = 0
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_psi = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    dM_dy = bspline_deriv1_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * dM_dy
                                end do
                            end do
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_psi = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    dM_dy = bspline_deriv1_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - N_val * dM_dy
                                end do
                            end do
                        else ! Interface with B3. B2 is MASTER -> C0
                            ! Enforce C0: psi_2 - psi_3 = 0
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_psi = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_psi = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    N_val = bspline_basis_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - N_val * M_val
                                end do
                            end do
                        endif
                    case (3) ! -- Block 3 is always SLAVE, enforces C1 --
                        if (y < CORNER_Y) then ! Interface with B1
                            ! Enforce C1: d(psi_1)/dx - d(psi_3)/dx = 0
                            do j_basis = 1, blocks(1)%N_y
                                do i_basis = 1, blocks(1)%N_x
                                    col_psi = get_global_psi_index(blocks, 1, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + dN_dx * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_psi = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - dN_dx * M_val
                                end do
                            end do
                        else ! Interface with B2
                            ! Enforce C1: d(psi_2)/dx - d(psi_3)/dx = 0
                            do j_basis = 1, blocks(2)%N_y
                                do i_basis = 1, blocks(2)%N_x
                                    col_psi = get_global_psi_index(blocks, 2, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + dN_dx * M_val
                                end do
                            end do
                            do j_basis = 1, blocks(3)%N_y
                                do i_basis = 1, blocks(3)%N_x
                                    col_psi = get_global_psi_index(blocks, 3, i_basis, j_basis)
                                    dN_dx = bspline_deriv1_physical(i_basis, blocks(3)%P_x, blocks(3)%knots_x, x, blocks(3)%xmin, blocks(3)%xmax)
                                    M_val = bspline_basis_physical(j_basis, blocks(3)%P_y, blocks(3)%knots_y, y, blocks(3)%ymin, blocks(3)%ymax)
                                    A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - dN_dx * M_val
                                end do
                            end do
                        endif
                    end select
                  end select
                end do
            global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
        end do
    end subroutine apply_psi_bcs

    ! =============================================================================
    ! Calculates the non-linear advection term -(u.grad)omega at interior points.
    ! =============================================================================
    subroutine calculate_advection_term(blocks, d_psi, d_omega, adv_term_vec)
        type(BLOCK_TYPE), intent(in) :: blocks(:)
        real(dp), intent(in) :: d_psi(:), d_omega(:)
        real(dp), intent(out) :: adv_term_vec(:)
        integer :: iblock, k, i_basis, j_basis, global_psi, global_omega
        integer :: row_index, global_point_offset
        real(dp) :: x, y, u_vel, v_vel, domega_dx, domega_dy
        real(dp) :: N_val, M_val, dN_dx, dM_dy
        type(block_type) :: current_block
        
        adv_term_vec = 0.0d0
        global_point_offset = 0
        
        do iblock = 1, NUM_BLOCKS
            current_block = blocks(iblock)
            do k = 1, current_block%N_x * current_block%N_y
                if (current_block%boundary_types(k) /= BTYPE_INTERIOR) CYCLE
                
                row_index = global_point_offset + k
                x = current_block%colloc_pts(k, 1)
                y = current_block%colloc_pts(k, 2)
                u_vel = 0.0d0
                v_vel = 0.0d0
                domega_dx = 0.0d0
                domega_dy = 0.0d0
                
                do j_basis = 1, current_block%N_y
                    do i_basis = 1, current_block%N_x
                        global_psi = get_global_psi_index(blocks, iblock, i_basis, j_basis)
                        global_omega = get_global_omega_index(blocks, iblock, i_basis, j_basis)
                        
                        N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                        M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                        dN_dx = bspline_deriv1_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                        dM_dy = bspline_deriv1_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                        
                        u_vel = u_vel + d_psi(global_psi) * (N_val * dM_dy)
                        v_vel = v_vel - d_psi(global_psi) * (dN_dx * M_val)
                        domega_dx = domega_dx + d_omega(global_omega) * (dN_dx * M_val)
                        domega_dy = domega_dy + d_omega(global_omega) * (N_val * dM_dy)
                    end do
                end do
                adv_term_vec(row_index) = - (u_vel * domega_dx + v_vel * domega_dy)
            end do
            global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
        end do
    end subroutine calculate_advection_term

end module mod_solver