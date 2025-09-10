module mod_grid
    use mod_params
    use mod_bspline
    implicit none
contains


subroutine initialize_block(block, block_id)
        type(block_type), intent(inout) :: block
        integer, intent(in)             :: block_id

        block%id = block_id
        block%P_x = P_x
        block%P_y = P_y

        select case (block_id)
        case (1)
                block%N_x = NX_B1
                block%N_y = NY_B1
                block%xmin = DOMAIN_XMIN
                block%xmax = CORNER_X
                block%ymin = DOMAIN_YMIN
                block%ymax = CORNER_Y
        case (2)
                block%N_x  = NX_B2
                block%N_y  = NY_B2
                block%xmin = DOMAIN_XMIN
                block%xmax = DOMAIN_XMAX
                block%ymin = CORNER_Y
                block%ymax = DOMAIN_YMAX
        case default
                write(*,*) "Error invalid block_id in initialize_block: ", block_id
                stop
                    
        end select

        allocate(block%colloc_pts(block%N_x * block%N_y, 2))
        allocate(block%boundary_types(block%N_x * block%N_y))
        allocate(block%knots_x(block%N_x + block%P_x + 1))
        allocate(block%knots_y(block%N_y + block%P_y + 1))


end subroutine initialize_block

    subroutine generate_grid(block)
    type(block_type), intent(inout) :: block

    integer :: i, j, k
    real(dp), allocatable :: x_norm(:), y_norm(:) ! Normalized coords [0,1]
    real(dp)              :: s
    real(dp)              :: break_ratio

    allocate(x_norm(block%N_x), y_norm(block%N_y))

    ! =================================================================
    ! CORRECTED X-COORDINATE GENERATION
    ! =================================================================
    ! The piecewise stretching logic was incorrect for the collocation points.
    ! Both blocks should use a smooth hyperbolic tangent stretching across their
    ! entire normalized [0,1] domain. The piecewise stretching of the KNOTS
    ! (handled below) is what ensures the points align at the interface.
    
    ! x direction (now the same simple logic for both blocks)
    do i = 1, block%N_x
        s = real(i - 1, dp) / real(block%N_x - 1, dp)
        x_norm(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
    enddo

    ! y direction (was already correct)
    do j = 1, block%N_y
        s = real(j - 1, dp) / real(block%N_y - 1, dp)
        y_norm(j) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
    enddo

    ! =================================================================
    ! KNOT VECTOR GENERATION (This part was always correct)
    ! =================================================================
    if (block%id == 2) then
        ! Use piecewise stretching for Block 2's x-knots to align with Block 1
        break_ratio = CORNER_X / DOMAIN_XMAX
        call generate_stretched_knots(block%P_x, block%N_x, block%knots_x, break_point_ratio=break_ratio)
    else
        ! Simple stretching for Block 1's x-knots
        call generate_stretched_knots(block%P_x, block%N_x, block%knots_x)
    endif
    ! y knots are simple for both blocks
    call generate_stretched_knots(block%P_y, block%N_y, block%knots_y)


    ! Map to physical domain and flag points (Unchanged)
    k = 0
    do j = 1, block%N_y
    do i = 1, block%N_x
            k = k + 1
            block%colloc_pts(k, 1) = block%xmin + x_norm(i) * (block%xmax - block%xmin)
            block%colloc_pts(k, 2) = block%ymin + y_norm(j) * (block%ymax - block%ymin)

            select case (block%id)
            case (1)
                if (j == block%N_y) then
                    block%boundary_types(k) = BTYPE_INTERFACE
                else if (i == 1 .OR. i == block%N_x .OR. j == 1) then
                    block%boundary_types(k) = BTYPE_WALL
                else
                    block%boundary_types(k) = BTYPE_INTERIOR
               endif
            case (2)
                if (j == 1 .AND. i <= NX_B2_PART1) then
                    block%boundary_types(k) = BTYPE_INTERFACE
                else if (j == block%N_y) then
                    block%boundary_types(k) = BTYPE_MOVING_LID
                else if (i == 1 .OR. i == block%N_x .OR. j == 1) then
                    block%boundary_types(k) = BTYPE_WALL
                else
                    block%boundary_types(k) = BTYPE_INTERIOR
                endif
            end select
    enddo
    enddo

    deallocate(x_norm, y_norm)
end subroutine generate_grid

SUBROUTINE generate_stretched_knots(p, num_basis, knots_out, break_point_ratio)
    integer, intent(in)     :: p, num_basis
    real(dp), intent(out)   :: knots_out(:)
    real(dp), intent(in), optional :: break_point_ratio

    integer  :: i, m, break_index
    integer  :: n_internal, n1, n2
    real(dp) :: s, stretched_s, tanh_val
    real(dp), PARAMETER :: stretch_factor = 2.5d0
    real(dp) :: tanh_norm

    m = p + num_basis + 1
    ! Set clamped boundaries
    knots_out(1 : p + 1) = 0.0_dp
    knots_out(num_basis + 1 : m) = 1.0_dp

    ! Pre-calculate the normalization factor for the stretching function
    tanh_norm = tanh(stretch_factor)
    ! Total number of internal knots to generate
    n_internal = num_basis - p - 1

    if (n_internal <= 0) return ! Nothing to do if there are no internal knots

    if (present(break_point_ratio)) then
        ! --- Two-part stretching for knots ---

        ! Correctly partition the number of internal knots
        n1 = NINT(REAL(n_internal, dp) * break_point_ratio)
        n2 = n_internal - n1
        break_index = p + 1 + n1

        ! Part 1: from 0 to break_point_ratio
        if (n1 > 0) then
            do i = 1, n1
                ! Correct 's' from 0 to 1 over n1 steps
                s = REAL(i - 1, dp) / REAL(n1 - 1, dp)
                if (n1 == 1) s = 0.5_dp ! Handle single-point case to avoid 0/0
                
                tanh_val = tanh(stretch_factor * (2.0d0 * s - 1.0d0))
                stretched_s = (tanh_val + tanh_norm) / (2.0_dp * tanh_norm)
                knots_out(p + 1 + i) = 0.0_dp + stretched_s * break_point_ratio
            enddo
        endif

        ! Part 2: from break_point_ratio to 1
        if (n2 > 0) then
            do i = 1, n2
                ! Correct 's' from 0 to 1 over n2 steps
                s = REAL(i - 1, dp) / REAL(n2 - 1, dp)
                if (n2 == 1) s = 0.5_dp ! Handle single-point case to avoid 0/0

                tanh_val = tanh(stretch_factor * (2.0d0 * s - 1.0d0))
                stretched_s = (tanh_val + tanh_norm) / (2.0_dp * tanh_norm)
                knots_out(break_index + i) = break_point_ratio + stretched_s * (1.0_dp - break_point_ratio)
            enddo
        endif
    else
        ! --- Simple one-part stretching for knots ---
        do i = 1, n_internal
            ! Correct 's' from 0 to 1 over n_internal steps
            s = REAL(i - 1, dp) / REAL(n_internal - 1, dp)
            if (n_internal == 1) s = 0.5_dp ! Handle single-point case to avoid 0/0

            tanh_val = tanh(stretch_factor * (2.0d0 * s - 1.0d0))
            stretched_s = (tanh_val + tanh_norm) / (2.0_dp * tanh_norm)
            knots_out(p + 1 + i) = stretched_s
        enddo
    endif
END SUBROUTINE generate_stretched_knots

end module mod_grid
