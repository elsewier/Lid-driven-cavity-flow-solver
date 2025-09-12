! ======================================================================
! Module mod_grid for a 3-Block L-Shaped Cavity (FINAL, SIMPLIFIED)
! ======================================================================
module mod_grid
    use mod_params
    implicit none
    
contains

    subroutine initialize_block(block, block_id)
        type(block_type), intent(inout) :: block
        integer, intent(in)             :: block_id

        block%id = block_id; block%P_x = P_x; block%P_y = P_y

        select case (block_id)
        case (1) ! Bottom-Left
            block%N_x = NX_B1; block%N_y = NY_B1
            block%xmin = DOMAIN_XMIN; block%xmax = CORNER_X
            block%ymin = DOMAIN_YMIN; block%ymax = CORNER_Y
        case (2) ! Top-Left
            block%N_x = NX_B2; block%N_y = NY_B2
            block%xmin = DOMAIN_XMIN; block%xmax = CORNER_X
            block%ymin = CORNER_Y; block%ymax = DOMAIN_YMAX
        case (3) ! Right
            block%N_x = NX_B3; block%N_y = NY_B3
            block%xmin = CORNER_X; block%xmax = DOMAIN_XMAX
            block%ymin = CORNER_Y; block%ymax = DOMAIN_YMAX
        case default
            write(*,*) "Error: invalid block_id in initialize_block: ", block_id
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
        real(dp), allocatable :: x_norm(:), y_norm(:)
        real(dp) :: s
        
        allocate(x_norm(block%N_x), y_norm(block%N_y))

        ! Generate normalized coordinates [0,1] for both directions.
        ! The complex piecewise logic is no longer needed.
        do i = 1, block%N_x
            s = real(i - 1, dp) / real(block%N_x - 1, dp)
            x_norm(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
        enddo
        do j = 1, block%N_y
            s = real(j - 1, dp) / real(block%N_y - 1, dp)
            y_norm(j) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
        enddo

        call generate_stretched_knots(block%P_x, block%N_x, block%knots_x)
        call generate_stretched_knots(block%P_y, block%N_y, block%knots_y)

        k = 0
        do j = 1, block%N_y
        do i = 1, block%N_x
            k = k + 1
            block%colloc_pts(k, 1) = block%xmin + x_norm(i) * (block%xmax - block%xmin)
            block%colloc_pts(k, 2) = block%ymin + y_norm(j) * (block%ymax - block%ymin)

            select case (block%id)
            case (1) ! Bottom-Left
                if (i == 1 .or. j == 1) then
                    block%boundary_types(k) = BTYPE_WALL
                else if (j == block%N_y) then
                    block%boundary_types(k) = BTYPE_INTERFACE
                else
                    block%boundary_types(k) = BTYPE_INTERIOR
                endif
            case (2) ! Top-Left
                if (i == 1) then
                    block%boundary_types(k) = BTYPE_WALL
                else if (j == block%N_y) then
                    block%boundary_types(k) = BTYPE_MOVING_LID
                else if (i == block%N_x .or. j == 1) then
                    block%boundary_types(k) = BTYPE_INTERFACE
                else
                    block%boundary_types(k) = BTYPE_INTERIOR
                endif
            case (3) ! Right
                if (i == block%N_x .or. j == 1) then
                    block%boundary_types(k) = BTYPE_WALL
                else if (j == block%N_y) then
                    block%boundary_types(k) = BTYPE_MOVING_LID
                else if (i == 1) then
                    block%boundary_types(k) = BTYPE_INTERFACE
                else
                    block%boundary_types(k) = BTYPE_INTERIOR
                endif
            end select
        enddo
        enddo
        deallocate(x_norm, y_norm)
    end subroutine generate_grid

    subroutine generate_stretched_knots(p, num_basis, knots_out)
        integer, intent(in)     :: p, num_basis
        real(dp), intent(out)   :: knots_out(:)
        integer  :: i, m
        real(dp) :: s
        
        m = p + num_basis + 1
        knots_out(1 : p + 1) = 0.0_dp
        knots_out(num_basis + 1 : m) = 1.0_dp
        do i = p + 2, num_basis
            s = REAL(i - p - 1, dp) / REAL(num_basis - p, dp)
            knots_out(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
        enddo
    end subroutine generate_stretched_knots

end module mod_grid
