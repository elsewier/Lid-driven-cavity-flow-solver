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
        real(dp), allocatable :: p1(:), p2(:)
        real(dp)              :: s
        real(dp)              :: break_ratio

        allocate(x_norm(block%N_x), y_norm(block%N_y))

        ! first generate coordinates 
        ! x direction
        if (block%id == 2) then ! for block2
            ! we will piecewise stretch to align the nodes
            allocate(p1(NX_B2_PART1), p2(NX_B2_PART2))
            do i = 1, NX_B2_PART1
                s = real(i - 1, dp) / real(NX_B2_PART1 - 1, dp)
                p1(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
            enddo
            do i = 1, NX_B2_PART2
                s = real(i - 1, dp) / real(NX_B2_PART2 - 1, dp)
                p2(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
            enddo
            ! Map the two normalized [0,1] arrays to a single normalized [0,1] array with a break
            break_ratio = CORNER_X / DOMAIN_XMAX
            x_norm(1:NX_B2_PART1) = p1 * break_ratio
            x_norm(NX_B2_PART1 + 1 : block%N_x) = break_ratio + p2(2:NX_B2_PART2) * (1.0_dp - break_ratio)

            deallocate(p1, p2)
        else ! block 1
            do i = 1, block%N_x
                s = real(i - 1, dp) / real(block%N_x - 1, dp)
                x_norm(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
            enddo
        endif

        ! y direction
        do j = 1, block%N_y
            s = real(j - 1, dp) / real(block%N_y - 1, dp)
            y_norm(j) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
        enddo

        ! generate knot vectors
        if (block%id == 2) then
            break_ratio = CORNER_X / DOMAIN_XMAX
            ! Create knots with C^0 continuity (multiplicity = p) at the break point
            call generate_stretched_knots(block%P_x, block%N_x, block%knots_x, &
                break_point_ratio=break_ratio, break_knot_multiplicity=block%P_x)
        else
            call generate_stretched_knots(block%P_x, block%N_x, block%knots_x)
        endif
        ! y knots are simple for both blocks
        call generate_stretched_knots(block%P_y, block%N_y, block%knots_y)


        ! convert to physical domain and flag the points
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


SUBROUTINE generate_stretched_knots(p, num_basis, knots_out, break_point_ratio, break_knot_multiplicity)
        integer, intent(in)     :: p, num_basis
        real(dp), intent(out)   :: knots_out(:)
        real(dp), intent(in), optional :: break_point_ratio
        integer, intent(in), optional  :: break_knot_multiplicity !

        integer  :: i, m, break_index, n1, n2, k_mult, n2_new, i_start_part2
        real(dp) :: s, stretched_s
        real(dp), PARAMETER :: stretch_factor = 2.5d0

        m = p + num_basis + 1
        knots_out(1 : p + 1) = 0.0_dp
        knots_out(num_basis + 1 : m) = 1.0_dp

        if (present(break_point_ratio)) then

            ! Set multiplicity. Default is 1 (a simple knot)
            k_mult = 1
            if (present(break_knot_multiplicity)) then
                k_mult = break_knot_multiplicity
            endif

            ! index for the breakpoint(corner)
            break_index = NINT(REAL(num_basis - p, dp) * break_point_ratio) + p + 1
            n1 = break_index - (p + 1)
            n2 = num_basis - break_index + 1 

            ! --- Part 1: Stretch from 0.0 to break_point_ratio 
            do i = p + 2, break_index
                s = REAL(i - (p + 1), dp) / REAL(n1, dp)
                stretched_s = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
                knots_out(i) = 0.0_dp + stretched_s * break_point_ratio
            enddo
            ! At this point, knots_out(break_index) = break_point_ratio

            ! --- Part 2: Add multiplicity ---
            ! Add (k_mult - 1) more knots at the same location
            ! This loop runs (k_mult - 1) times. If k_mult=1, it doesn't run.
            do i = break_index + 1, break_index + k_mult - 1
                knots_out(i) = break_point_ratio
            enddo

            ! --- Part 3: Stretch from break_point_ratio to 1.0 ---
            ! This loop must now start *after* the repeated knots
            i_start_part2 = break_index + k_mult
            
            ! Adjust original denominator 'n2' for knots "used" by multiplicity
            n2_new = n2 - (k_mult - 1) 
            
            if (n2_new > 0) then
                ! There are still knot slots left to fill with a stretch
                do i = i_start_part2, num_basis
                    ! Offset numerator to start from 1 (relative to the new start)
                    s = REAL(i - (i_start_part2 - 1), dp) / REAL(n2_new, dp)
                    stretched_s = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
                    knots_out(i) = break_point_ratio + stretched_s * (1.0_dp - break_point_ratio)
                enddo
            else
                ! Set all remaining internal knots to the break_point_ratio.
                do i = i_start_part2, num_basis
                    knots_out(i) = break_point_ratio
                enddo
            endif

        else
            ! --- Simple one-part stretching 
            do i = p + 2, num_basis
                s = REAL(i - p - 1, dp) / REAL(num_basis - p, dp)
                knots_out(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
            enddo
        endif
end subroutine generate_stretched_knots

end module mod_grid
