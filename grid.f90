! ======================================================================
! Module mod_grid for a 3-Block L-Shaped Cavity (Corrected Boundary Flagging)
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
    real(dp) :: s, x, y
    
    allocate(x_norm(block%N_x), y_norm(block%N_y))

    ! Generate normalized coordinates [0,1] for both directions
    do i = 1, block%N_x
        s = real(i - 1, dp) / real(block%N_x - 1, dp)
        x_norm(i) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
    enddo
    do j = 1, block%N_y
        s = real(j - 1, dp) / real(block%N_y - 1, dp)
        y_norm(j) = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s - 1.0d0)) + 1.0d0)
    enddo

    ! Force exact boundary values
    x_norm(1) = 0.0_dp
    x_norm(block%N_x) = 1.0_dp
    y_norm(1) = 0.0_dp
    y_norm(block%N_y) = 1.0_dp

    call generate_stretched_knots(block%P_x, block%N_x, block%knots_x)
    call generate_stretched_knots(block%P_y, block%N_y, block%knots_y)

    k = 0
    do j = 1, block%N_y
    do i = 1, block%N_x
        k = k + 1
        
        ! Directly set boundary points to exact values
        if (i == 1) then
            x = block%xmin
        else if (i == block%N_x) then
            x = block%xmax
        else
            x = block%xmin + x_norm(i) * (block%xmax - block%xmin)
        end if
        
        if (j == 1) then
            y = block%ymin
        else if (j == block%N_y) then
            y = block%ymax
        else
            y = block%ymin + y_norm(j) * (block%ymax - block%ymin)
        end if
        
        block%colloc_pts(k, 1) = x
        block%colloc_pts(k, 2) = y

        ! Boundary condition assignment - use coordinates instead of indices
! ======================================================================
!           <<<  DEFINITIVE BOUNDARY FLAGGING LOGIC  >>>
!   Replace your entire "select case (block%id)" block with this one.
! ======================================================================
select case (block%id)
case (1) ! Bottom-Left Block
    ! PRIORITY 1: Check for physical walls.
    ! A point is a wall if it's on the left (i=1), bottom (j=1),
    ! or right (i=N_x, the re-entrant corner) edges.
    if (i == 1 .or. j == 1 .or. i == block%N_x) then
        block%boundary_types(k) = BTYPE_WALL
    ! PRIORITY 2: Only if it's not a wall, check for an interface.
    else if (j == block%N_y) then
        block%boundary_types(k) = BTYPE_INTERFACE
    ! PRIORITY 3: Otherwise, it's an interior point.
    else
        block%boundary_types(k) = BTYPE_INTERIOR
    endif

case (2) ! Top-Left Block
    ! PRIORITY 1: Check for the moving lid.
    if (j == block%N_y) then
        block%boundary_types(k) = BTYPE_MOVING_LID
    ! PRIORITY 2: Check for other physical walls.
    ! The left side (i=1) and the single re-entrant corner point (i=N_x, j=1).
    else if (i == 1 .or. (i == block%N_x .and. j == 1)) then
        block%boundary_types(k) = BTYPE_WALL
    ! PRIORITY 3: Check for interfaces.
    else if (j == 1) then      ! Interface with Block 1
        block%boundary_types(k) = BTYPE_INTERFACE
    else if (i == block%N_x) then ! Interface with Block 3
        block%boundary_types(k) = BTYPE_INTERFACE
    ! PRIORITY 4: Otherwise, it's an interior point.
    else
        block%boundary_types(k) = BTYPE_INTERIOR
    endif

case (3) ! Right Block
    ! PRIORITY 1: Check for the moving lid.
    if (j == block%N_y) then
        block%boundary_types(k) = BTYPE_MOVING_LID
    ! PRIORITY 2: Check for physical walls.
    ! The bottom (j=1) and right (i=N_x) edges.
    else if (j == 1 .or. i == block%N_x) then
        block%boundary_types(k) = BTYPE_WALL
    ! PRIORITY 3: Check for interfaces.
    else if (i == 1) then
        block%boundary_types(k) = BTYPE_INTERFACE
    ! PRIORITY 4: Otherwise, it's an interior point.
    else
        block%boundary_types(k) = BTYPE_INTERIOR
    endif
end select

! Optional but recommended: Add this debug print to confirm the fix
if (block%id == 1 .and. i == 1 .and. j == block%N_y) then
    write(*,*) "DEBUG CHECK: Point (i,j)=", i, j, " at x=", block%colloc_pts(k,1), &
                " was flagged as type ", block%boundary_types(k), "(1=Wall)"
endif
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