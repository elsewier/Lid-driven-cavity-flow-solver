! =============================================================================
! PROGRAM: check_grid_and_knots
!
! PURPOSE:
!   Visualizes both the 2D collocation grid and the 1D knot vectors for
!   the multi-block setup. It creates 'grid_data.csv' and 'knots_data.csv'.
! =============================================================================
PROGRAM check_grid_and_knots
    USE mod_params
    USE mod_grid
    IMPLICIT NONE

    TYPE(block_type), ALLOCATABLE :: blocks(:)
    INTEGER :: iblock, k, i
    INTEGER, PARAMETER :: grid_unit = 10, knots_unit = 11
    REAL(dp) :: physical_knot_x
    REAL(dp) :: physical_knot_y

    WRITE(*,*) "================================================"
    WRITE(*,*) "   Multi-Block Grid & Knot Vector Checker"
    WRITE(*,*) "================================================"

    ! 1. ALLOCATE, INITIALIZE, AND GENERATE GRIDS
    ALLOCATE(blocks(NUM_BLOCKS))
    DO iblock = 1, NUM_BLOCKS
        CALL initialize_block(blocks(iblock), iblock)
        CALL generate_grid(blocks(iblock))
    END DO
    WRITE(*,*) "All blocks generated."

    ! 2. WRITE COLLOCATION GRID DATA TO FILE
    OPEN(unit=grid_unit, file='grid_data.csv', status='replace')
    WRITE(*,*) "Writing grid data to grid_data.csv..."
    WRITE(grid_unit, '(A,A,A,A,A,A,A)') '"block_id"', ',', '"x"', ',', '"y"', ',', '"type"'
    DO iblock = 1, NUM_BLOCKS
        DO k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            WRITE(grid_unit, '(I5, A, F10.6, A, F10.6, A, I3)') &
                blocks(iblock)%id,              ',', &
                blocks(iblock)%colloc_pts(k, 1), ',', &
                blocks(iblock)%colloc_pts(k, 2), ',', &
                blocks(iblock)%boundary_types(k)
        END DO
    END DO
    CLOSE(grid_unit)

    ! 3. WRITE KNOT VECTOR DATA TO FILE
    OPEN(unit=knots_unit, file='knots_data.csv', status='replace')
    WRITE(*,*) "Writing knot data to knots_data.csv..."
    WRITE(knots_unit, '(A,A,A,A,A)') '"block_id"', ',', '"direction"', ',', '"knot_value"'
    DO iblock = 1, NUM_BLOCKS
        ! Write x-knots for this block, scaled to its physical domain
        DO i = 1, SIZE(blocks(iblock)%knots_x)
            physical_knot_x = blocks(iblock)%xmin + blocks(iblock)%knots_x(i) * &
                              (blocks(iblock)%xmax - blocks(iblock)%xmin)
            WRITE(knots_unit, '(I5, A, A, A, F10.6)') blocks(iblock)%id, ',', '"x"', ',', physical_knot_x
        END DO
        ! Write y-knots for this block, scaled to its physical domain
        DO i = 1, SIZE(blocks(iblock)%knots_y)
            physical_knot_y = blocks(iblock)%ymin + blocks(iblock)%knots_y(i) * &
                              (blocks(iblock)%ymax - blocks(iblock)%ymin)
            WRITE(knots_unit, '(I5, A, A, A, F10.6)') blocks(iblock)%id, ',', '"y"', ',', physical_knot_y
        END DO
    END DO
    CLOSE(knots_unit)

    ! 4. CLEAN UP
    DO iblock = 1, NUM_BLOCKS
        DEALLOCATE(blocks(iblock)%colloc_pts, blocks(iblock)%boundary_types)
        DEALLOCATE(blocks(iblock)%knots_x, blocks(iblock)%knots_y)
    END DO
    DEALLOCATE(blocks)

    WRITE(*,*) "Done! You can now plot the results."
    WRITE(*,*) "================================================"

END PROGRAM check_grid_and_knots
