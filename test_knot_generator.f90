PROGRAM test_knot_generator
    USE mod_params
    ! We need to include the module that contains the subroutine.
    ! When compiling, you would compile grid.f90 first.
    USE mod_grid
    IMPLICIT NONE

    INTEGER, PARAMETER :: P = 3
    INTEGER, PARAMETER :: NUM_BASIS = 69
    REAL(dp), PARAMETER :: BREAK_RATIO = 0.7_dp
    INTEGER :: M
    REAL(dp), ALLOCATABLE :: knots(:)
    REAL(dp), PARAMETER :: TOLERANCE = 1.0E-12_dp

    M = P + NUM_BASIS + 1
    ALLOCATE(knots(M))

    ! Call the subroutine to be tested
    CALL generate_stretched_knots(P, NUM_BASIS, knots, break_point_ratio=BREAK_RATIO)

    ! After stretching, the last internal knot (at index NUM_BASIS) should have a value of exactly 1.0,
    ! because the stretching function maps s=1 to stretched_s=1, and the linear mapping
    ! then maps stretched_s=1 to a final value of 1.0.
    WRITE(*,*) "========================================"
    WRITE(*,*) "     Testing generate_stretched_knots"
    WRITE(*,*) "========================================"
    WRITE(*,'(A, F15.12)') "Value of last internal knot: ", knots(NUM_BASIS)
    WRITE(*,'(A, F15.12)') "Expected value             : ", 1.0_dp

    IF (ABS(knots(NUM_BASIS) - 1.0_dp) < TOLERANCE) THEN
       WRITE(*,*) "RESULT: PASSED"
    ELSE
       WRITE(*,*) "RESULT: FAILED"
       WRITE(*,*) "Reason: The normalization coordinate 's' did not reach 1.0."
    END IF
    WRITE(*,*) "========================================"

    DEALLOCATE(knots)
END PROGRAM test_knot_generator