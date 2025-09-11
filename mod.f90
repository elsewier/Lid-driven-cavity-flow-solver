module mod_params
        
        implicit none 
        
        ! precision
        integer, parameter :: dp = KIND(1.0d0) ! double precision for now 

        ! GEOMETRY DEFINITION
        real(dp), parameter :: DOMAIN_XMIN = 0.0d0
        real(dp), parameter :: DOMAIN_XMAX = 1.0d0
        real(dp), parameter :: DOMAIN_YMIN = 0.0d0
        real(dp), parameter :: DOMAIN_YMAX = 1.0d0

        ! LOCATION OF INNER CORNER 
        real(dp), parameter :: CORNER_X = 0.7d0
        real(dp), parameter :: CORNER_Y = 0.4d0

        ! GRID RESOLUTION - DISCRETIZATION PARAMETERS
        integer, parameter :: NUM_BLOCKS = 2

        ! =============================================================================
        ! CONFIGURATION 1: Degree p=3, Multiplicity m=3 (C^0 Continuity)
        ! =============================================================================
        integer, parameter :: P_x = 3, P_y = 3 ! B-spline degrees
        ! BLOCK 1 (bottom left block)
        integer, parameter :: NX_B1 = 49 ! (49 - 3 - 1) = 45, which is divisible by 3.
        integer, parameter :: NY_B1 = 16 ! (16 - 3 - 1) = 12, which is divisible by 3.
        ! BLOCK 2 (top block)
        integer, parameter :: NX_B2_PART1 = NX_B1
        integer, parameter :: NX_B2_PART2 = 22 ! NX_B2=70. (70-3-1)=66, divisible by 3.
        integer, parameter :: NX_B2 = NX_B2_PART1 + NX_B2_PART2 - 1 
        integer, parameter :: NY_B2 = 25 ! (25-3-1)=21, divisible by 3.
        ! --- Knot Multiplicity ---
        integer, parameter :: MULT_X_B1 = 3, MULT_Y_B1 = 3
        integer, parameter :: MULT_X_B2 = 3, MULT_Y_B2 = 3

        ! =============================================================================
        ! CONFIGURATION 2: Degree p=4, Multiplicity m=3 (C^1 Continuity)
        ! UNCOMMENT THE BLOCK BELOW TO USE
        ! =============================================================================
        ! integer, parameter :: P_x = 4, P_y = 4 ! B-spline degrees
        ! ! BLOCK 1 (bottom left block)
        ! integer, parameter :: NX_B1 = 50 ! (50 - 4 - 1) = 45, which is divisible by 3.
        ! integer, parameter :: NY_B1 = 17 ! (17 - 4 - 1) = 12, which is divisible by 3.
        ! ! BLOCK 2 (top block)
        ! integer, parameter :: NX_B2_PART1 = NX_B1
        ! integer, parameter :: NX_B2_PART2 = 22 ! NX_B2=71. (71-4-1)=66, divisible by 3.
        ! integer, parameter :: NX_B2 = NX_B2_PART1 + NX_B2_PART2 - 1 
        ! integer, parameter :: NY_B2 = 23 ! (23-4-1)=18, divisible by 3.
        ! ! --- Knot Multiplicity ---
        ! integer, parameter :: MULT_X_B1 = 3, MULT_Y_B1 = 3
        ! integer, parameter :: MULT_X_B2 = 3, MULT_Y_B2 = 3

        ! =============================================================================
        ! CONFIGURATION 3: Degree p=5, Multiplicity m=3 (C^2 Continuity)
        ! UNCOMMENT THE BLOCK BELOW TO USE
        ! =============================================================================
        ! integer, parameter :: P_x = 5, P_y = 5 ! B-spline degrees
        ! ! BLOCK 1 (bottom left block)
        ! integer, parameter :: NX_B1 = 52 ! (52 - 5 - 1) = 46 (ERROR). Let's fix. (N-6)%3=0. N must be mult of 3. Use 51.
        ! integer, parameter :: NX_B1 = 51 ! (51 - 5 - 1) = 45, which is divisible by 3.
        ! integer, parameter :: NY_B1 = 18 ! (18 - 5 - 1) = 12, which is divisible by 3.
        ! ! BLOCK 2 (top block)
        ! integer, parameter :: NX_B2_PART1 = NX_B1
        ! integer, parameter :: NX_B2_PART2 = 22 ! NX_B2=72. (72-5-1)=66, divisible by 3.
        ! integer, parameter :: NX_B2 = NX_B2_PART1 + NX_B2_PART2 - 1 
        ! integer, parameter :: NY_B2 = 24 ! (24-5-1)=18, divisible by 3.
        ! ! --- Knot Multiplicity ---
        ! integer, parameter :: MULT_X_B1 = 3, MULT_Y_B1 = 3
        ! integer, parameter :: MULT_X_B2 = 3, MULT_Y_B2 = 3

        ! Physical parameters
        ! ... (rest of the module is unchanged) ...
        TYPE :: block_type                               ! will be used for multiblock mesh 
                integer  :: id                           ! Block id 
                integer  :: N_x, N_y                     ! number of points in each direction
                integer  :: P_x, P_y                     ! B-spline degree in each direction

                real(dp) :: xmin, xmax, ymin, ymax       ! Pyhsical boundaries of the block 

                ! grid and basis information for this block
                real(dp), allocatable :: colloc_pts(:,:)     ! (nx*ny,2)
                integer, allocatable     :: boundary_types(:)   ! (nx * ny)
                real(dp), allocatable    :: knots_x(:), knots_y(:)
        end type block_type

        ! Discretization parameters 
        real(dp), parameter :: stretch_factor = 2.5d0

        ! Boundary types 
        integer, parameter :: BTYPE_INTERIOR            = 0
        integer, parameter :: BTYPE_MOVING_LID          = 1
        integer, parameter :: BTYPE_WALL                = 2  
        integer, parameter :: BTYPE_INTERFACE           = 3

end module mod_params

! =============================================================================
! FILE: mod_bspline.f90
! =============================================================================
module mod_bspline

    use mod_params
    implicit none

    ! The core recursive functions remain private to the module
    private :: bspline_basis, bspline_deriv1, bspline_deriv2

    ! We expose the new physical-domain functions to the user
    public :: bspline_basis_physical, bspline_deriv1_physical, bspline_deriv2_physical

contains

    !> Evaluates B-spline basis function in the physical domain [xmin, xmax].
    FUNCTION bspline_basis_physical(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u
        
        ! Map physical coordinate x to normalized coordinate u
        u = (x - xmin) / (xmax - xmin)
        val = bspline_basis(i, p, knots, u)
    END FUNCTION bspline_basis_physical

    !> Evaluates 1st derivative of B-spline w.r.t. the physical coordinate.
    FUNCTION bspline_deriv1_physical(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u, du_dx
        
        ! Map physical coordinate x to normalized coordinate u
        u = (x - xmin) / (xmax - xmin)
        
        ! Calculate derivative of the mapping for the chain rule: d/dx = (du/dx) * d/du
        du_dx = 1.0_dp / (xmax - xmin)
        
        val = bspline_deriv1(i, p, knots, u) * du_dx
    END FUNCTION bspline_deriv1_physical

    !> Evaluates 2nd derivative of B-spline w.r.t. the physical coordinate.
    FUNCTION bspline_deriv2_physical(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u, du_dx
        
        ! Map physical coordinate x to normalized coordinate u
        u = (x - xmin) / (xmax - xmin)
        
        ! Calculate derivative of the mapping for the chain rule: d^2/dx^2 = (du/dx)^2 * d^2/du^2
        du_dx = 1.0_dp / (xmax - xmin)
        
        val = bspline_deriv2(i, p, knots, u) * (du_dx**2)
    END FUNCTION bspline_deriv2_physical


    ! =========================================================================
    ! CORE RECURSIVE FUNCTIONS (UNCHANGED, NOW PRIVATE)
    ! =========================================================================

    recursive function bspline_basis(i, p, knots, u) result(val)
        ! ... (rest of the function is unchanged)
        ! Inputs
        integer, intent(in)     :: i            ! The index of basis function that we want to calculate
        integer, intent(in)     :: p            ! The polynomial degree of the bspline
        real(dp), intent(in)    :: knots(:)     ! knot values as array
        real(dp), intent(in)    :: u            ! the position where we want to evaluate the function

        ! output
        real(dp)                :: val          ! N_i,p(u)

        ! local vars
        real(dp) :: term1, term2
        real(dp) :: coeff1, coeff2


        ! Base case
        if (p == 0) then
                if (u >= knots(i) .and. u < knots(i + 1)) then
                        val = 1.0_dp ! value is 1 inside this knot
                else
                        val = 0.0_dp ! value is 0 outside this knot.
                end if
                return
        endif

        ! if p is not 0 we should divide it into smaller parts
        term1 = 0.0_dp; term2 = 0.0_dp
        if ((knots(i + p) - knots(i)) > 1.0e-12_dp) then
                coeff1 = (u - knots(i)) / (knots(i + p) - knots(i))
                term1 = coeff1 * bspline_basis(i, p - 1, knots, u)
        endif

        if ((knots(i + p + 1) - knots(i + 1)) > 1.0e-12_dp) then
                coeff2 = (knots(i + p + 1) - u) / (knots(i + p + 1) - knots(i + 1))
                term2 = coeff2 * bspline_basis(i + 1, p - 1, knots, u)
        endif

        val = term1 + term2 ! this is the final result N_i,p(u)

    end function bspline_basis


    function bspline_deriv1(i, p, knots, u) result(val)
        ! ... (function is unchanged)
        integer, intent(in)     :: i,p
        real(dp), intent(in)    :: knots(:), u
        real(dp)                :: val, term1, term2, coeff1, coeff2

        term1 = 0.0d0; term2 = 0.0d0
        if ((knots(i + p) - knots(i)) > 1.0e-12) then
                coeff1 = real(p, dp) / (knots(i + p) - knots(i))
                term1 = coeff1 * bspline_basis(i, p - 1, knots, u)
        endif

        if ((knots(i + p + 1) - knots(i + 1)) > 1.0e-12) then
                coeff2 = real(p, dp) / (knots(i + p + 1) - knots(i + 1))
                term2 = coeff2 * bspline_basis(i + 1, p - 1, knots, u)
        endif

        val = term1 - term2
    end function bspline_deriv1


    function bspline_deriv2(i, p, knots, u) result(val)
        ! ... (function is unchanged)
        integer, intent(in)     :: i,p
        real(dp), intent(in)    :: knots(:), u
        real(dp)                :: val, term1, term2, coeff1, coeff2
        real(dp)                :: deriv1_a, deriv1_b

        term1 = 0.0d0; term2 = 0.0d0
        if ((knots(i + p) - knots(i)) > 1.0e-12) then
                coeff1 = real(p, dp) / (knots(i + p) - knots(i))
                deriv1_a = bspline_deriv1(i, p - 1, knots, u)
                term1 = coeff1 * deriv1_a
        endif

        if ((knots(i + p + 1) - knots(i + 1)) > 1.0e-12) then
                coeff2 = real(p, dp) / (knots(i + p + 1) - knots(i + 1))
                deriv1_b = bspline_deriv1(i + 1, p - 1, knots, u)
                term2 = coeff2 * deriv1_b
        endif

        val = term1 - term2
    end function bspline_deriv2

end module mod_bspline