module mod_params
        
        implicit none 
        
        ! precision
        integer, parameter :: dp = KIND(1.0d0) ! double precision for now 

        ! GEOMETRY DEFINITION
        real(dp), parameter :: DOMAIN_XMIN = 0.0d0
        real(dp), parameter :: DOMAIN_XMAX = 1.0d0
        real(dp), parameter :: DOMAIN_YMIN = 0.0d0
        real(dp), parameter :: DOMAIN_YMAX = 1.0d0

        ! LOCATION OF THE RE-ENTRANT CORNER
        real(dp), parameter :: CORNER_X = 0.7d0
        real(dp), parameter :: CORNER_Y = 0.4d0

        ! =============================================================================
        ! MESH CONFIGURATION FOR A 3-BLOCK L-SHAPED DOMAIN
        ! =============================================================================
        integer, parameter :: NUM_BLOCKS = 3
        integer, parameter :: P_x = 3, P_y = 3

        ! Block 1 (Bottom-Left)
        integer, parameter :: NX_B1 = 35
        integer, parameter :: NY_B1 = 20

        ! Block 2 (Top-Left)
        integer, parameter :: NX_B2 = NX_B1  ! Must match Block 1 in x
        integer, parameter :: NY_B2 = 30

        ! Block 3 (Right)
        integer, parameter :: NX_B3 = 15
        ! --- THIS IS THE CORRECTED LINE ---
        ! Must match Block 2 in y for alignment.
        integer, parameter :: NY_B3 = NY_B2 

        ! =============================================================================
        ! SHARED DATA STRUCTURE
        ! =============================================================================
        TYPE :: block_type
                integer  :: id, N_x, N_y, P_x, P_y
                real(dp) :: xmin, xmax, ymin, ymax
                real(dp), allocatable :: colloc_pts(:,:), knots_x(:), knots_y(:)
                integer, allocatable  :: boundary_types(:)
        end type block_type

        ! Discretization parameters 
        real(dp), parameter :: stretch_factor = 2.5d0

        ! Boundary types 
        integer, parameter :: BTYPE_INTERIOR            = 0
        integer, parameter :: BTYPE_MOVING_LID          = 1
        integer, parameter :: BTYPE_WALL                = 2  
        integer, parameter :: BTYPE_INTERFACE           = 3

end module mod_params

module mod_bspline

    use mod_params
    implicit none

    private :: bspline_basis, bspline_deriv1, bspline_deriv2

    public :: bspline_basis_physical, bspline_deriv1_physical, bspline_deriv2_physical

contains

    ! Evaluates B-spline basis function in the physical domain [xmin, xmax].
    FUNCTION bspline_basis_physical(i, p, knots, x, xmin, xmax) RESULT(val)
        USE mod_params, ONLY: dp
        INTEGER, INTENT(IN) :: i, p
        REAL(dp), INTENT(IN) :: knots(:), x, xmin, xmax
        REAL(dp) :: val, u
        
        ! Map physical coordinate x to normalized coordinate u
        u = (x - xmin) / (xmax - xmin)
        val = bspline_basis(i, p, knots, u)
    END FUNCTION bspline_basis_physical

    ! Evaluates 1st derivative of B-spline w.r.t. the physical coordinate.
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

    ! Evaluates 2nd derivative of B-spline w.r.t. the physical coordinate.
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



    recursive function bspline_basis(i, p, knots, u) result(val)
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
