! Defines global constants and data types 
! 09/04/2015
! Ali Yesildag

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
        integer, parameter :: P_x = 3, P_y = 3 ! B-spline degrees

        ! BLOCK 1 (bottom left block)
        integer, parameter :: NX_B1 = 49 ! number of points in x-direction for block 1 
        integer, parameter :: NY_B1 = 16 ! number of points in x-direction for block 1

        ! NOTE: I think NX_B1 must have 7x node and NX_B2_PART2 = 3x based on the ratio
        ! BLOCK 2 (top block)
        integer, parameter :: NX_B2_PART1 = NX_B1 ! points in [0, corner_x]
        integer, parameter :: NX_B2_PART2 = 21    ! points in [corner_x, 1] 
        integer, parameter :: NX_B2 = NX_B2_PART1 + NX_B2_PART2 - 1 
        integer, parameter :: NY_B2 = 24

        ! Physical parameters
        real(dp), parameter :: REYNOLDS_NUMBER = 100.0d0

        ! L-S RK3 Coefficients from Spalart, Moser & Rogers (1991)
        real(dp), dimension(3)  :: rk_alpha = (/ 29.0d0 / 60.0d0, -3.0d0 / 40.0d0, 1.0d0 / 6.0d0 /)
        real(dp), dimension(3)  :: rk_beta  = (/ 37.0d0 / 160.0d0, 5.0d0 / 24.0d0, 1.0d0 / 6.0d0 /)
        real(dp), dimension(3)  :: rk_gamma = (/ 8.0d0  / 15.0d0, 5.0d0 / 12.0d0, 3.0d0 / 4.0d0 /)
        real(dp), dimension(3)  :: rk_zeta  = (/ 0.0d0, -17.0d0 / 60.0d0, -5.0d0 / 12.0d0 /)


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
        ! integer, parameter :: N_x = 100  ! Number of basis functions in x-direction
        ! integer, parameter :: N_y = 100  ! Number of basis functions in y-direction
        ! integer, parameter :: P_x = 5   ! Degree of b-splines in x-direction 
        ! integer, parameter :: P_y = 5   ! Degree of b-splines in y-direction 

        ! Boundary types 
        integer, parameter :: BTYPE_INTERIOR            = 0
        integer, parameter :: BTYPE_MOVING_LID          = 1
        integer, parameter :: BTYPE_WALL                = 2  
        integer, parameter :: BTYPE_INTERFACE           = 3

        ! LAPACK variables 
        integer, allocatable  :: piv(:) ! Pivot indices for dgesv 
        integer               :: linfo  ! status flag from lapack 
        real(dp), allocatable :: A_psi_copy(:,:) ! to save A_psi during RK3
end module mod_params



module mod_bspline

        use mod_params
        implicit none 

contains 

        recursive function bspline_basis(i, p, knots, u) result(val)
        ! calculates the value of a single b-spline basis function using Cox-de Boor recursion formula
        
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




