! Defines global constants and data types 
! 09/04/2015
! Ali Yesildag

module mod_params
        
        implicit none 
        
        ! precision
        integer, parameter :: dp = KIND(1.0d0) ! double precision for now 

        ! Physical parameters


        ! Discretization parameters 
        integer, parameter :: N_x = 15  ! Number of basis functions in x-direction
        integer, parameter :: N_y = 15  ! Number of basis functions in y-direction
        integer, parameter :: P_x = 3   ! Degree of b-splines in x-direction 
        integer, parameter :: P_y = 3   ! Degree of b-splines in y-direction 

        ! Boundary types 


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

end module mod_bspline

!TODO: bspline_deriv1 and bspline_deriv2 will be needed
! function bspline_deriv1(i, p, knots, u) result(val)
! function bspline_deriv2(i, p, knots, u) result(val)





