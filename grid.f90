! sets up the grid, bcs and assembles the matrix 

module mod_setup
        use mod_params
        use mod_bspline
        implicit none 
contains

        ! generate a uniform knot vector

subroutine generate_knot_vector(p, num_basis, stretching_factor, knots)

        ! inputs
        integer, intent(in)     :: p, num_basis
        real(dp), intent(in)    :: stretching_factor
        real(dp), intent(out)   :: knots(:)

        ! local vars 
        integer                 :: i,m 
        real(dp)                :: s, points

        m = p + num_basis + 1 ! total number of knots

        do i = 1, p + 1
                knots(i) = 0.0d0
        enddo 

        do i = p + 2, num_basis 
                s       = real(i - p -1, dp) / real(num_basis - p, dp) ! this is uniform points right [0,1]
                points  = tanh(stretching_factor * (2.0d0 * s - 1.0d0)) ! to use the full range of tanh we need to scale it to [-1,1]

                knots(i) = 0.5d0 * (points + 1.0d0) ! scale back to [0,1]
        enddo

        do i = num_basis + 1, m
                knots(i) = 1.0d0
        enddo 

end subroutine generate_knot_vector

end module mod_setup

