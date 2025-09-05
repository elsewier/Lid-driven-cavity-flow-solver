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

subroutine setup_collocation_grid(colloc_pts, boundary_types)
        

        ! outputs
        real(dp), intent(out)   :: colloc_pts(:,:)      ! array of (x,y) coords stored as (Nx*Ny,2)
        integer, intent(out)    :: boundary_types(:)    ! (Nx*Ny) array which stores the type of each point 

        ! local vars
        integer                 :: i, j, k 
        real(dp)                :: x, y 
        real(dp), parameter     :: thr = 1.0e-6         ! will be used while comparing floats
        real(dp)                :: s_x, s_y 
        real(dp), parameter     :: stretch_factor = 2.5d0

        k = 0
        do j = 1, N_y           ! loop over wall-normal direction
                do i = 1, N_x   ! loop over streamwise direction 
                        k = k + 1 ! global point index (1 to Nx*Ny)
                        
                        ! we need to map (i,j) indices to (x,y) coordinates
                        ! since bsplines are constructed using stretched mesh, we need to stretch collocation pts accordingly
                        s_x = real(i - 1, dp) / real(N_x - 1, dp)
                        s_y = real(j - 1, dp) / real(N_y - 1, dp) ! these are uniform between [0, 1]

                        ! stretching
                        if (stretch_factor < 1e-06) then ! no stretch 
                                x = s_x
                                y = s_y
                        else 
                                x = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s_x - 1.0d0)) + 1.0d0)
                                y = 0.5d0 * (tanh(stretch_factor * (2.0d0 * s_y - 1.0d0)) + 1.0d0)
                        endif 
                        
                        colloc_pts(k, 1) = x; colloc_pts(k, 2) = y

                        ! boundary types and bc checks 

                        ! removed part
                        if (x > 0.7d0 + thr .and. y < 0.4d0 + thr) then 
                                boundary_types(k) = BTYPE_OUT_OF_DOMAIN 

                        ! moving lid 
                        else if (abs(y - 1.0d0) < thr) then 
                                boundary_types(k) = BTYPE_MOVING_LID 

                        ! walls 
                        else if ((abs(y - 0.0d0) < thr .and. x < 0.7 + thr)     .or. &          ! bottom wall (x from 0 to 0.7)
                                 (abs(x - 0.0d0) < thr)                         .or. &          ! left wall   (y from 0 to 1)
                                 (abs(x - 1.0d0) < thr .and. y > 0.4d0 - thr)   .or. &          ! right wall  (y from 0.4 to 1)
                                 (abs(x - 0.7d0) < thr .and. y < 0.4d0 + thr)   .or. &          ! inner vertical wall (y 0 to 0.4)
                                 (abs(y - 0.4d0) < thr .and. x > 0.7d0 - thr)) then             ! inner horizontal wall (x 0.7 to 1)
                                boundary_types(k) = BTYPE_WALL 

                        ! interior points 
                        else 
                                boundary_types(k) = BTYPE_INTERIOR 
                        endif 
                enddo 
        enddo 

endsubroutine setup_collocation_grid
        


                        
                        

end module mod_setup

