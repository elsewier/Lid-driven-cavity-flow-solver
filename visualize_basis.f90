
program visualize_basis
        use mod_params
        use mod_bspline
        use mod_grid
        implicit none 

        integer :: num_basis
        integer :: num_knots
        real(dp), allocatable :: knots(:) 

        integer :: i,j  ! loop counter
        real(dp) :: u   ! the position where we evaluate the functions 
        real(dp) :: basis_value ! calculated value of basis function
        real(dp), parameter :: stretch_factor = 2.5d0


        num_basis = N_x 
        num_knots = num_basis + P_x + 1 

        allocate(knots(num_knots))

        call generate_knot_vector(P_x, num_basis, stretch_factor, knots)
        
        open(11, file='basis_functions.csv', status='replace')
        write(*,*) "Writing data to basis_functions.csv ..."

        ! write(11, '(A, ..., A)')
        write(11, '(A, G0, A)', advance='no') '"u"', ',', '"N_1(u)"'

        do j = 2, num_basis
                write(11, '(A, G0, A)', advance='no') ',', '"N_', j, '(u)"'
        enddo 
        write(11, *)


        do i = 0, 201 - 1
                u = real(i,dp) / real(201 - 1, dp)

        write(11, '(F8.4, A)', advance='no') u, ','

        do j = 1, num_basis
                basis_value = bspline_basis(j, P_x, knots, u)

                if (j < num_basis) then 
                        write(11, '(F8.4, A)', advance='no') basis_value, ','
                else 
                        write(11, '(F8.4)') basis_value
                endif 
        enddo 
        enddo

        close(11)
        deallocate(knots)



end program visualize_basis
