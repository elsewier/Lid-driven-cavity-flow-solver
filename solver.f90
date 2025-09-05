module mod_solver 

        use mod_params
        use mod_bspline
        implicit none 
contains 
        
        subroutine assemble_system(A, b, colloc_pts, boundary_types, knots_x, knots_y)

        real(dp), intent(out)   :: A(:,:), b(:)
        real(dp), intent(in)    :: colloc_pts(:,:)
        integer, intent(in)     :: boundary_types(:)
        real(dp), intent(in)    :: knots_(x), knots_y(:)

        ! local vars
        integer                 :: k, row_idx, i_basis, j_basis, col_idx_psi, col_idx_omega
        real(dp)                :: x, y, N_val, M_val, d2N_dx2, d2M_dy2 ! N(i), M(j) 
        integer, parameter      :: NxNy = N_x * N_y

        A = 0.0d0; b = 0.0d0 

        do k = 1, size(boundary_types)
                if (boundary_types(k) == BTYPE_OUT_OF_DOMAIN) CYCLE ! do nothing 
                x       = colloc_pts(k, 1)
                y       = colloc_pts(k, 2)
                row_idx = k

                select case (boundary_types(k))
                case (BTYPE_INTERIOR)
                        do j_basis = 1, N_y     
                                do i_basis = 1, N_x
                                        col_idx_psi     = (j_basis - 1) * N_x + i_basis ! 2d index to 1d index
                                        col_idx_omega   = NxNy + col_idx_psi
                                        N_val           = bspline_basis(i_basis, P_x, knots_x, x)
                                        M_val           = bspline_basis(j_basis, P_y, knots_y, y)
                                        d2N_dx2         = bspline_deriv2(i_basis, P_x, knots_x, x)
                                        d2M_dy2         = bspline_deriv2(j_basis, P_y, knots_y, y)
                                        A(row_idx, col_idx_psi) = d2N_dx2 * M_val + N_val * d2M_dy2
                                        A(row_idx, col_idx_omega) = N_val * M_val
                                enddo
                        enddo
                        b(row_idx) = 0.0d0

                case (BTYPE_WALL)
                        ! TODO
                case (BTYPE_MOVING_LID)
                        ! TODO
                case DEFAULT

                end select 

        enddo

        end subroutine assemble_system




end module mod_solver
