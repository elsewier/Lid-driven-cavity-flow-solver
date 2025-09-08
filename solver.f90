module mod_solver 
  use mod_params
  use mod_bspline
  use mod_grid
  implicit none 

  private 
  integer :: num_coeffs_b1 = 0
  integer :: num_coeffs_b2 = 0
  integer :: psi_offset_b1 = 0
  integer :: psi_offset_b2 = 0
  integer :: omega_offset_b1 = 0
  integer :: omega_offset_b2 = 0

  public :: initialize_solver, get_global_psi_index, get_global_omega_index, assemble_constant_matrices, apply_psi_bcs

    
  contains 
    ! calculates and stores the number of coefficients and the starting index/offset for each block
    ! 1 - all psi coeffs from block1 
    ! 2 - all psi coeffs from block2 
    ! 3 - all omega coeffs from block1
    ! 4 - all omega coeffs from block2
    ! {x} = [C_B1 | C_B2 | D_B1 | D_B2]
    subroutine initialize_solver(blocks)
      type(block_type), intent(in) :: blocks(:)
      integer :: total_psi_coeffs

      num_coeffs_b1 = blocks(1)%N_x * blocks(1)%N_y
      num_coeffs_b2 = blocks(2)%N_x * blocks(2)%N_y

      ! block1 psi coeffs starts at 0 offset, block2 start after block1 
      psi_offset_b1 = 0 
      psi_offset_b2 = num_coeffs_b1
      ! block1 omega coeffs start after all psi coeffs 
      omega_offset_b1 = total_psi_coeffs 
      omega_offset_b2 = total_psi_coeffs + num_coeffs_b1
    end subroutine initialize_solver

    function get_global_psi_index(block_id, i, j) result(global_index)
      integer, intent(in) :: block_id, i, j 
      integer             :: global_index
      integer             :: local_index 

      if (block_id == 1) then 
        local_index = (j - 1) * NX_B1 + i ! 2d index to 1d
        global_index = psi_offset_b1 + local_index
      else if (block_id == 2) then 
        local_index = (j - 1) * NX_B2 + i 
        global_index = psi_offset_b2 + local_index
      end if 
    end function get_global_psi_index

      
    function get_global_omega_index(block_id, i, j) result(global_index)
      integer, intent(in) :: block_id, i, j 
      integer             :: global_index
      integer             :: local_index 

      if (block_id == 1) then 
        local_index = (j - 1) * NX_B1 + i ! 2d index to 1d
        global_index = omega_offset_b1 + local_index
      else if (block_id == 2) then 
        local_index = (j - 1) * NX_B2 + i 
        global_index = omega_offset_b2 + local_index
      end if 
    end function get_global_omega_index



    ! calculates the time-independent matrices.
    ! Outputs:
    !   A_psi :  the global laplacian matrix for the streamfunction
    !   M_psi :  the global mass matrix for the streamfunction
    !   M_omega: the global mass matrix for the vorticity 
    subroutine assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(out)         :: A_psi(:,:), M_psi(:,:), M_omega(:,:)

      integer                       :: iblock, k, row_index, i_basis, j_basis 
      integer                       :: col_psi, col_omega 
      real(dp)                      :: x, y
      real(dp)                      :: N_val, M_val, d2N_dx2, d2M_dy2
      integer                       :: global_point_offset

      A_psi = 0.0d0; M_psi = 0.0d0; M_omega = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
          row_index = global_point_offset + k 
          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)

          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
              col_omega = get_global_psi_index(iblock, i_basis, j_basis)
              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)

              M_psi(row_index, col_psi) = N_val * M_val 
              M_omega(row_index, col_omega) = N_val * M_val 
            enddo 
          enddo 

          ! for the A_psi, we only apply the nabla^2 operator at interior points. the boundary rows will be overwritten with BC 
          if (blocks(iblock)%boundary_types(k) == BTYPE_INTERIOR) then 
          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)
              d2N_dx2   = bspline_deriv2(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
              d2M_dy2   = bspline_deriv2(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)

              A_psi(row_index, col_psi) = d2N_dx2 * M_val + N_val * d2M_dy2
            enddo 
          enddo 
        endif
        enddo 
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
      enddo 
    end subroutine assemble_constant_matrices

    subroutine apply_psi_bcs(blocks, A_psi, b_psi_bc, current_time)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(inout)       :: A_psi(:,:) ! matrix to be modified 
      real(dp), intent(out)         :: b_psi_bc(:) ! rhs vector for bcs 
      real(dp), intent(in)          :: current_time 

      integer                       :: iblock, k, row_index, i_basis, j_basis 
      integer                       :: col_psi
      integer                       :: global_point_offset
      real(dp)                      :: x, y
      real(dp)                      :: N_val, M_val, dM_dy 
      real(dp)                      :: u_lid
      

      b_psi_bc = 0.0d0; global_point_offset = 0

      ! loop through all blocks and all points 
      do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
          row_index = global_point_offset + k 
          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)

          select case (blocks(iblock)%boundary_types(k)) ! we only care about boundaries here
          case (BTYPE_INTERIOR)
            cycle 
          
          case (BTYPE_WALL) ! psi = 0 
            A_psi(row_index, :) = 0.0d0 
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
                M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)
                A_psi(row_index, col_psi) = N_val * M_val
              enddo 
            enddo 
            b_psi_bc(row_index) = 0.0d0

          case (BTYPE_MOVING_LID) ! d(psi)/dy = u_lid(x,t)
            ! TODO: CHANGE THE VELOCITY PROFILE LATER FOR NOW I WILL USE SINUSOIDAL VELOCITY
            u_lid = sin(acos(-1.0d0) * x) * (0.5d0 * (1.0d0 - cos(acos(-1.0d0) * current_time)))

            A_psi(row_index,:) = 0.0d0 
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
                dM_dy     = bspline_deriv1(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)
                A_psi(row_index, col_psi) = N_val * M_val
              enddo 
            enddo 
            b_psi_bc(row_index) = u_lid

          case (BTYPE_INTERFACE) ! enforce continuity and derivative equivalance
            A_psi(row_index, :) = 0.0d0 

            if (iblock == 1) then ! psi_1 - psi_2 = 0, coefficient multiplied by +1 
              do j_basis = 1, blocks(1)%N_y 
                do i_basis = 1, blocks(1)%N_x
                  col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x)
                  M_val     = bspline_basis(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * M_val)
                enddo 
              enddo 

              ! coeffs from block2 multiplied by -1
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x
                  col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x)
                  M_val     = bspline_basis(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - (N_val * M_val)
                enddo 
              enddo 

            else if (iblock == 2) then 
              ! we enforce derivative continuity here 
              ! d(psi_B1)/dy - d(psi_B2)/dy = 0
             
              ! again add contibutions from block1 with +1 coeff 
              do j_basis = 1, blocks(1)%N_y 
                do i_basis = 1, blocks(1)%N_x 
                  col_psi   = get_global_psi_index(1, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x)
                  M_val     = bspline_basis(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * M_val) 
                enddo 
              enddo 

              ! negative contribution from block2
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x 
                  col_psi   = get_global_psi_index(2, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x)
                  M_val     = bspline_basis(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - (N_val * M_val) 
                enddo 
              enddo 
            end if 
            b_psi_bc(row_index) = 0.0d0
          end select 
        enddo 
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
      enddo 
    end subroutine apply_psi_bcs

    subroutine calculate_advection_term(blocks, d_psi, d_omega, adv_term_vec)
      ! calculates the non-linear advection term at all interior collocation points (N(omega))
      type(BLOCK_TYPE), intent(in)  :: blocks(:)
      real(dp), intent(in)          :: d_psi(:), d_omega(:)
      real(dp), intent(out)         :: adv_term_vec(:)

      integer   :: iblock, k, i_basis, j_basis
      integer   :: global_psi, global_omega
      integer   :: row_index
      integer   :: global_point_offset
      real(dp)  :: x, y, u_vel, v_vel, domega_dx, domega_dy 
      real(dp)  :: N_val, M_val, dN_dx, dM_dy 

      adv_term_vec = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
        
          ! we only need to calculate this at interior points 
          if (blocks(iblock)%boundary_types(k) /= BTYPE_INTERIOR) CYCLE 

          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)

          u_vel = 0.0d0; v_vel = 0.0d0; domega_dx = 0.0d0; domega_dy = 0.0d0

          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              global_psi    = get_global_psi_index(iblock, i_basis, j_basis)
              global_omega  = get_global_omega_index(iblock, i_basis, j_basis)

              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)
              dN_dx     = bspline_deriv1(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, x)
              dM_dy     = bspline_deriv1(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, y)

              ! calculate velocities: u = d(psi)/dy, v = -d(psi)/dx 
              u_vel = u_vel + d_psi(global_psi) * (N_val * dM_dy)
              v_vel = v_vel - d_psi(global_psi) * (dN_dx * M_val)

              ! calculate vorticity gradients: d(omega)/dx, d(omega)/dy 
              domega_dx = domega_dx + d_omega(global_omega) * (dN_dx * M_val)
              domega_dy = domega_dy + d_omega(global_omega) * (N_val * M_val)
            enddo 
          enddo 

          ! calculate the advection term and store it 
          row_index = global_point_offset + k 
          adv_term_vec(row_index) = -(u_vel * domega_dx + v_vel * domega_dy)
        enddo 
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
      enddo 
    end subroutine calculate_advection_term

              







end module mod_solver
