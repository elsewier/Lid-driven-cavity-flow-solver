! =============================================================================
! FILE: solver.f90 (Corrected)
! =============================================================================
module mod_solver 
  use mod_params
  use mod_bspline
  use mod_grid
  implicit none 

  private 
  
  ! Module-level variables to hold solver configuration
  integer, save :: num_coeffs_b1 = 0
  integer, save :: num_coeffs_b2 = 0
  integer, save :: psi_offset_b1 = 0
  integer, save :: psi_offset_b2 = 0
  integer, save :: omega_offset_b1 = 0
  integer, save :: omega_offset_b2 = 0
  
  ! NEW: Module-level array to store the number of basis functions in x for each block.
  ! This solves the "blocks_info" compiler error.
  integer, allocatable, private, save :: nx_per_block(:)

  public :: initialize_solver, get_global_psi_index, get_global_omega_index, assemble_constant_matrices, apply_psi_bcs, calculate_advection_term
    
  contains 
    
    !> Calculates and stores the number of coefficients and the starting index/offset for each block.
    !> Also stores the N_x dimension of each block for use by the index functions.
    subroutine initialize_solver(blocks)
      type(block_type), intent(in) :: blocks(:)
      integer :: total_psi_coeffs, num_blocks_in, i

      num_blocks_in = size(blocks)
      
      ! Allocate and store the N_x values
      if (allocated(nx_per_block)) deallocate(nx_per_block)
      allocate(nx_per_block(num_blocks_in))
      do i = 1, num_blocks_in
        nx_per_block(i) = blocks(i)%N_x
      end do

      num_coeffs_b1 = blocks(1)%N_x * blocks(1)%N_y
      num_coeffs_b2 = blocks(2)%N_x * blocks(2)%N_y
      total_psi_coeffs = num_coeffs_b1 + num_coeffs_b2

      psi_offset_b1 = 0 
      psi_offset_b2 = num_coeffs_b1
      omega_offset_b1 = total_psi_coeffs 
      omega_offset_b2 = total_psi_coeffs + num_coeffs_b1
    end subroutine initialize_solver

    !> Converts a local (block_id, i, j) index to a global index for the psi vector.
    function get_global_psi_index(block_id, i, j) result(global_index)
      integer, intent(in) :: block_id, i, j 
      integer             :: global_index, local_index 

      if (block_id == 1) then
        ! CORRECTED: Use the stored N_x value for the block
        local_index = (j - 1) * nx_per_block(1) + i
        global_index = psi_offset_b1 + local_index
      else if (block_id == 2) then
        ! CORRECTED: Use the stored N_x value for the block
        local_index = (j - 1) * nx_per_block(2) + i
        global_index = psi_offset_b2 + local_index
      end if 
    end function get_global_psi_index

    !> Converts a local (block_id, i, j) index to a global index for the omega vector.
    function get_global_omega_index(block_id, i, j) result(global_index)
      integer, intent(in) :: block_id, i, j 
      integer             :: global_index, local_index 

      if (block_id == 1) then
        ! CORRECTED: Use the stored N_x value for the block
        local_index = (j - 1) * nx_per_block(1) + i
        global_index = omega_offset_b1 + local_index
      else if (block_id == 2) then
        ! CORRECTED: Use the stored N_x value for the block
        local_index = (j - 1) * nx_per_block(2) + i
        global_index = omega_offset_b2 + local_index
      end if 
    end function get_global_omega_index


    !> Assembles the time-independent matrices A_psi, M_psi, and M_omega.
    subroutine assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(out)         :: A_psi(:,:), M_psi(:,:), M_omega(:,:)
      ! ... (This subroutine is now correct from the previous response) ...
      integer                       :: iblock, k, row_index, i_basis, j_basis 
      integer                       :: col_psi, col_omega 
      real(dp)                      :: x, y
      real(dp)                      :: N_val, M_val, d2N_dx2, d2M_dy2
      integer                       :: global_point_offset
      type(block_type)              :: current_block

      A_psi = 0.0d0; M_psi = 0.0d0; M_omega = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        current_block = blocks(iblock)
        do k = 1, current_block%N_x * current_block%N_y 
          row_index = global_point_offset + k 
          x = current_block%colloc_pts(k, 1)
          y = current_block%colloc_pts(k, 2)

          do j_basis = 1, current_block%N_y 
            do i_basis = 1, current_block%N_x
              col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
              col_omega = col_psi 
              N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
              M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
              M_psi(row_index, col_psi)     = N_val * M_val 
              M_omega(row_index, col_omega) = N_val * M_val 
            enddo 
          enddo 

          if (current_block%boundary_types(k) == BTYPE_INTERIOR) then 
            do j_basis = 1, current_block%N_y 
              do i_basis = 1, current_block%N_x
                col_psi = get_global_psi_index(iblock, i_basis, j_basis)
                N_val   = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                M_val   = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                d2N_dx2 = bspline_deriv2_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                d2M_dy2 = bspline_deriv2_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                A_psi(row_index, col_psi) = d2N_dx2 * M_val + N_val * d2M_dy2
              enddo 
            enddo 
          endif
        enddo 
        global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
      enddo 
    end subroutine assemble_constant_matrices

    !> Modifies A_psi and builds b_psi_bc to enforce boundary conditions for the streamfunction.
    subroutine apply_psi_bcs(blocks, A_psi, b_psi_bc, current_time)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(inout)       :: A_psi(:,:) 
      real(dp), intent(out)         :: b_psi_bc(:) 
      real(dp), intent(in)          :: current_time 
      ! ... (This subroutine is now correct from the previous response) ...
      integer                       :: iblock, k, row_index, i_basis, j_basis 
      integer                       :: col_psi
      integer                       :: global_point_offset
      real(dp)                      :: x, y
      real(dp)                      :: N_val, M_val, dM_dy 
      real(dp)                      :: u_lid
      type(block_type)              :: current_block

      b_psi_bc = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        current_block = blocks(iblock)
        do k = 1, current_block%N_x * current_block%N_y 
          row_index = global_point_offset + k 
          x = current_block%colloc_pts(k, 1)
          y = current_block%colloc_pts(k, 2)

          select case (current_block%boundary_types(k))
          case (BTYPE_INTERIOR)
            cycle 
          case (BTYPE_WALL)
            A_psi(row_index, :) = 0.0d0 
            do j_basis = 1, current_block%N_y 
              do i_basis = 1, current_block%N_x
                col_psi = get_global_psi_index(iblock, i_basis, j_basis)
                N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                A_psi(row_index, col_psi) = N_val * M_val
              enddo 
            enddo 
            b_psi_bc(row_index) = 0.0d0
          case (BTYPE_MOVING_LID)
            u_lid = sin(acos(-1.0d0) * x) * (0.5d0 * (1.0d0 - cos(acos(-1.0d0) * current_time)))
            A_psi(row_index,:) = 0.0d0 
            do j_basis = 1, current_block%N_y 
              do i_basis = 1, current_block%N_x
                col_psi = get_global_psi_index(iblock, i_basis, j_basis)
                N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
                dM_dy = bspline_deriv1_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
                A_psi(row_index, col_psi) = N_val * dM_dy
              enddo 
            enddo 
            b_psi_bc(row_index) = u_lid
          case (BTYPE_INTERFACE)
            A_psi(row_index, :) = 0.0d0 
            b_psi_bc(row_index) = 0.0d0
            if (iblock == 1) then 
              do j_basis = 1, blocks(1)%N_y 
                do i_basis = 1, blocks(1)%N_x
                  col_psi = get_global_psi_index(1, i_basis, j_basis)
                  N_val   = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                  M_val   = bspline_basis_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * M_val)
                enddo 
              enddo 
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x
                  col_psi = get_global_psi_index(2, i_basis, j_basis)
                  N_val   = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                  M_val   = bspline_basis_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - (N_val * M_val)
                enddo 
              enddo 
            else if (iblock == 2) then
              do j_basis = 1, blocks(1)%N_y 
                do i_basis = 1, blocks(1)%N_x 
                  col_psi = get_global_psi_index(1, i_basis, j_basis)
                  N_val   = bspline_basis_physical(i_basis, blocks(1)%P_x, blocks(1)%knots_x, x, blocks(1)%xmin, blocks(1)%xmax)
                  dM_dy   = bspline_deriv1_physical(j_basis, blocks(1)%P_y, blocks(1)%knots_y, y, blocks(1)%ymin, blocks(1)%ymax)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * dM_dy) 
                enddo 
              enddo 
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x 
                  col_psi = get_global_psi_index(2, i_basis, j_basis)
                  N_val   = bspline_basis_physical(i_basis, blocks(2)%P_x, blocks(2)%knots_x, x, blocks(2)%xmin, blocks(2)%xmax)
                  dM_dy   = bspline_deriv1_physical(j_basis, blocks(2)%P_y, blocks(2)%knots_y, y, blocks(2)%ymin, blocks(2)%ymax)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - (N_val * dM_dy) 
                enddo 
              enddo 
            end if 
          end select 
        enddo 
        global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
      enddo 
    end subroutine apply_psi_bcs

    !> Calculates the non-linear advection term -(u.grad)omega at all interior points.
    subroutine calculate_advection_term(blocks, d_psi, d_omega, adv_term_vec)
      type(BLOCK_TYPE), intent(in)  :: blocks(:)
      real(dp), intent(in)          :: d_psi(:), d_omega(:)
      real(dp), intent(out)         :: adv_term_vec(:)

      integer   :: iblock, k, i_basis, j_basis
      integer   :: global_psi, global_omega
      integer   :: row_index, global_point_offset
      real(dp)  :: x, y, u_vel, v_vel, domega_dx, domega_dy 
      real(dp)  :: N_val, M_val, dN_dx, dM_dy 
      type(block_type) :: current_block

      adv_term_vec = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        current_block = blocks(iblock)
        do k = 1, current_block%N_x * current_block%N_y 
        
          if (current_block%boundary_types(k) /= BTYPE_INTERIOR) CYCLE 
          row_index = global_point_offset + k 
          x = current_block%colloc_pts(k, 1)
          y = current_block%colloc_pts(k, 2)
          u_vel = 0.0d0; v_vel = 0.0d0; domega_dx = 0.0d0; domega_dy = 0.0d0

          do j_basis = 1, current_block%N_y 
            do i_basis = 1, current_block%N_x
              global_psi    = get_global_psi_index(iblock, i_basis, j_basis)
              global_omega  = get_global_omega_index(iblock, i_basis, j_basis)

              N_val = bspline_basis_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
              M_val = bspline_basis_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)
              dN_dx = bspline_deriv1_physical(i_basis, current_block%P_x, current_block%knots_x, x, current_block%xmin, current_block%xmax)
              dM_dy = bspline_deriv1_physical(j_basis, current_block%P_y, current_block%knots_y, y, current_block%ymin, current_block%ymax)

              u_vel = u_vel + d_psi(global_psi) * (N_val * dM_dy)
              v_vel = v_vel - d_psi(global_psi) * (dN_dx * M_val)
              domega_dx = domega_dx + d_omega(global_omega) * (dN_dx * M_val)
              ! CORRECTED: Was N_val * M_val, now correctly uses dM_dy
              domega_dy = domega_dy + d_omega(global_omega) * (N_val * dM_dy)
            enddo 
          enddo 
          
          adv_term_vec(row_index) = -(u_vel * domega_dx + v_vel * domega_dy)
        enddo 
        global_point_offset = global_point_offset + (current_block%N_x * current_block%N_y)
      enddo 
    end subroutine calculate_advection_term

end module mod_solver