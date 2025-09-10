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

  public  :: initialize_solver, get_global_psi_index, get_global_omega_index, assemble_constant_matrices, apply_psi_bcs
  public  :: calculate_advection_term, apply_omega_bcs, assemble_rk3_step
  public  :: extract_laplacian_operator 
  public  :: output_results

    
  contains 
    ! calculates and stores the number of coefficients and the starting index/offset for each block
    subroutine initialize_solver(blocks)
      type(block_type), intent(in) :: blocks(:)
      integer :: total_psi_coeffs

      num_coeffs_b1 = blocks(1)%N_x * blocks(1)%N_y
      num_coeffs_b2 = blocks(2)%N_x * blocks(2)%N_y

      total_psi_coeffs = num_coeffs_b1 + num_coeffs_b2 

      ! block1 psi coeffs starts at 0 offset, block2 start after block1 
      psi_offset_b1 = 0 
      psi_offset_b2 = num_coeffs_b1
      ! block1 omega coeffs start after all psi coeffs 
      omega_offset_b1 = 0
      omega_offset_b2 = num_coeffs_b1
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
    !   NOTE: bspline operates on normalized coord system. we need to scale it
    subroutine assemble_constant_matrices(blocks, A_psi, M_psi, M_omega)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(out)         :: A_psi(:,:), M_psi(:,:), M_omega(:,:)

      integer                       :: iblock, k, row_index, i_basis, j_basis 
      integer                       :: col_psi, col_omega 
      real(dp)                      :: x, y, u_norm, v_norm 
      real(dp)                      :: N_val, M_val, d2N_dx2_scaled, d2M_dy2_scaled
      real(dp)                      :: jac_x, jac_y, jac_x2, jac_y2
      integer                       :: global_point_offset

      A_psi = 0.0d0; M_psi = 0.0d0; M_omega = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        jac_x   = 1.0d0 / (blocks(iblock)%xmax - blocks(iblock)%xmin)
        jac_y   = 1.0d0 / (blocks(iblock)%ymax - blocks(iblock)%ymin)
        jac_x2  = jac_x * jac_x
        jac_y2  = jac_y * jac_y

        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
          row_index = global_point_offset + k 
          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)
          u_norm = (x - blocks(iblock)%xmin) * jac_x
          v_norm = (y - blocks(iblock)%ymin) * jac_y


          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
              col_omega = get_global_omega_index(iblock, i_basis, j_basis)
              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)

              M_psi(row_index, col_psi) = M_psi(row_index, col_psi) + N_val * M_val 
              M_omega(row_index, col_omega) = M_omega(row_index, col_omega) + N_val * M_val 
            enddo 
          enddo 

          ! ! for the A_psi, we only apply the nabla^2 operator at interior points. the boundary rows will be overwritten with BC 
          ! if (blocks(iblock)%boundary_types(k) == BTYPE_INTERIOR) then 
          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
              d2N_dx2_scaled   = bspline_deriv2(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm) * jac_x2
              d2M_dy2_scaled   = bspline_deriv2(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm) * jac_y2

              A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + d2N_dx2_scaled * M_val + N_val * d2M_dy2_scaled
            enddo 
          enddo 
        ! endif
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
      real(dp)                      :: x, y, u_norm, v_norm
      real(dp)                      :: N_val, M_val, dM_dy_scaled
      real(dp)                      :: u_lid
      real(dp)                      :: jac_x, jac_y
      real(dp)                      :: u1, v1, jac_y1
      real(dp)                      :: u2, v2, jac_y2
      

      b_psi_bc = 0.0d0; global_point_offset = 0

      ! loop through all blocks and all points 
      do iblock = 1, NUM_BLOCKS
        jac_x   = 1.0d0 / (blocks(iblock)%xmax - blocks(iblock)%xmin)
        jac_y   = 1.0d0 / (blocks(iblock)%ymax - blocks(iblock)%ymin)
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
          row_index = global_point_offset + k 
          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)
          u_norm = (x - blocks(iblock)%xmin) * jac_x
          v_norm = (y - blocks(iblock)%ymin) * jac_y

          select case (blocks(iblock)%boundary_types(k)) ! we only care about boundaries here
          case (BTYPE_INTERIOR)
            cycle 
          
          case (BTYPE_WALL) ! psi = 0 
            A_psi(row_index, :) = 0.0d0 
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
                M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
                A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * M_val
              enddo 
            enddo 
            b_psi_bc(row_index) = 0.0d0

          case (BTYPE_MOVING_LID) ! d(psi)/dy = u_lid(x,t)
            ! TODO: CHANGE THE VELOCITY PROFILE LATER FOR NOW I WILL USE SINUSOIDAL VELOCITY
            ! u_lid = sin(acos(-1.0d0) * x) * (0.5d0 * (1.0d0 - cos(acos(-1.0d0) * current_time)))
            u_lid = 1.0d0

            A_psi(row_index,:) = 0.0d0 
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_psi   = get_global_psi_index(iblock, i_basis, j_basis)
                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
                dM_dy_scaled     = bspline_deriv1(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm) * jac_y
                A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + N_val * dM_dy_scaled
              enddo 
            enddo 
            b_psi_bc(row_index) = u_lid

          case (BTYPE_INTERFACE) ! enforce continuity and derivative equivalance
            A_psi(row_index, :) = 0.0d0 

            ! NOTE: the coordinates (x,y) on the interface are physically same for both blocks
            ! we need to map (x,y) to the normalized coordinates of each block 
            ! 
            ! block 1 
            jac_y1  = 1.0d0 / (blocks(1)%ymax - blocks(1)%ymin)
            u1      = (x - blocks(1)%xmin) / (blocks(1)%xmax - blocks(1)%xmin)
            v1      = (y - blocks(1)%ymin) / (blocks(1)%ymax - blocks(1)%ymin)

            ! block 2 
            jac_y2  = 1.0d0 / (blocks(2)%ymax - blocks(2)%ymin)
            u2      = (x - blocks(2)%xmin) / (blocks(2)%xmax - blocks(2)%xmin)
            v2      = (y - blocks(2)%ymin) / (blocks(2)%ymax - blocks(2)%ymin)

            if (iblock == 1) then ! psi_1 - psi_2 = 0, coefficient multiplied by +1 
              do j_basis = 1, blocks(1)%N_y 
                do i_basis = 1, blocks(1)%N_x
                  col_psi   = get_global_psi_index(1, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(1)%P_x, blocks(1)%knots_x, u1)
                  M_val     = bspline_basis(j_basis, blocks(1)%P_y, blocks(1)%knots_y, v1)
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * M_val)
                enddo 
              enddo 

              ! coeffs from block2 multiplied by -1
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x
                  col_psi   = get_global_psi_index(2, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(2)%P_x, blocks(2)%knots_x, u2)
                  M_val     = bspline_basis(j_basis, blocks(2)%P_y, blocks(2)%knots_y, v2)
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
                  N_val     = bspline_basis(i_basis, blocks(1)%P_x, blocks(1)%knots_x, u1)
                  dM_dy_scaled     = bspline_deriv1(j_basis, blocks(1)%P_y, blocks(1)%knots_y, v1) * jac_y1
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) + (N_val * dM_dy_scaled) 
                enddo 
              enddo 

              ! negative contribution from block2
              do j_basis = 1, blocks(2)%N_y 
                do i_basis = 1, blocks(2)%N_x 
                  col_psi   = get_global_psi_index(2, i_basis, j_basis)
                  N_val     = bspline_basis(i_basis, blocks(2)%P_x, blocks(2)%knots_x, u2)
                  dM_dy_scaled     = bspline_deriv1(j_basis, blocks(2)%P_y, blocks(2)%knots_y, v2) * jac_y2
                  A_psi(row_index, col_psi) = A_psi(row_index, col_psi) - (N_val * dM_dy_scaled) 
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
      real(dp)  :: x, y, u_norm, v_norm, u_vel, v_vel, domega_dx, domega_dy 
      real(dp)  :: N_val, M_val, dN_dx_scaled, dM_dy_scaled 
      real(dp)  :: jac_x, jac_y

      adv_term_vec = 0.0d0; global_point_offset = 0

      do iblock = 1, NUM_BLOCKS
        jac_x   = 1.0d0 / (blocks(iblock)%xmax - blocks(iblock)%xmin)
        jac_y   = 1.0d0 / (blocks(iblock)%ymax - blocks(iblock)%ymin)
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
        
          ! we only need to calculate this at interior points 
          if (blocks(iblock)%boundary_types(k) /= BTYPE_INTERIOR) CYCLE 

          x = blocks(iblock)%colloc_pts(k, 1)
          y = blocks(iblock)%colloc_pts(k, 2)
          u_norm = (x - blocks(iblock)%xmin) * jac_x
          v_norm = (y - blocks(iblock)%ymin) * jac_y

          u_vel = 0.0d0; v_vel = 0.0d0; domega_dx = 0.0d0; domega_dy = 0.0d0

          do j_basis = 1, blocks(iblock)%N_y 
            do i_basis = 1, blocks(iblock)%N_x
              global_psi    = get_global_psi_index(iblock, i_basis, j_basis)
              global_omega  = get_global_omega_index(iblock, i_basis, j_basis)

              N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
              M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
              dN_dx_scaled     = bspline_deriv1(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm) * jac_x
              dM_dy_scaled     = bspline_deriv1(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm) * jac_y

              ! calculate velocities: u = d(psi)/dy, v = -d(psi)/dx 
              u_vel = u_vel + d_psi(global_psi) * (N_val * dM_dy_scaled)
              v_vel = v_vel - d_psi(global_psi) * (dN_dx_scaled * M_val)

              ! calculate vorticity gradients: d(omega)/dx, d(omega)/dy 
              domega_dx = domega_dx + d_omega(global_omega) * (dN_dx_scaled * M_val)
              domega_dy = domega_dy + d_omega(global_omega) * (N_val * dM_dy_scaled)
            enddo 
          enddo 

          ! calculate the advection term and store it 
          row_index = global_point_offset + k 
          adv_term_vec(row_index) = -(u_vel * domega_dx + v_vel * domega_dy)
        enddo 
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
      enddo 
    end subroutine calculate_advection_term


    subroutine apply_omega_bcs(blocks, d_psi_n, A_step, b_step)
      ! enforces the non-slip condition by setting the vorticity 0 at the wall. 
      ! LHS and RHS modifies for an implicit solve 
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(in)          :: d_psi_n(:)
      real(dp), intent(inout)       :: A_step(:,:), b_step(:)

      integer   :: iblock, k, row_index, i_basis, j_basis
      integer   :: col_psi, col_omega
      integer   :: global_point_offset 
      integer   :: i_local, j_local
      real(dp)  :: x, y, u_norm, v_norm, N_val, M_val, d2N_dx2_scaled, d2M_dy2_scaled
      real(dp)  :: omega_wall_val, laplacian_psi
      real(dp)  :: jac_x, jac_y, jac_x2, jac_y2

      global_point_offset = 0 

      do iblock = 1, NUM_BLOCKS
        jac_x   = 1.0d0 / (blocks(iblock)%xmax - blocks(iblock)%xmin)
        jac_y   = 1.0d0 / (blocks(iblock)%ymax - blocks(iblock)%ymin)
        jac_x2  = jac_x * jac_x
        jac_y2  = jac_y * jac_y
        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y 
          row_index = global_point_offset + k 

          ! we only want walls and lid points 
          select case (blocks(iblock)%boundary_types(k))
          case (BTYPE_WALL, BTYPE_MOVING_LID)
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)
            u_norm = (x - blocks(iblock)%xmin) * jac_x
            v_norm = (y - blocks(iblock)%ymin) * jac_y

            ! calculate required wall vorticity
            ! omega = -laplacian(psi)
            omega_wall_val = 0.0d0 

            ! omega_at_point = omega_wall_val
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_psi   = get_global_psi_index(iblock, i_basis, j_basis)

                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
                M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
                d2N_dx2_scaled   = bspline_deriv2(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm) * jac_x2
                d2M_dy2_scaled   = bspline_deriv2(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm) * jac_y2

                laplacian_psi = d_psi_n(col_psi) * (d2N_dx2_scaled * M_val + N_val * d2M_dy2_scaled)

                omega_wall_val = omega_wall_val - laplacian_psi
              enddo 
            enddo 

            A_step(row_index, :)  = 0.0d0 
            j_local = (k - 1) / blocks(iblock)%N_x + 1 
            i_local = mod(k - 1, blocks(iblock)%N_x) + 1
            A_step(row_index, get_global_omega_index(iblock, i_local, j_local)) = 1.0d0 ! dirichlet bc
            b_step(row_index)     = omega_wall_val  ! this is the RHS of the equation = omega_wall_val

            ! now we need to enforce it to the LHS too 
            do j_basis = 1, blocks(iblock)%N_y 
              do i_basis = 1, blocks(iblock)%N_x
                col_omega = get_global_omega_index(iblock, i_basis, j_basis)

                N_val     = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
                M_val     = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
                A_step(row_index, col_omega) = A_step(row_index, col_omega) + N_val * M_val ! this is the LHS = d_omega * N * M  

              enddo 
            enddo 
          end select 
        enddo 
        global_point_offset = global_point_offset + (blocks(iblock)%N_x * blocks(iblock)%N_y)
      enddo 
    end subroutine apply_omega_bcs 

    subroutine extract_laplacian_operator(blocks, A_psi, K_omega)
      type(block_type), intent(in)  :: blocks(:)
      real(dp), intent(in)          :: A_psi(:,:)
      real(dp), intent(out)         :: K_omega(:,:)

      ! for now we assume A_psi = k_omega 
      K_omega = A_psi 
    end subroutine extract_laplacian_operator

    subroutine assemble_rk3_step(blocks, step_num, A_step, b_step, dt, K, M, &
                                  d_omega_n, d_psi_n, d_omega_s1, d_omega_s2, &
                                  N_n, N_s1, N_s2)

      type(block_type), intent(in)  :: blocks(:)
      integer, intent(in)           :: step_num
      real(dp), intent(out)         :: A_step(:,:), b_step(:)
      real(dp), intent(in)          :: dt
      real(dp), intent(in)          :: K(:,:), M(:,:)
      real(dp), intent(in)          :: d_omega_n(:), d_psi_n(:), d_omega_s1(:), d_omega_s2(:), N_n(:), N_s1(:), N_s2(:)
      real(dp)                      :: diffusion_coeff = 1.0d0 / REYNOLDS_NUMBER 


      ! Select coefficients and RHS terms based on the step number 
      select case (step_num)
      case (1)
        ! {b1} = M*d_omega_n + dt * alpha1*(1/Re)*K*d_omega_n + dt*gamma1*N_n 
        b_step = matmul(M, d_omega_n) + (dt * rk_alpha(1) * diffusion_coeff) * matmul(K, d_omega_n) + dt * rk_gamma(1) * N_n 
        ! [A1] = [M] - dt * beta1 * (1/Re) * [K]
        A_step = M - (dt * rk_beta(1) * diffusion_coeff) * K 
      case (2)
        ! {b2} = M*d_omega_s1 + dt * alpha2*(1/Re)*K*d_omega_s1 + dt*gamma2*N_s1 + dt*zeta2*N_n 
        b_step = matmul(M, d_omega_s1) + (dt * rk_alpha(2) * diffusion_coeff) * matmul(K, d_omega_s1) + & 
                dt * rk_gamma(2) * N_s1 + dt * rk_zeta(2) * N_n 
        ! [A2] = [M] - dt * beta2 * (1/Re) * [K]
        A_step = M - (dt * rk_beta(2) * diffusion_coeff) * K 
      case (3)
        ! {b3} = M*d_omega_s2 + dt * alpha3*(1/Re)*K*d_omega_s2 + dt*gamma3*N_s2 + dt*zeta3*N_s1 
        b_step = matmul(M, d_omega_s2) + (dt * rk_alpha(3) * diffusion_coeff) * matmul(K, d_omega_s2) + & 
                dt * rk_gamma(3) * N_s2 + dt * rk_zeta(3) * N_s1
        ! [A3] = [M] - dt * beta3 * (1/Re) * [K]
        A_step = M - (dt * rk_beta(3) * diffusion_coeff) * K 
      end select 

      ! assemble the implicit LHS matrix 

      ! overwrite boundary rows with wall vorticity condition 
      call apply_omega_bcs(blocks, d_psi_n, A_step, b_step)


    end subroutine assemble_rk3_step

    ! subroutine assemble_rk3_step1(blocks, A_step, b_step, dt, K, M, d_omega_n, d_psi_n, N_n)
    !   type(block_type), intent(in)  :: blocks(:)
    !   real(dp), intent(out)         :: A_step(:,:), b_step(:)
    !   real(dp), intent(in)          :: dt
    !   real(dp), intent(in)          :: K(:,:), M(:,:)
    !   real(dp), intent(in)          :: d_omega_n(:), d_psi_n(:), N_n(:)
    !
    !   real(dp), parameter :: alpha1 = 29.0d0 / 60.0d0 
    !   real(dp), parameter :: beta1  = 37.0d0 / 160.0d0
    !   real(dp), parameter :: gamma1 = 8.0d0 / 15.0d0 
    !
    !   real(dp)            :: diffusion_coeff = 1.0d0 / REYNOLDS_NUMBER 
    !
    !   ! assemble the implicit LHS matrix for interior points 
    !   ! [A1] = [M] - dt * beta1 * (1/Re) * [K]
    !   A_step = M - (dt * beta1 * diffusion_coeff) * K 
    !
    !   ! assemble the explicit RHS vector 
    !   ! {b1} = M*d_omega_n + dt * alpha1*(1/Re)*K*d_omega_n + dt*gamma1*N_n 
    !   b_step = matmul(M, d_omega_n) + (dt * alpha1 * diffusion_coeff) * matmul(K, d_omega_n) + dt * gamma1 * N_n 
    !
    !   ! overwrite boundary rows with wall vorticity condition 
    !   call apply_omega_bcs(blocks, d_psi_n, A_step, b_step)
    !
    ! end subroutine assemble_rk3_step1
    !
    ! subroutine assemble_rk3_step2(blocks, A_step, b_step, dt, K, M, d_omega_s1, d_psi_n, N_n, N_s1)
    !   type(block_type), intent(in)  :: blocks(:)
    !   real(dp), intent(out)         :: A_step(:,:), b_step(:)
    !   real(dp), intent(in)          :: dt
    !   real(dp), intent(in)          :: K(:,:), M(:,:)
    !   real(dp), intent(in)          :: d_omega_s1(:), d_psi_n(:), N_n(:), N_s1(:)
    !
    !   real(dp), parameter :: alpha2 = -3.0d0 / 40.0d0
    !   real(dp), parameter :: beta2  = 5.0d0 / 24.0d0 
    !   real(dp), parameter :: gamma2 = 5.0d0 / 12.0d0 
    !   real(dp), parameter :: zeta2  = -17.0d0 / 60.0d0
    !
    !   real(dp)            :: diffusion_coeff = 1.0d0 / REYNOLDS_NUMBER 
    !
    !   ! LHS vector 
    !   ! [A2] = [M] - dt * beta2 * (1/Re) * [K]
    !   A_step = M - (dt * beta2 * diffusion_coeff) * K 
    !
    !   ! RHS vector 
    !   ! {b2} = M*d_omega_s1 + dt * alpha2*(1/Re)*K*d_omega_s1 + dt*gamma2*N_s1 + dt*zeta2*N_n 
    !   b_step = matmul(M, d_omega_s1) + (dt * alpha2 * diffusion_coeff) * matmul(K, d_omega_s1) + & 
    !            dt * gamma2 * N_s1 + dt * zeta2 * N_n 
    !   ! overwrite boundary rows with wall vorticity condition 
    !   call apply_omega_bcs(blocks, d_psi_n, A_step, b_step)
    !
    ! end subroutine assemble_rk3_step2
    !
    ! subroutine assemble_rk3_step3(blocks, A_step, b_step, dt, K, M, d_omega_s2, d_psi_n, N_s1, N_s2)
    !   type(block_type), intent(in)  :: blocks(:)
    !   real(dp), intent(out)         :: A_step(:,:), b_step(:)
    !   real(dp), intent(in)          :: dt
    !   real(dp), intent(in)          :: K(:,:), M(:,:)
    !   real(dp), intent(in)          :: d_omega_s2(:), d_psi_n(:), N_s1(:), N_s2(:)
    !
    !   real(dp), parameter :: alpha3 = 1.0d0 / 6.0d0
    !   real(dp), parameter :: beta3  = 1.0d0 / 6.0d0
    !   real(dp), parameter :: gamma3 = 3.0d0 / 4.0d0 
    !   real(dp), parameter :: zeta3  = -5.0d0 / 12.0d0
    !
    !   real(dp)            :: diffusion_coeff = 1.0d0 / REYNOLDS_NUMBER 
    !
    !   ! LHS vector 
    !   ! [A3] = [M] - dt * beta3 * (1/Re) * [K]
    !   A_step = M - (dt * beta3 * diffusion_coeff) * K 
    !
    !   ! RHS vector 
    !   ! {b3} = M*d_omega_s2 + dt * alpha3*(1/Re)*K*d_omega_s2 + dt*gamma3*N_s2 + dt*zeta3*N_s1 
    !   b_step = matmul(M, d_omega_s2) + (dt * alpha3 * diffusion_coeff) * matmul(K, d_omega_s2) + & 
    !            dt * gamma3 * N_s2 + dt * zeta3 * N_s1
    !   ! overwrite boundary rows with wall vorticity condition 
    !   call apply_omega_bcs(blocks, d_psi_n, A_step, b_step)
    !
    ! end subroutine assemble_rk3_step3

      
	  ! Add this subroutine to the end of the mod_solver module, replacing the old one.

subroutine output_results(blocks, d_psi, d_omega, timestep, file_prefix)
    type(block_type), intent(in) :: blocks(:)
    real(dp), intent(in)         :: d_psi(:), d_omega(:)
    integer, intent(in)          :: timestep
    character(len=*), intent(in) :: file_prefix

    character(len=256) :: filename
    integer, parameter :: out_unit = 20
    integer :: iblock, k, i_basis, j_basis
    real(dp) :: x, y, u_norm, v_norm, psi_val, omega_val, u_vel, v_vel
    real(dp) :: N_val, M_val, dN_dx, dM_dy
    real(dp) :: jac_x, jac_y
    integer :: g_psi_idx, g_omega_idx

    write(filename, '(A, I0.6, A)') trim(file_prefix), timestep, '.dat'
    open(unit=out_unit, file=filename, status='replace', action='write')

    write(out_unit, '(A)') '# Columns: X, Y, Streamfunction(psi), Vorticity(omega), U_Velocity, V_Velocity'

    do iblock = 1, NUM_BLOCKS
        ! --- DEFINE SCALING FACTORS (JACOBIANS) FOR THIS BLOCK ---
        jac_x = 1.0_dp / (blocks(iblock)%xmax - blocks(iblock)%xmin)
        jac_y = 1.0_dp / (blocks(iblock)%ymax - blocks(iblock)%ymin)

        do k = 1, blocks(iblock)%N_x * blocks(iblock)%N_y
            x = blocks(iblock)%colloc_pts(k, 1)
            y = blocks(iblock)%colloc_pts(k, 2)

            ! --- MAP PHYSICAL COORDS (x,y) TO NORMALIZED COORDS (u_norm, v_norm) ---
            u_norm = (x - blocks(iblock)%xmin) * jac_x
            v_norm = (y - blocks(iblock)%ymin) * jac_y

            psi_val = 0.0_dp; omega_val = 0.0_dp; u_vel = 0.0_dp; v_vel = 0.0_dp

            do j_basis = 1, blocks(iblock)%N_y
                do i_basis = 1, blocks(iblock)%N_x
                    ! --- EVALUATE USING NORMALIZED COORDS AND SCALE DERIVATIVES ---
                    N_val = bspline_basis(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm)
                    dN_dx = bspline_deriv1(i_basis, blocks(iblock)%P_x, blocks(iblock)%knots_x, u_norm) * jac_x
                    M_val = bspline_basis(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm)
                    dM_dy = bspline_deriv1(j_basis, blocks(iblock)%P_y, blocks(iblock)%knots_y, v_norm) * jac_y

                    g_psi_idx   = get_global_psi_index(iblock, i_basis, j_basis)
                    g_omega_idx = get_global_omega_index(iblock, i_basis, j_basis)

                    psi_val   = psi_val   + d_psi(g_psi_idx) * N_val * M_val
                    omega_val = omega_val + d_omega(g_omega_idx) * N_val * M_val
                    u_vel     = u_vel     + d_psi(g_psi_idx) * N_val * dM_dy
                    v_vel     = v_vel     - d_psi(g_psi_idx) * dN_dx * M_val
                enddo
            enddo
            write(out_unit, '(6E22.12)') x, y, psi_val, omega_val, u_vel, v_vel
        enddo
    enddo
    close(out_unit)
end subroutine output_results


end module mod_solver
