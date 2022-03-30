module route
    use constants,  only: dp, pi
    use parameters, only: num_k_length, length_scale
    use tool,       only: tool_count_unique, tool_location
    use type,       only: type_crystal
    use spglib_f08
    implicit none
        
    private
    
    public  :: route_reduced_k_mesh

contains
    subroutine route_reduced_k_mesh(Kxs, Kys, weight, k_path_map, crystal)
        type(type_crystal(*, *, *))                 :: crystal
        real(dp), intent(out), allocatable :: kxs(:), kys(:)
        integer, intent(out), allocatable  :: weight(:)
        integer, intent(out), allocatable  :: k_path_map(:)
        real(dp)                           :: positions_2d(size(crystal%atom_positions(:, 1)), size(crystal%atom_positions(1, :)))
        integer                            :: map(floor(num_k_length / length_scale)**2), &
                                              grid_point(3, &
                                                         floor(num_k_length &
                                                           / length_scale)**2), & 
                                              num_ir_grid, i, j, n
        integer, allocatable               :: used_map_to_map(:), used_grid_point(:, :), used_map(:)
        positions_2d = crystal%atom_positions
        positions_2d(3, :) = 0.0_dp
        num_ir_grid = spg_get_ir_reciprocal_mesh(grid_point, map, &
                      (/ floor(num_k_length / length_scale), &
                         floor(num_k_length / length_scale), &
                         1 /), &
                      (/ 0, 0, 0 /), 1, &
                      crystal%lattice_vectors, positions_2d, crystal%atom_types, size(crystal%atom_types), 1e-5_dp)
        allocate(used_map(num_k_length**2))
        used_map = -1
        do i = 1, num_k_length
            do j = 1, num_k_length
                if (i <= num_k_length / 2 &
                    .and. j >= num_k_length / 2 + mod(num_k_length, 2)) then
                    used_map((i - 1) * num_k_length + j) = map(&
                        ((i + floor(dble(num_k_length) / length_scale) &
                            - num_k_length / 2) - 1) &
                          * floor(dble(num_k_length / length_scale)) &
                        + (j - num_k_length / 2 - mod(num_k_length, 2))) + 1
                end if
                if (i >= num_k_length / 2 + mod(num_k_length, 2) &
                    .and. j <= num_k_length / 2) then
                    used_map((i - 1) * num_k_length + j) = map(&
                        ((i - num_k_length / 2) - mod(num_k_length, 2)) &
                            * floor(dble(num_k_length / length_scale)) &
                        + (j + floor(dble(num_k_length) / length_scale) &
                             - num_k_length / 2)) + 1
                end if
                if (i <= num_k_length / 2 .and. j <= num_k_length / 2) then
                    used_map((i - 1) * num_k_length + j) = map(&
                        ((i + floor(dble(num_k_length) / length_scale) &
                            - num_k_length / 2) - 1) &
                          * floor(dble(num_k_length / length_scale)) &
                        + (j + floor(dble(num_k_length) / length_scale) &
                             - num_k_length / 2)) + 1
                end if
                if (i > num_k_length / 2 .and. j > num_k_length / 2) then
                    used_map((i - 1) * num_k_length + j) = map(&
                        ((i - num_k_length / 2) - 1) &
                                * floor(num_k_length / length_scale) &
                            + (j - num_k_length / 2)) + 1
                end if
            end do
        end do
        allocate(used_map_to_map(tool_count_unique(used_map)))
        allocate(weight(size(used_map_to_map)))
        allocate(used_grid_point(3, size(used_map_to_map)))
        allocate(kxs(size(used_map_to_map)))
        allocate(kys(size(used_map_to_map)))
        used_map_to_map = -1
        weight = 0
        n = 0
        do i = 1, size(used_map)
            if (used_map(i) == -1) then
                cycle
            end if
            if (any(used_map_to_map == used_map(i))) then
                used_map(i) = tool_location(used_map_to_map, used_map(i))
            else ! Unique
                n = n + 1
                used_map_to_map(n)    = used_map(i)
                used_grid_point(:, n) = grid_point(:, used_map(i))
                used_map(i) = n
            end if
        end do
        do i = 1, size(used_map)
            if (used_map(i) == -1) then
                cycle
            end if
            weight(used_map(i)) = weight(used_map(i)) + 1
        end do
        do i = 1, size(used_grid_point(1, :))
            kxs(i) = 2 * pi * dble(used_grid_point(1, i)) &
                     / floor(num_k_length / length_scale)
            kys(i) = 2 * pi * dble(used_grid_point(2, i)) &
                     / floor(num_k_length / length_scale)
        end do
        allocate(k_path_map(2 * floor((num_k_length - 1) / 2.0_dp) + 1))
        do i = 1, floor((num_k_length - 1) / 2.0_dp)
            k_path_map(i) = used_map((num_k_length - i) * num_k_length &
                                + num_k_length - i + 1)
        end do
        do i = floor((num_k_length - 1) / 2.0_dp) + 2, size(k_path_map)
            k_path_map(i) = used_map(&
                (i + mod(num_k_length, &
                 floor((num_k_length - 1) / 2.0_dp)) - 2) * num_k_length &
                + num_k_length - floor((num_k_length - 1) / 2.0_dp))
        end do
        k_path_map(floor((num_k_length - 1) / 2.0_dp) + 1) = used_map(&
            (num_k_length - floor((num_k_length - 1) / 2.0_dp) - 1) &
                * num_k_length &
            + num_k_length - floor((num_k_length - 1) / 2.0_dp))
        
    end subroutine route_reduced_k_mesh

end module route
