module hqt_path
use hqt_constants, only: hqt_dp
use hqt_type,      only: hqt_k_mesh
use hqt_tool,      only: hqt_range
implicit none

private

public :: hqt_path_map, hqt_closest_point

contains
    subroutine hqt_path_map(k_mesh, nodes, length, map, node_indexs)
        type(hqt_k_mesh), intent(in)  :: k_mesh
        real(hqt_dp),     intent(in)  :: nodes(:, :)
        integer,          intent(in)  :: length
        integer                       :: edge_lengths(size(nodes(1, :)) - 1)
        real(hqt_dp)                  :: edges(size(nodes(1, :)) - 1)
        real(hqt_dp), allocatable     :: x_range(:)
        real(hqt_dp), allocatable     :: y_range(:)
        integer                       :: point(2)
        integer                       :: i
        integer                       :: j
        integer                       :: k
        integer,          intent(out) :: map(length)
        integer,          intent(out) :: node_indexs(size(nodes(1, :)))
        do i = 1, size(edges)
            edges(i) = sqrt((nodes(1, i + 1) - nodes(1, i))**2 + (nodes(2, i + 1) - nodes(2, i))**2)
        end do
        do i = 1, size(edges)
            edge_lengths(i) = floor((length - 1) * edges(i) / sum(edges))
        end do
        map = -1
        k = 1
        allocate(x_range(maxval(edge_lengths)))
        allocate(y_range(maxval(edge_lengths)))
        do i = 1, size(edges)
            x_range = hqt_range(nodes(1, i), nodes(1, i + 1), edge_lengths(i) + 1)
            y_range = hqt_range(nodes(2, i), nodes(2, i + 1), edge_lengths(i) + 1)
            node_indexs(i) = k
            do j = 1, edge_lengths(i)
                point = hqt_closest_point(k_mesh, x_range(j), y_range(j))
                map(k) = k_mesh%map(point(1), point(2))
                k = k + 1
            end do
        end do
        node_indexs(size(node_indexs)) = k
        do i = 1, size(map)
            if (map(i) == -1) then
                point = hqt_closest_point(k_mesh, nodes(1, size(nodes(1, :))), nodes(2, size(nodes(1, :))))
                map(i) = k_mesh%map(point(1), point(2))
            end if
        end do
    
    end subroutine hqt_path_map
    
    function hqt_closest_point(k_mesh, kx, ky) result(point)
        type(hqt_k_mesh), intent(in) :: k_mesh
        real(hqt_dp),     intent(in) :: kx
        real(hqt_dp),     intent(in) :: ky
        integer                      :: point(2)
        point(1) = nint(((kx + maxval(k_mesh%reduced_points(1, :))) * k_mesh%length) / (2 * maxval(k_mesh%reduced_points(1, :))))
        point(2) = nint(((ky + maxval(k_mesh%reduced_points(2, :))) * k_mesh%length) / (2 * maxval(k_mesh%reduced_points(2, :))))
        if (point(1) < 1) point(1) = 1
        if (point(2) < 1) point(2) = 1
        if (point(1) > k_mesh%length) point(1) = k_mesh%length
        if (point(2) > k_mesh%length) point(2) = k_mesh%length
        
    end function hqt_closest_point

end module hqt_path
