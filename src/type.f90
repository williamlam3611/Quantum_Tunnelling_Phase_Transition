module hqt_type
use hqt_spglib_f08, only: spg_get_ir_reciprocal_mesh
use hqt_constants,  only: hqt_dp, hqt_pi
use hqt_tool,       only: operator(+), hqt_token_write, hqt_token_read, hqt_index, &
                          hqt_to_string, hqt_to_integer, hqt_to_complex, hqt_to_real
implicit none

private

public :: hqt_crystal, hqt_heterostructure, hqt_k_mesh, hqt_serialise, &
          hqt_deserialise_crystal, hqt_deserialise_heterostructure, hqt_deserialise_k_mesh


type hqt_crystal
    real(hqt_dp)                 :: crystal_length
    real(hqt_dp)                 :: lattice_vectors(3, 3)
    real(hqt_dp), allocatable    :: atom_positions(:, :)
    integer,  allocatable    :: atom_types(:)
    
    ! Tunneling Potential
    complex(hqt_dp), allocatable :: tunneling(:, :, :, :, :)
    integer, allocatable     :: tunneling_weight(:, :, :)
    
    ! Band Structure
    integer                  :: num_bands
    integer                  :: max_hopping
    
    ! Conducting band minimum
    real(hqt_dp)                 :: conducting_band_minimum_kx = 0.0_hqt_dp
    real(hqt_dp)                 :: conducting_band_minimum_ky = 0.0_hqt_dp

end type

type hqt_heterostructure
    type(hqt_crystal)    :: material
    integer               :: num_layers
    real(hqt_dp), allocatable :: potential(:)
    
    integer               :: min_layer
    integer               :: max_layer
    integer               :: min_state
    integer               :: max_state

end type

type hqt_k_mesh
    integer               :: length
    real(hqt_dp)              :: scale
    integer               :: num_reduced_points
    real(hqt_dp), allocatable :: reduced_points(:, :)
    integer,  allocatable :: weight(:)
    integer,  allocatable :: map(:, :)
    
end type

interface hqt_crystal
    module procedure hqt_crystal_arguments
end interface hqt_crystal

interface hqt_k_mesh
    module procedure hqt_k_mesh_crystal
end interface hqt_k_mesh

interface hqt_heterostructure
    module procedure hqt_heterostructure_arguments
end interface hqt_heterostructure

interface hqt_serialise
    module procedure hqt_serialise_crystal, hqt_serialise_heterostructure, hqt_serialise_k_mesh
end interface hqt_serialise


contains
    function hqt_crystal_arguments(crystal_length, lattice_vectors, atom_positions, atom_types, &
                                    tunneling, tunneling_weight, num_bands, max_hopping, &
                                    conducting_band_minimum_kx, conducting_band_minimum_ky) result (material)
        real(hqt_dp), intent(in) :: crystal_length
        real(hqt_dp), intent(in) :: lattice_vectors(3, 3)
        real(hqt_dp), intent(in) :: atom_positions(:, :)
        integer,  intent(in) :: atom_types(:)
        real(hqt_dp), intent(in) :: tunneling(:, :, :, :, :)
        integer,  intent(in) :: tunneling_weight(:, :, :)
        integer,  intent(in) :: num_bands
        integer,  intent(in) :: max_hopping
        real(hqt_dp), intent(in), optional :: conducting_band_minimum_kx
        real(hqt_dp), intent(in), optional :: conducting_band_minimum_ky
        type(hqt_crystal)   :: material
        material%crystal_length = crystal_length
        material%lattice_vectors = lattice_vectors
        material%atom_positions = atom_positions
        material%atom_types = atom_types
        material%tunneling = tunneling
        material%tunneling_weight = tunneling_weight
        material%num_bands = num_bands
        material%max_hopping = max_hopping
        if (present(conducting_band_minimum_kx)) material%conducting_band_minimum_kx = conducting_band_minimum_kx
        if (present(conducting_band_minimum_ky)) material%conducting_band_minimum_ky = conducting_band_minimum_ky
        
    end function hqt_crystal_arguments
    
    function hqt_heterostructure_arguments(crystal, num_layers) result(structure)
        type(hqt_crystal)         :: crystal
        integer               :: num_layers
        type(hqt_heterostructure) :: structure
        allocate(structure%potential(num_layers))
        structure%potential  = 0.0_hqt_dp
        structure%num_layers = num_layers
        structure%material   = crystal
        structure%min_layer  = 1
        structure%max_layer  = structure%num_layers
        structure%min_state  = 1
        structure%max_state  = structure%material%num_bands * structure%num_layers
    
    end function hqt_heterostructure_arguments
    
    function hqt_k_mesh_crystal(length, crystal, scale) result(mesh)
        integer,            intent(in) :: length
        type(hqt_crystal), intent(in) :: crystal
        real(hqt_dp),           intent(in) :: scale
        real(hqt_dp)                       :: positions(size(crystal%atom_positions(:, 1)), size(crystal%atom_positions(1, :)))
        integer                        :: grid_point(3, floor(length / scale)**2)
        integer                        :: map(floor(length / scale)**2)
        integer                        :: i
        integer                        :: j
        integer, allocatable           :: set(:)
        type(hqt_k_mesh)              :: mesh
        mesh%length = length
        mesh%scale = scale
        allocate(mesh%map(length, length))
        
        positions       = crystal%atom_positions
        positions(3, :) = 0.0_hqt_dp
        mesh%num_reduced_points = spg_get_ir_reciprocal_mesh(grid_point, map, &
                                                        (/ floor(length / scale), floor(length / scale), 1 /), &
                                                        (/ 0, 0, 0 /), 1, crystal%lattice_vectors, positions, crystal%atom_types, &
                                                        size(crystal%atom_types), 1e-5_hqt_dp)
        allocate(set(mesh%num_reduced_points))
        allocate(mesh%reduced_points(2, mesh%num_reduced_points))
        allocate(mesh%weight(mesh%num_reduced_points))
        
        set = -1
        do i = 1, size(map)
            if (.not. any(map(i) == set)) then
                do j = 1, size(set)
                    if (set(j) == -1) then
                        set(j) = map(i)
                        exit
                    end if
                end do
            end if
        end do
                
        mesh%weight = 0              
        do i = 1, length
            do j = 1, length
                if (i <= length / 2 .and. j >= length / 2 + mod(length, 2)) then
                    mesh%map(i, j) = hqt_index(set, map(((i + floor(length / scale) - length / 2) - 1) * floor(length / scale) &
                                         + (j - length / 2 - mod(length, 2))))
                else if (i >= length / 2 + mod(length, 2) .and. j <= length / 2) then
                    mesh%map(i, j) = hqt_index(set, map(((i - length / 2) - mod(length, 2)) * floor(length / scale) &
                                         + (j + floor(length / scale) - length / 2)))
                else if (i <= length / 2 .and. j <= length / 2) then
                    mesh%map(i, j) = hqt_index(set, map(((i + floor(length / scale) - length / 2) - 1) * floor(length / scale) &
                                         + (j + floor(length / scale) - length / 2)))
                else if (i > length / 2 .and. j > length / 2) then
                    mesh%map(i, j) = hqt_index(set, map(((i - length / 2) - 1) * floor(length / scale) &
                                         + (j - length / 2)))
                end if
                mesh%weight(mesh%map(i, j)) = mesh%weight(mesh%map(i, j)) + 1
            end do
        end do
        
        do i = 1, mesh%num_reduced_points
            mesh%reduced_points(1, i) = 2 * hqt_pi * dble(grid_point(1, set(i) + 1)) / floor(length / scale)
            mesh%reduced_points(2, i) = 2 * hqt_pi * dble(grid_point(2, set(i) + 1)) / floor(length / scale)
        end do
    
    end function hqt_k_mesh_crystal
    
    
    function hqt_serialise_crystal(crystal) result(string)
        type(hqt_crystal), intent(in) :: crystal
        character(:), allocatable      :: string
        string = string + hqt_token_write("crystal_length",          hqt_to_string(crystal%crystal_length)) + "\n"
        string = string + hqt_token_write("num_bands",               hqt_to_string(crystal%num_bands)) + "\n"
        string = string + hqt_token_write("max_hopping",             hqt_to_string(crystal%max_hopping)) + "\n"
        string = string + hqt_token_write("conducting_band_minimum_kx", hqt_to_string(crystal%conducting_band_minimum_kx)) + "\n"
        string = string + hqt_token_write("conducting_band_minimum_ky", hqt_to_string(crystal%conducting_band_minimum_ky)) + "\n"
        
        string = string + hqt_token_write("lattice_vectors",         hqt_to_string(pack(crystal%lattice_vectors, .true.))) + "\n"
        string = string + hqt_token_write("atom_types",              hqt_to_string(     crystal%atom_types)) + "\n"
        string = string + hqt_token_write("atom_positions",          hqt_to_string(pack(crystal%atom_positions, .true.))) + "\n"
        string = string + hqt_token_write("tunneling_weight",        hqt_to_string(pack(crystal%tunneling_weight, .true.))) + "\n"
        string = string + hqt_token_write("tunneling",               hqt_to_string(pack(crystal%tunneling, .true.)))
        
    end function hqt_serialise_crystal
    
    function hqt_deserialise_crystal(string) result(crystal)
        character(*), intent(in) :: string
        integer, allocatable     :: temp_integer(:)
        type(hqt_crystal)       :: crystal
        crystal%crystal_length             = minval(hqt_to_real(   hqt_token_read(string, "crystal_length")), dim = 1)
        crystal%num_bands                  = minval(hqt_to_integer(hqt_token_read(string, "num_bands")), dim = 1)
        crystal%max_hopping                = minval(hqt_to_integer(hqt_token_read(string, "max_hopping")), dim = 1)
        crystal%conducting_band_minimum_kx = minval(hqt_to_real(   hqt_token_read(string, "conducting_band_minimum_kx")), dim = 1)
        crystal%conducting_band_minimum_ky = minval(hqt_to_real(   hqt_token_read(string, "conducting_band_minimum_ky")), dim = 1)
        
        crystal%lattice_vectors            = reshape(hqt_to_real(   hqt_token_read(string, "lattice_vectors")), (/ 3, 3 /))
        temp_integer                       =         hqt_to_integer(hqt_token_read(string, "atom_types"))
        allocate(crystal%atom_types(      size(temp_integer)))
        crystal%atom_types                 = temp_integer 
        allocate(crystal%atom_positions(  3, size(crystal%atom_types)))            
        crystal%atom_positions             = reshape(hqt_to_real(   hqt_token_read(string, "atom_positions")), &
                                                     (/ 3, size(crystal%atom_types) /))
        allocate(crystal%tunneling_weight(2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1))
        crystal%tunneling_weight           = reshape(hqt_to_integer(hqt_token_read(string, "tunneling_weight")), &
                                                     (/ 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                                        2 * crystal%max_hopping + 1 /))
        allocate(crystal%tunneling(       2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                          crystal%num_bands, crystal%num_bands))                          
        crystal%tunneling                  = reshape(hqt_to_complex(hqt_token_read(string, "tunneling")), &
                                                     (/ 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                                        2 * crystal%max_hopping + 1, crystal%num_bands, crystal%num_bands /))
                                                        
    end function hqt_deserialise_crystal
    
    function hqt_serialise_heterostructure(structure) result(string)
        type(hqt_heterostructure), intent(in) :: structure
        character(:), allocatable              :: string
        string = string + hqt_token_write("num_layers", hqt_to_string(structure%num_layers)) + "\n"
        string = string + hqt_token_write("min_layer",  hqt_to_string(structure%min_layer)) + "\n"
        string = string + hqt_token_write("max_layer",  hqt_to_string(structure%max_layer)) + "\n"
        string = string + hqt_token_write("min_state",  hqt_to_string(structure%min_state)) + "\n"
        string = string + hqt_token_write("max_state",  hqt_to_string(structure%max_state)) + "\n"
        
        string = string + hqt_token_write("potential",  hqt_to_string(structure%potential)) + "\n"
        string = string + hqt_token_write("material",   "{" + hqt_serialise(structure%material) + "}")
    
    end function hqt_serialise_heterostructure
    
    function hqt_deserialise_heterostructure(string) result(structure)
        character(*), intent(in)   :: string
        type(hqt_heterostructure) :: structure
        structure%num_layers = minval(hqt_to_real(hqt_token_read(string, "num_layers")), dim = 1)
        structure%min_layer = minval(hqt_to_real(hqt_token_read(string, "min_layer")), dim = 1)
        structure%max_layer = minval(hqt_to_real(hqt_token_read(string, "max_layer")), dim = 1)
        structure%min_state = minval(hqt_to_real(hqt_token_read(string, "min_state")), dim = 1)
        structure%max_state = minval(hqt_to_real(hqt_token_read(string, "max_state")), dim = 1)
        
        allocate(structure%potential(structure%num_layers))
        structure%potential = hqt_to_real(hqt_token_read(string, "potential"))
        structure%material  = hqt_deserialise_crystal(hqt_token_read(string, "material"))
                                                        
    end function hqt_deserialise_heterostructure
    
    
    function hqt_serialise_k_mesh(mesh) result(string)
        type(hqt_k_mesh), intent(in) :: mesh
        character(:), allocatable     :: string
        string = string + hqt_token_write("length",             hqt_to_string(mesh%length)) + "\n"
        string = string + hqt_token_write("scale",              hqt_to_string(mesh%scale)) + "\n"
        string = string + hqt_token_write("num_reduced_points", hqt_to_string(mesh%num_reduced_points)) + "\n"
        
        string = string + hqt_token_write("reduced_points",     hqt_to_string(pack(mesh%reduced_points, .true.))) + "\n"
        string = string + hqt_token_write("weight",             hqt_to_string(mesh%weight)) + "\n"
        string = string + hqt_token_write("map",                hqt_to_string(pack(mesh%map, .true.)))
    
    end function hqt_serialise_k_mesh
    
    function hqt_deserialise_k_mesh(string) result(mesh)
        character(*), intent(in)   :: string
        type(hqt_k_mesh)          :: mesh
        mesh%length             = minval(hqt_to_integer(hqt_token_read(string, "length")), dim = 1)
        mesh%scale              = minval(hqt_to_real(   hqt_token_read(string, "scale")), dim = 1)
        mesh%num_reduced_points = minval(hqt_to_integer(hqt_token_read(string, "num_reduced_points")), dim = 1)
        
        allocate(mesh%reduced_points(2, mesh%num_reduced_points))
        allocate(mesh%weight(mesh%num_reduced_points))
        allocate(mesh%map(mesh%length, mesh%length))
        
        mesh%reduced_points = reshape(hqt_to_real(hqt_token_read(   string, "reduced_points")), shape(mesh%reduced_points))
        mesh%weight         = reshape(hqt_to_integer(hqt_token_read(string, "weight")), shape(mesh%weight))
        mesh%map            = reshape(hqt_to_integer(hqt_token_read(string, "map")), shape(mesh%map))
                                                        
    end function hqt_deserialise_k_mesh
    

end module hqt_type
