module type
use spglib_f08, only: spg_get_ir_reciprocal_mesh
use constants,  only: dp, pi
use tool,       only: operator(+), tool_token_write, tool_token_read, tool_index, &
                      tool_to_string, tool_to_integer, tool_to_complex, tool_to_real
implicit none

private

public :: type_crystal, type_heterostructure, type_k_mesh, type_serialise, &
          type_deserialise_crystal, type_deserialise_heterostructure, type_deserialise_k_mesh


type type_crystal
    real(dp)                 :: crystal_length
    real(dp)                 :: lattice_vectors(3, 3)
    real(dp), allocatable    :: atom_positions(:, :)
    integer,  allocatable    :: atom_types(:)
    
    ! Tunneling Potential
    complex(dp), allocatable :: tunneling(:, :, :, :, :)
    integer, allocatable     :: tunneling_weight(:, :, :)
    
    ! Band Structure
    integer                  :: num_bands
    integer                  :: max_hopping
    
    ! Conducting band minimum
    real(dp)                 :: conducting_band_minimum_kx = 0.0_dp
    real(dp)                 :: conducting_band_minimum_ky = 0.0_dp

end type

type type_heterostructure
    type(type_crystal)    :: material
    integer               :: num_layers
    real(dp), allocatable :: potential(:)
    
    integer               :: min_layer
    integer               :: max_layer
    integer               :: min_state
    integer               :: max_state

end type

type type_k_mesh
    integer               :: length
    real(dp)              :: scale
    integer               :: num_reduced_points
    real(dp), allocatable :: reduced_points(:, :)
    integer,  allocatable :: weight(:)
    integer,  allocatable :: map(:, :)
    
end type

interface type_crystal
    module procedure type_crystal_arguments
end interface type_crystal

interface type_k_mesh
    module procedure type_k_mesh_crystal
end interface type_k_mesh

interface type_heterostructure
    module procedure type_heterostructure_arguments
end interface type_heterostructure

interface type_serialise
    module procedure type_serialise_crystal, type_serialise_heterostructure, type_serialise_k_mesh
end interface type_serialise


contains
    function type_crystal_arguments(crystal_length, lattice_vectors, atom_positions, atom_types, &
                                    tunneling, tunneling_weight, num_bands, max_hopping, &
                                    conducting_band_minimum_kx, conducting_band_minimum_ky) result (material)
        real(dp), intent(in) :: crystal_length
        real(dp), intent(in) :: lattice_vectors(3, 3)
        real(dp), intent(in) :: atom_positions(:, :)
        integer,  intent(in) :: atom_types(:)
        real(dp), intent(in) :: tunneling(:, :, :, :, :)
        integer,  intent(in) :: tunneling_weight(:, :, :)
        integer,  intent(in) :: num_bands
        integer,  intent(in) :: max_hopping
        real(dp), intent(in), optional :: conducting_band_minimum_kx
        real(dp), intent(in), optional :: conducting_band_minimum_ky
        type(type_crystal)   :: material
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
        
    end function type_crystal_arguments
    
    function type_heterostructure_arguments(crystal, num_layers) result(structure)
        type(type_crystal)         :: crystal
        integer               :: num_layers
        type(type_heterostructure) :: structure
        allocate(structure%potential(num_layers))
        structure%potential  = 0.0_dp
        structure%num_layers = num_layers
        structure%material   = crystal
        structure%min_layer  = 1
        structure%max_layer  = structure%num_layers
        structure%min_state  = 1
        structure%max_state  = structure%material%num_bands * structure%num_layers
    
    end function type_heterostructure_arguments
    
    function type_k_mesh_crystal(length, crystal, scale) result(mesh)
        integer,            intent(in) :: length
        type(type_crystal), intent(in) :: crystal
        real(dp),           intent(in) :: scale
        real(dp)                       :: positions(size(crystal%atom_positions(:, 1)), size(crystal%atom_positions(1, :)))
        integer                        :: grid_point(3, floor(length / scale)**2)
        integer                        :: map(floor(length / scale)**2)
        integer                        :: i
        integer                        :: j
        integer, allocatable           :: set(:)
        type(type_k_mesh)              :: mesh
        mesh%length = length
        mesh%scale = scale
        allocate(mesh%map(length, length))
        
        positions       = crystal%atom_positions
        positions(3, :) = 0.0_dp
        mesh%num_reduced_points = spg_get_ir_reciprocal_mesh(grid_point, map, &
                                                        (/ floor(length / scale), floor(length / scale), 1 /), &
                                                        (/ 0, 0, 0 /), 1, crystal%lattice_vectors, positions, crystal%atom_types, &
                                                        size(crystal%atom_types), 1e-5_dp)
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
                    mesh%map(i, j) = tool_index(set, map(((i + floor(length / scale) - length / 2) - 1) * floor(length / scale) &
                                         + (j - length / 2 - mod(length, 2))))
                else if (i >= length / 2 + mod(length, 2) .and. j <= length / 2) then
                    mesh%map(i, j) = tool_index(set, map(((i - length / 2) - mod(length, 2)) * floor(length / scale) &
                                         + (j + floor(length / scale) - length / 2)))
                else if (i <= length / 2 .and. j <= length / 2) then
                    mesh%map(i, j) = tool_index(set, map(((i + floor(length / scale) - length / 2) - 1) * floor(length / scale) &
                                         + (j + floor(length / scale) - length / 2)))
                else if (i > length / 2 .and. j > length / 2) then
                    mesh%map(i, j) = tool_index(set, map(((i - length / 2) - 1) * floor(length / scale) &
                                         + (j - length / 2)))
                end if
                mesh%weight(mesh%map(i, j)) = mesh%weight(mesh%map(i, j)) + 1
            end do
        end do
        
        do i = 1, mesh%num_reduced_points
            mesh%reduced_points(1, i) = 2 * pi * dble(grid_point(1, set(i) + 1)) / floor(length / scale)
            mesh%reduced_points(2, i) = 2 * pi * dble(grid_point(2, set(i) + 1)) / floor(length / scale)
        end do
    
    end function type_k_mesh_crystal
    
    
    function type_serialise_crystal(crystal) result(string)
        type(type_crystal), intent(in) :: crystal
        character(:), allocatable      :: string
        string = string + tool_token_write("crystal_length",          tool_to_string(crystal%crystal_length)) + "\n"
        string = string + tool_token_write("num_bands",               tool_to_string(crystal%num_bands)) + "\n"
        string = string + tool_token_write("max_hopping",             tool_to_string(crystal%max_hopping)) + "\n"
        string = string + tool_token_write("conducting_band_minimum_kx", tool_to_string(crystal%conducting_band_minimum_kx)) + "\n"
        string = string + tool_token_write("conducting_band_minimum_ky", tool_to_string(crystal%conducting_band_minimum_ky)) + "\n"
        
        string = string + tool_token_write("lattice_vectors",         tool_to_string(pack(crystal%lattice_vectors, .true.))) + "\n"
        string = string + tool_token_write("atom_types",              tool_to_string(     crystal%atom_types)) + "\n"
        string = string + tool_token_write("atom_positions",          tool_to_string(pack(crystal%atom_positions, .true.))) + "\n"
        string = string + tool_token_write("tunneling_weight",        tool_to_string(pack(crystal%tunneling_weight, .true.))) + "\n"
        string = string + tool_token_write("tunneling",               tool_to_string(pack(crystal%tunneling, .true.)))
        
    end function type_serialise_crystal
    
    function type_deserialise_crystal(string) result(crystal)
        character(*), intent(in) :: string
        integer, allocatable     :: temp_integer(:)
        type(type_crystal)       :: crystal
        crystal%crystal_length             = minval(tool_to_real(   tool_token_read(string, "crystal_length")), dim = 1)
        crystal%num_bands                  = minval(tool_to_integer(tool_token_read(string, "num_bands")), dim = 1)
        crystal%max_hopping                = minval(tool_to_integer(tool_token_read(string, "max_hopping")), dim = 1)
        crystal%conducting_band_minimum_kx = minval(tool_to_real(   tool_token_read(string, "conducting_band_minimum_kx")), dim = 1)
        crystal%conducting_band_minimum_ky = minval(tool_to_real(   tool_token_read(string, "conducting_band_minimum_ky")), dim = 1)
        
        crystal%lattice_vectors            = reshape(tool_to_real(   tool_token_read(string, "lattice_vectors")), (/ 3, 3 /))
        temp_integer                       =         tool_to_integer(tool_token_read(string, "atom_types"))
        allocate(crystal%atom_types(      size(temp_integer)))
        crystal%atom_types                 = temp_integer 
        allocate(crystal%atom_positions(  3, size(crystal%atom_types)))            
        crystal%atom_positions             = reshape(tool_to_real(   tool_token_read(string, "atom_positions")), &
                                                     (/ 3, size(crystal%atom_types) /))
        allocate(crystal%tunneling_weight(2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1))
        crystal%tunneling_weight           = reshape(tool_to_integer(tool_token_read(string, "tunneling_weight")), &
                                                     (/ 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                                        2 * crystal%max_hopping + 1 /))
        allocate(crystal%tunneling(       2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                          crystal%num_bands, crystal%num_bands))                          
        crystal%tunneling                  = reshape(tool_to_complex(tool_token_read(string, "tunneling")), &
                                                     (/ 2 * crystal%max_hopping + 1, 2 * crystal%max_hopping + 1, &
                                                        2 * crystal%max_hopping + 1, crystal%num_bands, crystal%num_bands /))
                                                        
    end function type_deserialise_crystal
    
    function type_serialise_heterostructure(structure) result(string)
        type(type_heterostructure), intent(in) :: structure
        character(:), allocatable              :: string
        string = string + tool_token_write("num_layers", tool_to_string(structure%num_layers)) + "\n"
        string = string + tool_token_write("min_layer",  tool_to_string(structure%min_layer)) + "\n"
        string = string + tool_token_write("max_layer",  tool_to_string(structure%max_layer)) + "\n"
        string = string + tool_token_write("min_state",  tool_to_string(structure%min_state)) + "\n"
        string = string + tool_token_write("max_state",  tool_to_string(structure%max_state)) + "\n"
        
        string = string + tool_token_write("potential",  tool_to_string(structure%potential)) + "\n"
        string = string + tool_token_write("material",   "{" + type_serialise(structure%material) + "}")
    
    end function type_serialise_heterostructure
    
    function type_deserialise_heterostructure(string) result(structure)
        character(*), intent(in)   :: string
        type(type_heterostructure) :: structure
        structure%num_layers = minval(tool_to_real(tool_token_read(string, "num_layers")), dim = 1)
        structure%min_layer = minval(tool_to_real(tool_token_read(string, "min_layer")), dim = 1)
        structure%max_layer = minval(tool_to_real(tool_token_read(string, "max_layer")), dim = 1)
        structure%min_state = minval(tool_to_real(tool_token_read(string, "min_state")), dim = 1)
        structure%max_state = minval(tool_to_real(tool_token_read(string, "max_state")), dim = 1)
        
        allocate(structure%potential(structure%num_layers))
        structure%potential = tool_to_real(tool_token_read(string, "potential"))
        structure%material  = type_deserialise_crystal(tool_token_read(string, "material"))
                                                        
    end function type_deserialise_heterostructure
    
    
    function type_serialise_k_mesh(mesh) result(string)
        type(type_k_mesh), intent(in) :: mesh
        character(:), allocatable     :: string
        string = string + tool_token_write("length",             tool_to_string(mesh%length)) + "\n"
        string = string + tool_token_write("scale",              tool_to_string(mesh%scale)) + "\n"
        string = string + tool_token_write("num_reduced_points", tool_to_string(mesh%num_reduced_points)) + "\n"
        
        string = string + tool_token_write("reduced_points",     tool_to_string(pack(mesh%reduced_points, .true.))) + "\n"
        string = string + tool_token_write("weight",             tool_to_string(mesh%weight)) + "\n"
        string = string + tool_token_write("map",                tool_to_string(pack(mesh%map, .true.)))
    
    end function type_serialise_k_mesh
    
    function type_deserialise_k_mesh(string) result(mesh)
        character(*), intent(in)   :: string
        type(type_k_mesh)          :: mesh
        mesh%length             = minval(tool_to_integer(tool_token_read(string, "length")), dim = 1)
        mesh%scale              = minval(tool_to_real(   tool_token_read(string, "scale")), dim = 1)
        mesh%num_reduced_points = minval(tool_to_integer(tool_token_read(string, "num_reduced_points")), dim = 1)
        
        allocate(mesh%reduced_points(2, mesh%num_reduced_points))
        allocate(mesh%weight(mesh%num_reduced_points))
        allocate(mesh%map(mesh%length, mesh%length))
        
        mesh%reduced_points = reshape(tool_to_real(tool_token_read(   string, "reduced_points")), shape(mesh%reduced_points))
        mesh%weight         = reshape(tool_to_integer(tool_token_read(string, "weight")), shape(mesh%weight))
        mesh%map            = reshape(tool_to_integer(tool_token_read(string, "map")), shape(mesh%map))
                                                        
    end function type_deserialise_k_mesh
    

end module type
