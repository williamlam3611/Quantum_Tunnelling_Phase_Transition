program main
    use timer
    use cpu
    use hqtlib
    use import
    use route
    use potential
    use transform
    use bulk
    use matrix
    use sort
    use spectral
    use density
    use export
    use graph
    use mpi
    implicit none
    
    integer,        parameter   :: num_layers         = 50
    
    integer                     :: well_start         = 1
    integer                     :: well_stop          = 10
    real*8                      :: depth              = -0.1d0
        
    integer,        parameter   :: num_k_length       = 50
    integer,        parameter   :: num_energy         = 256
    real*8,         parameter   :: length_scale       = 1d0
    real*8,         parameter   :: broadening         = 0.0025d0    
    
    integer                     :: min_layer          = -1
    integer                     :: max_layer          = -1
    integer                     :: min_band           = -1
    integer                     :: max_band           = -1
    integer                     :: min_state          = -1
    integer                     :: max_state          = -1
    integer                     :: temp_max_states    = -1
    real*8                      :: energy_min         = -1d0
    real*8                      :: energy_max         = 0d0
    
    character(*),   parameter   :: input_file         = "SrTiO3_hr.dat"
    character(*),   parameter   :: out_dir            = "./out/poisson/"
    real*8,         parameter   :: crystal_length     = 3.905d-10
    
    real*8,         parameter   :: lattice(3, 3)      = reshape((/ 1d0, 0d0, 0d0, &
                                                                0d0, 1d0, 0d0, &
                                                                0d0, 0d0, 1d0 /), &
                                                                shape(lattice))
    real*8,         parameter   :: positions(3, 5)    = reshape((/ 0d0  , 0d0  , 0d0  , &
                                                                0.5d0, 0.5d0, 0.5d0, &
                                                                0.5d0, 0.5d0, 0d0  , &
                                                                0.5d0, 0d0  , 0.5d0, &
                                                                0d0  , 0.5d0, 0.5d0 /), &
                                                               shape(positions))
    integer,        parameter   :: atom_types(*)      = (/ 1, 2, 3 ,3, 3/)
    
    complex*16,     allocatable :: hr(:, :, :, :, :)
    integer,        allocatable :: hrw(:, :, :)
    
    integer                     :: max_hopping, num_states, num_bands
    
    real*8,         allocatable :: kx(:), ky(:)
    integer,        allocatable :: kw(:)
    integer,        allocatable :: kp(:), kl(:)
    integer                     :: k_start, k_stop
    
    real*8                      :: pot(num_layers)
    
    integer                     :: num_found
    real*8                      :: cbm
    
    real*8,         allocatable :: energy(:)
    complex*16,     allocatable :: weight(:, :)
    real*8                      :: energy_range(num_energy)
    
    real*8,         allocatable :: den_cpu(:, :, :), den(:, :, :)
    real*8,         allocatable :: contribution_cpu(:, :, :) , contribution(:, :, :)                             
    real*8,         allocatable :: energymap_cpu(:, :), energymap(:, :)
    
    character(256)              :: path, variation_dir, energymap_dir
    
    integer                     :: i, j, k, l, m, e, n, start_time
    
    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Initilise..."
        
    end if
    
    ! Set Up potential
    pot = 0d0
    call potential_add_well(pot, well_start, well_stop, depth, 0d0)
    
    ! Extract Data and build route
    call import_data(input_file, hr, hrw, max_hopping, num_bands)
    call route_build(kx, ky, kw, kp, num_k_length, length_scale, &
                     lattice, positions, atom_types)
    
    ! Produce inverse path mapping kw -> kp
    allocate(kl(size(kw)))
    kl = -1
    do i = 1, size(kp)
        kl(kp(i)) = i
    end do
    
    ! Determine Conductiong Band Minimum
    cbm = hqtlib_find_energy_min(hr, hrw, num_layers)
    energy_min = energy_min + cbm
    energy_max = energy_max + cbm
    
    ! Determine Maximum number of bands
    if (energy_min == -1d0 + cbm) energy_min = hqtlib_find_energy_min(hr, hrw, pot, broadening)
    if (energy_max == -1d0 + cbm) energy_max = hqtlib_find_energy_max(hr, hrw, pot, length_scale, broadening)
    energy_range = route_range_double(energy_min, energy_max, num_energy) - cbm
    
    ! Determine bounds
    if (max_layer == -1) max_layer = num_layers
    if (min_layer == -1) min_layer = 1
    if (max_band == -1) max_band = num_bands
    if (min_band == -1) min_band = 1
    if (max_state == -1) max_state = -1
    if (min_state == -1) min_state = 1
    
    ! Allocate Dynamic Arrays
    allocate(energy(             num_bands * num_layers))
    allocate(weight(             num_bands * num_layers, num_bands * num_layers))
    
    allocate(den_cpu(            num_layers, max_band - min_band + 1, max_state - min_state + 1))
    allocate(energymap_cpu(      size(kp),   max_state - min_state + 1))
    allocate(contribution_cpu(   size(kp),   max_band - min_band + 1, max_state - min_state + 1))
    if (cpu_is_master()) then
        allocate(den(            num_layers, max_band - min_band + 1, max_state - min_state + 1))
        allocate(energymap(      size(kp),   max_state - min_state + 1))
        allocate(contribution(   size(kp),   max_band - min_band + 1, max_state - min_state + 1))
    end if
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir, "Run ")
        call export_vstack(trim(path)//"density.dat", &
            "#          Total    Quantum Well")
    end if
    
    ! Determine CPU Work
    call cpu_split_work(k_start, k_stop, size(kw))
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    
    do i = 1, 10
        den_cpu          = 0d0
        energymap_cpu    = 0d0
        contribution_cpu = 0d0
        if (cpu_is_master()) then
            den          = 0d0
            energymap    = 0d0
            contribution = 0d0
        end if
        
        ! Spectral into Density and Heatmap
        if (cpu_is_master()) then
            write(*, fmt = "(A)", advance = "no") "  Generate Data: "
            call system_clock(start_time)
        end if
        
        temp_max_states = max_state
        if (temp_max_states == -1) temp_max_states = hqtlib_find_max_num_states(hr, hrw, pot, energy_min, energy_max)
        do k = k_start, k_stop
            call hqtlib_find_energy_and_weight(energy, weight, num_found, kx(k), ky(k), hr, hrw, pot, min_state, temp_max_states)
            energy(:num_found) = energy(:num_found) - cbm
            
            ! Density
            den_cpu = den_cpu + kw(k) * hqtlib_find_density(energy(:num_found), weight(:, :num_found), energy_range, &
                                                            broadening, num_bands, crystal_length, sum(kw), &
                                                            min_layer, max_layer, min_band, max_band, min_state, temp_max_states)
            
            
            if (kl(k) .ne. -1) then ! True if on Path
                
                ! Energy Structure
                energymap_cpu(kl(k) , :num_found)  = energy(:num_found)
                
                ! Band Contribution
                contribution_cpu(kl(k), :, :) = sum(hqtlib_find_band_contribution(energy(:num_found), weight(:, :num_found), &
                                                                                  num_bands, &
                                                            min_layer, max_layer, min_band, max_band, min_state, temp_max_states),1)
            end if
            
            if (cpu_is_master()) call main_status(k - k_start + 1, k_stop - k_start + 1, start_time)
        end do
        if (cpu_is_master()) then
            call main_clear_status()
            write(*, fmt = "(A)") "100%"
        end if
        
        den_cpu        = den_cpu * 1d-18
        call cpu_sum(den_cpu, den)
        call cpu_sum(energymap_cpu, energymap)
        call cpu_sum(contribution_cpu, contribution)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            variation_dir = export_create_dir(trim(path), "Variation "//export_to_string(i))
            ! Export and Plot Variation Data
            call export_hstack(trim(variation_dir)//"potential.dat", pot)
            call graph_basic_plot("potential", "potential", 1, "Layer", "Potential [eV]", 1, trim(variation_dir))
            call export_hstack(trim(variation_dir)//"density.dat",         sum(sum(den, 2), 2))
            call export_hstack(trim(variation_dir)//"density.dat",         transpose(sort_normalise(sum(den, 3))))
            call graph_colour(data_folder = trim(variation_dir)//"density.dat", &
                              output_file = trim(variation_dir)//"density", &
                              x_label = "Layer", &
                              y_label = "Carrier Density [nm^{-2}]", &
                              x_triangle = 0.75d0, &
                              y_triangle = 0.75d0, &
                              triangle_size = 0.15d0, &
                              column_label_1 = "d_{xy}", &
                              column_label_2 = "d_{xz}", &
                              column_label_3 = "d_{yz}")
            energymap_dir = export_create_dir(trim(variation_dir), "Energy_map")
            do j = 1, temp_max_states
                call export_hstack(trim(energymap_dir)//"Energy_level"//export_to_string(j)//".dat", energymap(:, j))
                call export_hstack(trim(energymap_dir)//"Energy_level"//export_to_string(j)//".dat", &
                                   transpose(contribution(:, :, j)))
            end do
            call graph_colour_MGX(data_folder = trim(energymap_dir), &
                              output_file = trim(energymap_dir), &
                              x_label = "K_{M Γ X} [π/a]", &
                              y_label = "E - E_{cbm} [ev]", &
                              x_triangle = 0.2d0, &
                              y_triangle = 0.2d0, &
                              triangle_size = 0.15d0, &
                              column_label_1 = "d_{xy}", &
                              column_label_2 = "d_{xz}", &
                              column_label_3 = "d_{yz}", &
                              k_scale = length_scale, &
                              min_y = minval(energy_range), &
                              max_y = maxval(energy_range))
            
            ! Export and Plot Data
            call export_vstack(trim(path)//"density.dat", reshape((/ &
                sum(den), sum(den(well_start:well_stop, :, :)) /), (/ 1, 2 /)))
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do
    
    call cpu_stop()
    
contains
    subroutine main_status(i, total, start_time)
        integer, intent(in) :: i
        integer, intent(in) :: total
        integer, intent(in) :: start_time
        integer             :: tic, tps
        call system_clock(tic, tps)
        tic = tic - start_time
        if (i > 1) write(*, fmt = "(A27)", advance = "no") "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
        write(*, fmt = "(I3, A7, I16, A1)", advance = "no") &
            floor(dble(i * 100) / total), "% ETR: ", floor(dble(tic * (total - i)) / (tps * i)), "s"
    
    end subroutine main_status
        
    subroutine main_clear_status()
        write(*, fmt = "(A27)", advance = "no") "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
        write(*, fmt = "(A27)", advance = "no") "                           "
        write(*, fmt = "(A27)", advance = "no") "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
    end subroutine main_clear_status

end program main
