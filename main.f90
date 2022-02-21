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
    
    integer,        parameter   :: num_variation      = 20 !25
    integer,        parameter   :: num_layers         = 50 !150
    
    integer                     :: well_1_start       = 1 
    integer                     :: well_1_stop        = 1 
    real*8                      :: well_1_start_depth = 0d0
    real*8                      :: well_1_stop_depth  = 0d0
    
    integer                     :: well_2_start       = 1!6
    integer                     :: well_2_stop        = 1!45!20
    real*8                      :: well_2_start_depth = -0.5d0 !-0.5d0
    real*8                      :: well_2_stop_depth  = 0d0

    real*8                      :: c_parameter = 0.09 !0.01 +0.5
    
    ! for varying depth
    real*8                      :: b_parameter = 0.645!0.6
    
    ! for varying width 
    real*8                      :: a_parameter = 2.5
        
    
    integer,        parameter   :: num_k_length       = 100
    integer,        parameter   :: num_energy         = 100
    real*8,         parameter   :: length_scale       = 1d0
    real*8,         parameter   :: broadening         = 0.0025d0    
    
    integer                     :: min_layer          = -1
    integer                     :: max_layer          = -1
    integer                     :: min_band           = -1
    integer                     :: max_band           = -1
    integer                     :: min_state          = -1
    integer                     :: max_state          = 10
    real*8                      :: energy_min         = -1d0
    real*8                      :: energy_max         = 0d0
    
    
    integer                     :: x_offset           = 0
    character(*),   parameter   :: x_label            = "Surface Layer Doping [-0.1eV]"
    
    character(*),   parameter   :: input_file         = "SrTiO3_hr.dat"
    character(*),   parameter   :: out_dir            = "./out/"
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
    
    integer                     :: max_hopping, num_states, max_num_states(num_variation), num_bands
    
    real*8,         allocatable :: kx(:), ky(:)
    integer,        allocatable :: kw(:)
    integer,        allocatable :: kp(:), kl(:)
    integer                     :: k_start, k_stop
    
    real*8                      :: pot(num_variation, num_layers)
    
    real*8,         allocatable :: energy(:)
    complex*16,     allocatable :: weight(:, :)
    integer                     :: num_found
    real*8                      :: cbm
    
    real*8                      :: energy_range(num_variation, num_energy)
    
    real*8,         allocatable :: heatmap_cpu(:, :, :, :), heatmap(:, :, :, :)
    real*8,         allocatable :: den_cpu(:, :, :), den(:, :, :)
    real*8,         allocatable :: contribution_cpu(:, :, :) , contribution(:, :, :)                             
    real*8,         allocatable :: energymap_cpu(:, :), energymap(:, :)
    
    character(256)              :: path, variation_dir, band_dir, energy_state_dir, energymap_dir
    
    integer                     :: i, j, k, l, m, e, n, start_time
    
    integer                     :: gamma_k_pos, cpu_from, stat
    
    real*8,         allocatable :: gamma_energy(:), gamma_energy_list(:,:)
    
    logical                     :: exists
    
    CHARACTER(len=255) :: cwd
    
    real*8,         allocatable :: variation_parameter_list(:)
    
    real*8                      :: width, depth
    
    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Initilise..."
        
    end if

    allocate(variation_parameter_list(num_variation))
    variation_parameter_list = 0d0
    
    ! Set Up potentials
    pot = 0d0
    do i = 1, num_variation
        ! Manipulate Wells
        ! -------------------------------------------------------------------------------------------------------------------------
        !well_1_start_depth = (i ) * -0.1d0
        
        
        !!! variation of potential (depth)
        !well_2_start_depth = i * -0.025d0
        !variation_parameter_list(i) = well_2_start_depth  
        
        !!! variation of width 
        well_2_stop = i
        variation_parameter_list(i) = well_2_stop
        
        
        ! for curve potential (vary width)
        ! a_parameter = 2.0 +  i * 0.028  ! 1.9 - 2.6  (12 - 42 layers ) for depth = -0.3
        !a_parameter = 1.5 +  i * 0.028  ! 1.9 - 2.6  (12 - 42 layers ) for depth = -0.4
        !a_parameter = 1.3 +  i * 0.028  ! 1.9 - 2.6  (12 - 42 layers ) for depth = -0.5
        
        !a_parameter = 1.5 +  i * 0.015  ! 1.9 - 2.6  (12 - 42 layers ) for depth = -0.5
        
        ! for curve potential (vary depth)
        
        !b_parameter = 0.6 + i * 0.00225
        
        
        
        ! -------------------------------------------------------------------------------------------------------------------------
        
        ! Set Up Potential
        ! for normal potential
        !call potential_add_well(pot(i, :), well_1_start, well_1_stop, well_1_start_depth, well_1_stop_depth)
        call potential_add_well(pot(i, :), well_2_start, well_2_stop, well_2_start_depth, well_2_stop_depth)
        
        ! for curve potential (vary width)
        !call potential_add_well_curve_width(pot(i, :), well_1_start, well_1_start_depth, a_parameter,c_parameter, width)
        !call potential_add_well_curve_width (pot(i, :), well_2_start, well_2_start_depth, a_parameter, c_parameter, width)  
        
        ! for curve potential (vary depth)
        !call potential_add_well_curve_depth(pot(i, :), well_1_start, well_1_stop, b_parameter,c_parameter, depth)
        !call potential_add_well_curve_depth (pot(i, :), well_2_start, well_2_stop, b_parameter, c_parameter, depth)       

        ! for varying width --> curve potential 
        !if (cpu_is_master()) then
        !    variation_parameter_list(i) = width 
        !end if
        
        ! for varying depth --> curve potential 
        !if (cpu_is_master()) then
        !    variation_parameter_list(i) = depth
        !end if
             
    end do
    
    ! Extract Data
    if (cpu_is_master()) then
        call import_data(input_file, hr, hrw, max_hopping, num_bands)
    end if
    call cpu_broadcast(max_hopping, 1)
    call cpu_broadcast(num_bands, 1)
    if (.not. cpu_is_master()) then
        allocate(hr(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1, num_bands, num_bands))
        allocate(hrw(2 * max_hopping + 1, 2 * max_hopping + 1, 2 * max_hopping + 1))
    end if
    call cpu_broadcast(hr, size(hr))
    call cpu_broadcast(hrw, size(hrw))
    
    ! Build Route
    if (cpu_is_master()) then
        call route_build(kx, ky, kw, kp, num_k_length, length_scale, &
            lattice, positions, atom_types)
        i = size(kw)
        j = size(kp)
    end if

    call cpu_broadcast(i, 1)
    call cpu_broadcast(j, 1)
    if (.not. cpu_is_master()) then
        allocate(kx(i))
        allocate(ky(i))
        allocate(kw(i))
        allocate(kp(j))
    end if
    call cpu_broadcast(kx, size(kx))
    call cpu_broadcast(ky, size(ky))
    call cpu_broadcast(kw, size(kw))
    call cpu_broadcast(kp, size(kp))
    
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
    do i = 1, num_variation
        if (energy_min == -1d0 + cbm) energy_min = hqtlib_find_energy_min(hr, hrw, pot(i, :), broadening)
        if (energy_max == -1d0 + cbm) energy_max = hqtlib_find_energy_max(hr, hrw, pot(i, :), length_scale, broadening)
        max_num_states(i) = hqtlib_find_max_num_states(hr, hrw, pot(i, :), energy_min, energy_max)
        if (max_num_states(i) == 0) then
            max_num_states(i) = 1
        end if
        if (max_state .ne. -1 .and. max_num_states(i) > max_state) then
            max_num_states(i) = max_state
        end if
        energy_range(i, :) = route_range_double(energy_min, energy_max, num_energy) - cbm
    end do
    
    ! Determine bounds
    if (max_state > maxval(max_num_states)) max_state = maxval(max_num_states)
    if (max_layer == -1) max_layer = num_layers
    if (min_layer == -1) min_layer = 1
    if (max_band == -1) max_band = num_bands
    if (min_band == -1) min_band = 1
    if (max_state == -1) max_state = maxval(max_num_states)
    if (min_state == -1) min_state = 1
    
    ! Allocate Dynamic Arrays
    allocate(energy(             num_bands * num_layers))
    allocate(weight(             num_bands * num_layers, num_bands * num_layers))
    
    allocate(heatmap_cpu(        size(kp),   max_band - min_band + 1, max_state - min_state + 1, num_energy))
    allocate(den_cpu(            num_layers, max_band - min_band + 1, max_state - min_state + 1))
    allocate(energymap_cpu(      size(kp),   max_state - min_state + 1))
    allocate(contribution_cpu(   size(kp),   max_band - min_band + 1, max_state - min_state + 1))
    if (cpu_is_master()) then
        allocate(heatmap(        size(kp),   max_band - min_band + 1, max_state - min_state + 1, num_energy))
        allocate(den(            num_layers, max_band - min_band + 1, max_state - min_state + 1))
        allocate(energymap(      size(kp),   max_state - min_state + 1))
        allocate(contribution(   size(kp),   max_band - min_band + 1, max_state - min_state + 1))
    end if
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir, "Run ")
        call export_vstack(trim(path)//"density.dat", &
            "#          Total  Quantum Well 1  Quantum Well 2             Sum"// &
            "       Avg Total      Avg Well 1      Avg Well 2    Avg Well Sum")
        call export_vstack(trim(path)//"variation_parameter_list.dat", variation_parameter_list)
    end if
    
    ! Determine CPU Work
    call cpu_split_work(k_start, k_stop, size(kw))
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (cpu_is_master()) then
        allocate(gamma_energy_list(num_variation, num_layers * num_bands))
        gamma_energy_list = 0d0
    end if 
    
    allocate(gamma_energy(num_layers* num_bands))
    
    do i = 1, num_variation
        den_cpu          = 0d0
        heatmap_cpu      = 0d0
        energymap_cpu    = 0d0
        contribution_cpu = 0d0
        if (cpu_is_master()) then
            den          = 0d0
            heatmap      = 0d0
            energymap    = 0d0
            contribution = 0d0
        end if
        
        ! determine gamma point 
        gamma_k_pos  = size(kw)
        gamma_energy = 0d0
        ! Spectral into Density and Heatmap
        if (cpu_is_master()) then
            write(*, fmt = "(A)", advance = "no") "  Generate Data: "
            call system_clock(start_time)
        end if
        do k = k_start, k_stop
            call hqtlib_find_energy_and_weight(energy, weight, num_found, kx(k), ky(k), hr, hrw, pot(i, :), 1, max_num_states(i))
            energy(:num_found) = energy(:num_found) - cbm
            
            if (k == size(kw)) then     
                do j = 1, num_found
                    gamma_energy(j) = energy(j)    
                end do
                call cpu_send_double(gamma_energy, size(gamma_energy), 0, i)
                cpu_from = cpu_get_id()
                
            end if 
            
            ! Density
            den_cpu = den_cpu + kw(k) * hqtlib_find_density(energy(:num_found), weight(:, :num_found), energy_range(i, :), &
                                                            broadening, num_bands, crystal_length, sum(kw), &
                                                            min_layer, max_layer, min_band, max_band, min_state, max_state)
            
            if (kl(k) .ne. -1) then ! True if on Path
                ! Heatmap
                heatmap_cpu(kl(k), :, :, :) = sum(hqtlib_find_spectral_distribution(energy(:num_found), weight(:, :num_found), &
                                                                                    energy_range(i, :), broadening, num_bands, &
                                                              min_layer, max_layer, min_band, max_band, min_state, max_state), 1)
                
                ! Energy Structure
                energymap_cpu(kl(k), :num_found)  = energy(:num_found)
                ! Band Contribution
                contribution_cpu(kl(k), :, :) = sum(hqtlib_find_band_contribution(energy(:num_found), weight(:, :num_found), &
                                                                                  num_bands, &
                                                            min_layer, max_layer, min_band, max_band, min_state, max_state), 1)
            end if
            
            if (cpu_is_master()) call main_status(k - k_start + 1, k_stop - k_start + 1, start_time)
        end do
        if (cpu_is_master()) then
            call main_clear_status()
            write(*, fmt = "(A)") "100%"
        end if
        
        if (cpu_is_master()) then
            call cpu_recv_double(gamma_energy, size(gamma_energy),  MPI_ANY_SOURCE , i)
            gamma_energy_list(i,:) = gamma_energy
        end if 
        
        den_cpu        = den_cpu * 1d-18
        call cpu_sum(den_cpu, den)
        call cpu_sum(heatmap_cpu, heatmap)
        call cpu_sum(energymap_cpu, energymap)
        call cpu_sum(contribution_cpu, contribution)
        
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            if (cpu_is_master()) then
                write(*, fmt = "(A)", advance = "no") "  Plot Data:     "
                call system_clock(start_time)
            end if
            ! Export Meta Data
            variation_dir = export_create_dir(trim(path), "Variation "//export_to_string(i))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                "#      Num K Len      Num Layers      Num Path K Num Path Energy")
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ num_k_length, num_layers,     size(kp),       num_energy /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                "#     Broadening Crystal Len (a)   K Len [2pi*a]")
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ broadening, crystal_length * 1d10, length_scale /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                "#            Min             Max")
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ minval(pot(i, :)),                    maxval(pot(i, :)) /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ minval(den),                          maxval(den) /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ minval(energy_range(i, :)),           maxval(energy_range(i, :)) /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/ minval(log(sum(sum(heatmap, 2), 2))), maxval(log(sum(sum(heatmap, 2), 2))) /))
            call export_hstack(trim(variation_dir)//"meta.dat", &
                (/ "            ", "            ", "            ", "            ", "            ", &
                " # Potential", " # Density  ", " # Heatmap y", " # Heatmap c" /))
              
            ! Export and Plot Variation Data
            call export_hstack(trim(variation_dir)//"gamma_energy.dat",    gamma_energy)
            call export_hstack(trim(variation_dir)//"potential.dat",       pot(i, :))
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
                              column_label_1 = "1", &
                              column_label_2 = "2", &
                              column_label_3 = "3")
            call export_hstack(trim(variation_dir)//"band structure.dat",  log(sum(sum(heatmap, 2), 2)))                  
            call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(variation_dir))
            energymap_dir = export_create_dir(trim(variation_dir), "Energy_map")
            do j = 1, maxval(max_num_states)
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
                              column_label_1 = "1", &
                              column_label_2 = "2", &
                              column_label_3 = "3", &
                              k_scale = length_scale, &
                              min_y = minval(energy_range(i, :)), &
                              max_y = maxval(energy_range(i, :)))
            
            ! Export and Plot Variation Band Data
            do j = 1, max_band - min_band + 1
                band_dir = export_create_dir(trim(variation_dir)//"Bands/", "Band "//export_to_string(j))
                call export_hstack(trim(band_dir)//"density"//".dat", sum(den(:, j, :), 2))
                call export_hstack(trim(band_dir)//"band structure"//".dat", log(sum(heatmap(:, j, :, :), 2)))
                call graph_basic_plot("density", "density", 1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(band_dir))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(band_dir))
            end do
            
            ! Export and Plot Variation Energy State Data
            do j = 1, max_num_states(i) - min_state + 1
                energy_state_dir = export_create_dir(trim(variation_dir)//"Energy States/", "Energy State "//export_to_string(j))
                call export_hstack(trim(energy_state_dir)//"density.dat", sum(den(:, :, j), 2))
                call export_hstack(trim(energy_state_dir)//"density.dat", transpose(sort_normalise(den(:, :, j))))
                call graph_colour(data_folder = trim(energy_state_dir)//"density.dat", &
                                  output_file = trim(energy_state_dir)//"density", &
                                  x_label = "Layer", &
                                  y_label = "Carrier Density [nm^{-2}]", &
                                  x_triangle = 0.75d0, &
                                  y_triangle = 0.75d0, &
                                  triangle_size = 0.15d0, &
                                  column_label_1 = "1", &
                                  column_label_2 = "2", &
                                  column_label_3 = "3")
                call export_hstack(trim(energy_state_dir)//"band structure"//".dat", log(sum(heatmap(:, :, j, :), 2)))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(energy_state_dir))
                ! Bands
                do n = 1, max_band - min_band + 1
                    band_dir = export_create_dir(trim(energy_state_dir)//"Bands/", "Band " &
                        //export_to_string(n))
                    call export_hstack(trim(band_dir)//"density"//".dat", &
                        den(:, n, j))
                    call export_hstack(trim(band_dir)//"band structure"//".dat", &
                        log(heatmap(:, n, j, :)))
                    call graph_basic_plot("density", "density", &
                        1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(band_dir))
                    call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(band_dir))
                end do
                ! Energy Map
                energymap_dir = export_create_dir(trim(energy_state_dir), "Energy_map")
                call export_hstack(trim(energymap_dir)//"Energy_level.dat", energymap(:, j))
                call export_hstack(trim(energymap_dir)//"Energy_level.dat", transpose(contribution(:, :, j)))
                call graph_colour_MGX(data_folder = trim(energymap_dir), &
                              output_file = trim(energymap_dir), &
                              x_label = "K_{M Γ X} [π/a]", &
                              y_label = "E - E_{cbm} [ev]", &
                              x_triangle = 0.2d0, &
                              y_triangle = 0.2d0, &
                              triangle_size = 0.15d0, &
                              column_label_1 = "1", &
                              column_label_2 = "2", &
                              column_label_3 = "3", &
                              k_scale = length_scale, &
                              min_y = minval(energy_range(i, :)), &
                              max_y = maxval(energy_range(i, :)))
                if (cpu_is_master()) call main_status(j, max_state - min_state + 1, start_time)
            end do
            if (cpu_is_master()) then
                call main_clear_status()
                write(*, fmt = "(A)") "100%"
            end if
            
            ! Export and Plot Data
            call export_vstack(trim(path)//"density.dat", reshape((/ &
                sum(den), &
                sum(den(well_1_start:well_1_stop, :, :)), &
                sum(den(well_2_start:well_2_stop, :, :)), &
                sum(den(well_1_start:well_1_stop, :, :)) + &
                    sum(den(well_2_start:well_2_stop, :, :)), &
                sum(den) / size(den), &
                sum(den(well_1_start:well_1_stop, :, :)) / (abs(well_1_stop - well_1_start) + 1), &
                sum(den(well_2_start:well_2_stop, :, :)) / (abs(well_2_stop - well_2_start) + 1), &
                sum(den(well_1_start:well_1_stop, :, :)) + &
                    sum(den(well_2_start:well_2_stop, :, :)) / (abs(well_1_stop + well_2_stop - well_1_start - well_2_start) + 2) &
                /), (/ 1, 8 /)))
            if (i > 2) then
                call graph_basic_plot("density", "Total_Density", &
                    1, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "reservoir_density", &
                    2, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "transport_density", &
                    3, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "reservoir_plus_transport_density", &
                    4, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_Total_Density", &
                    5, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_reservoir_density", &
                    6, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_transport_density", &
                    7, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                call graph_basic_plot("density", "Average_reservoir_plus_transport_density", &
                    8, x_label, "Carrier Density [nm^{-2}]", x_offset, trim(path))
                


                do j = 1, i 
                    gamma_energy_list(j, 11) = variation_parameter_list(j)
                end do 
                call export_vstack(trim(path)//"gamma_energy_list.dat", gamma_energy_list(1:i,1:11))

                !!! for variation of potential (depth)
                !call graph_multiple_plot("gamma_energy_list", "gamma_energy_plot", "Qunatum well potential [eV]", &
                !"energy states [eV]", 1, trim(path))
                
                call graph_multiple_plot("gamma_energy_list", "gamma_energy_plot", "Qunatum well potential [eV]", &
                "energy states [eV]", 0, trim(path))
                
                
                
                
                inquire(file = trim(path)//"gamma_energy_list.dat", exist = exists)
                if(exists) then 

                    open(unit=1234, file=trim(path)//"gamma_energy_list.dat", status='old')
                    close(1234, status='delete')
                end if 
                
                if (i== num_variation) then 
                
                    do j = 1, num_variation
                        gamma_energy_list(j, 11) = variation_parameter_list(j)
                    end do 
                    
                    call export_vstack(trim(path)//"gamma_energy_list.dat", gamma_energy_list(1:i,1:11))              
                
                end if 
                
                
                
            end if
            
            write(*, fmt = "(A10, I4.1, A9)") "Variation ", i, " complete"
            call timer_split()
        end if
    end do
    
    !if (cpu_is_master()) then

        !call export_data(trim(path)//"gamma_energy_list.dat", gamma_energy_list(:,1:5))

        !call graph_multiple_plot("gamma_energy_list", "gamma_energy_plot", "potential [eV]", "energy states [eV]", 1, trim(path))

    !end if 
    
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
    end subroutine main_clear_status

end program main
