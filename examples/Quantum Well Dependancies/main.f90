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
    
    integer,        parameter   :: num_variation      = 20 !50
    integer,        parameter   :: num_layers         = 60! 120 ! 150
    
    integer                     :: well_start         = 1!6
    integer                     :: well_stop          = 10!45!20
    real*8                      :: well_start_depth   = -0.2d0 !-0.5d0
    real*8                      :: well_stop_depth    = 0! -0.1d0! 0d0

    real*8                      :: c_parameter = 0.09 !0.01 +0.5
    
    ! for varying depth
    real*8                      :: b_parameter = 0.645!0.6
    
    ! for varying width 
    real*8                      :: a_parameter = 2.5
        
    
    integer,        parameter   :: num_k_length       = 128! 256
    integer,        parameter   :: num_energy         = 128! 256
    real*8,         parameter   :: length_scale       = 1d0
    real*8,         parameter   :: broadening         = 0.0025d0    
    
    integer                     :: min_layer          = -1
    integer                     :: max_layer          = -1
    integer                     :: min_band           = -1
    integer                     :: max_band           = -1
    integer                     :: min_state          = -1
    integer                     :: max_state          = -1
    real*8                      :: energy_min         = -1d0
    real*8                      :: energy_max         = 0d0
    
    
    integer                     :: x_offset           = 0
    character(*),   parameter   :: x_label            = "Surface Layer Doping [-0.1eV]"
    
    character(*),   parameter   :: input_file         = "SrTiO3_hr.dat"
    character(*),   parameter   :: out_dir            = "./out/3d/"
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
    
    integer                     :: max_hopping, num_states, max_num_states(num_variation* num_variation), num_bands
    
    real*8,         allocatable :: kx(:), ky(:)
    integer,        allocatable :: kw(:)
    integer,        allocatable :: kp(:), kl(:)
    integer                     :: k_start, k_stop
    
    real*8                      :: pot(num_variation * num_variation , num_layers)
    
    real*8,         allocatable :: energy(:)
    complex*16,     allocatable :: weight(:, :)
    integer                     :: num_found
    real*8                      :: cbm
    
    real*8                      :: energy_range(num_variation* num_variation, num_energy)
    
    real*8,         allocatable :: heatmap_cpu(:, :, :, :), heatmap(:, :, :, :)
    real*8,         allocatable :: den_cpu(:, :, :), den(:, :, :)
    real*8,         allocatable :: contribution_cpu(:, :, :) , contribution(:, :, :)                             
    real*8,         allocatable :: energymap_cpu(:, :), energymap(:, :)
    
    character(256)              :: path, variation_dir, band_dir, energy_state_dir, energymap_dir, gamma_dir, &
                                   m_dir, x_dir, gx_dir, gm_dir
    
    integer                     :: i, j, k, l, m, e, n, start_time
    
    !!!
    ! for 3D plotting 
    real*8,         allocatable :: variation_parameter_width(:)
    real*8,         allocatable :: variation_parameter_depth(:)
    
    integer                     :: loop_count_i, nlines, io, file_number
    
    character(:), allocatable :: dir, dir_sub
    
    logical  , allocatable :: plotted 
    
    real*8,         allocatable :: data_copied(:,:,:)
    
    call cpu_start()
    if (cpu_is_master()) then
        call timer_start()
        print *, "Initilise..."
        
    end if
    
    allocate(variation_parameter_width(num_variation * num_variation))
    allocate(variation_parameter_depth(num_variation * num_variation))
    
    variation_parameter_width = 0d0
    variation_parameter_depth = 0d0
    
    ! Set Up potentials
    pot = 0d0
    

    do loop_count_i = 1, num_variation
        !well_start_depth = (loop_count_i-1) * -0.05!-0.005d0
        well_start_depth = -0.0125d0 * loop_count_i !-0.05!-0.005d0
        

        do i = 1, num_variation
            ! Manipulate Wells
            ! -------------------------------------------------------------------------------------------------------------------------
            !well_1_start_depth = (i ) * -0.1d0
            
            
            !!! variation of potential (depth)
            !well_start_depth = (i-1) * -0.005d0
            !variation_parameter_list(i) = well_start_depth  
            
            !!! variation of width 
            !well_stop = i -1 
            !variation_parameter_list(i) = well_stop
                        
            !! for 3D plotting
            well_stop = i + 3
            
            
            variation_parameter_depth( (loop_count_i -1) * num_variation + i  ) = well_start_depth
            variation_parameter_width( (loop_count_i -1) * num_variation + i ) = well_stop

            ! -------------------------------------------------------------------------------------------------------------------------
            
            ! Set Up Potential
            !call potential_add_well(pot(i, :), well_1_start, well_1_stop, well_1_start_depth, well_1_stop_depth)
            call potential_add_well(pot( (loop_count_i-1)*num_variation + i , :), well_start, well_stop, &
            well_start_depth, well_stop_depth)
                 
        end do
    end do 
        
    ! Extract Data
    call import_data(input_file, hr, hrw, max_hopping, num_bands)
    
    ! Build Route
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
    j = 0
    k = 0
    if (energy_min == -1d0 + cbm) j = 1
    if (energy_max == -1d0 + cbm) k = 1
    
    ! Determine Maximum number of bands
    do i = 1, num_variation* num_variation
        if (j == 1) energy_min = hqtlib_find_energy_min(hr, hrw, pot(i, :), broadening)
        if (k == 1) energy_max = hqtlib_find_energy_max(hr, hrw, pot(i, :), length_scale, broadening)
        energy_range(i, :) = route_range_double(&
                                merge(hqtlib_find_energy_min(hr, hrw, pot(i, :), broadening), energy_min, j == 1), &
                                merge(hqtlib_find_energy_max(hr, hrw, pot(i, :), length_scale, broadening), energy_max, k == 1), &
                                num_energy) - cbm
        max_num_states(i) = hqtlib_find_max_num_states(hr, hrw, pot(i, :), energy_min, energy_max)
        if (max_state .ne. -1 .and. max_num_states(i) > max_state) then
            max_num_states(i) = max_state
        end if
    end do
    
    ! Determine bounds
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
            "#          Total    Quantum Well")
        
        gamma_dir = export_create_dir(trim(path), "Gamma")
        m_dir     = export_create_dir(trim(path), "M")
        x_dir     = export_create_dir(trim(path), "X")
        gx_dir     = export_create_dir(trim(path), "GX")
        gm_dir     = export_create_dir(trim(path), "GM")
    end if
    
    ! Determine CPU Work
    call cpu_split_work(k_start, k_stop, size(kw))
    
    ! Analyse
    if (cpu_is_master()) then
        call timer_split()
        print *, "Begin Analysis..."
    end if
    
    do i = 1, num_variation * num_variation
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
        
        ! Spectral into Density and Heatmap
        if (cpu_is_master()) then
            write(*, fmt = "(A)", advance = "no") "  Generate Data: "
            call system_clock(start_time)
                        
        end if
        
        do k = k_start, k_stop
            if (max_num_states(i) .ne. 0)  then
                call hqtlib_find_energy_and_weight(energy, weight, num_found, kx(k), ky(k), &
                    hr, hrw, pot(i, :), min_state, max_num_states(i))
                energy(:num_found) = energy(:num_found) - cbm
                
                ! Density
                den_cpu = den_cpu + kw(k) * hqtlib_find_density(energy(:num_found), weight(:, :num_found), energy_range(i, :), &
                                                                broadening, num_bands, crystal_length, sum(kw), &
                                                                min_layer, max_layer, min_band, max_band, min_state, max_state)
                
                
                if (kl(k) .ne. -1) then ! True if on Path
                    ! Heatmap
                    heatmap_cpu(kl(k), :, :, :) = sum(hqtlib_find_spectral_distribution(energy(:num_found), &
                        weight(:, :num_found), energy_range(i, :), broadening, num_bands, &
                        min_layer, max_layer, min_band, max_band, min_state, max_state), 1)
                    
                    ! Energy Structure
                    do j = 1, num_found
                        if (energy(j) > maxval(energy_range(i, :))) then
                            energymap_cpu(kl(k), j) = maxval(energy_range(i, :))
                        else if (energy(j) < minval(energy_range(i, :))) then
                            energymap_cpu(kl(k), j) = minval(energy_range(i, :))
                        else
                            energymap_cpu(kl(k), j) = energy(j)
                        end if
                    end do
                    
                    ! Band Contribution
                    contribution_cpu(kl(k), :, :) = sum(hqtlib_find_band_contribution(energy(:num_found), weight(:, :num_found), &
                                                                                      num_bands, &
                                                                min_layer, max_layer, min_band, max_band, min_state, max_state), 1)
                end if
            end if
            if (cpu_is_master()) call main_status(k - k_start + 1, k_stop - k_start + 1, start_time)
        end do
        if (cpu_is_master()) then
            call main_clear_status()
            write(*, fmt = "(A)") "100%"
        end if
        
        den_cpu     = den_cpu * 1d-18
        heatmap_cpu = heatmap_cpu + 1d0
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
            
            ! Variation Meta Data
            variation_dir = export_create_dir(trim(path), "Variation "//export_to_string(i))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                "#      Num K Len      Num Layers      Num Path K Num Path Energy")
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/   num_k_length,     num_layers,       size(kp),     num_energy /))
            call export_vstack(trim(variation_dir)//"meta.dat", &
                "#     Broadening Crystal Len (a)   K Len [2pi*a]")
            call export_vstack(trim(variation_dir)//"meta.dat", &
                (/     broadening, crystal_length * 1d10, length_scale /))
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
              
            ! Variation Potential
            call export_hstack(trim(variation_dir)//"potential.dat", pot(i, :))
            call graph_basic(data_folder = trim(variation_dir)//"potential.dat", &
                             output_file = trim(variation_dir)//"potential", &
                             x_label = "Layer", &
                             y_label = "Potential [eV]")
            
            ! Variation Density
            call export_hstack(trim(variation_dir)//"density.dat", sum(sum(den, 2), 2))
            call export_hstack(trim(variation_dir)//"density.dat", transpose(sort_normalise(sum(den, 3))))
            call graph_colour(data_folder    = trim(variation_dir)//"density.dat", &
                              output_file    = trim(variation_dir)//"density", &
                              x_label        = "Layer", &
                              y_label        = "Carrier Density [nm^{-2}]", &
                              x_triangle     = 0.75d0, &
                              y_triangle     = 0.75d0, &
                              triangle_size  = 0.15d0, &
                              column_label_1 = "d_{xy}", &
                              column_label_2 = "d_{xz}", &
                              column_label_3 = "d_{yz}")
                              
            ! Variation Heatmap                  
            call export_hstack(trim(variation_dir)//"band structure.dat",  log(sum(sum(heatmap, 2), 2)))                  
            call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(variation_dir))
            
            ! Variation Energymap
            energymap_dir = export_create_dir(trim(variation_dir), "Energymap")
            do j = 1, maxval(max_num_states) - min_state + 1
                call export_hstack(trim(energymap_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", energymap(:, j))
                call export_hstack(trim(energymap_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                                   transpose(contribution(:, :, j)))
            end do
            call graph_colour_MGX(data_folder    = trim(energymap_dir), &
                                  output_file    = trim(energymap_dir), &
                                  x_label        = "K_{M Γ X} [π/a]", &
                                  y_label        = "E - E_{cbm} [ev]", &
                                  x_triangle     = 0.2d0, &
                                  y_triangle     = 0.2d0, &
                                  triangle_size  = 0.15d0, &
                                  column_label_1 = "d_{xy}", &
                                  column_label_2 = "d_{xz}", &
                                  column_label_3 = "d_{yz}", &
                                  k_scale        = length_scale, &
                                  min_y          = minval(energy_range(i, :)), &
                                  max_y          = maxval(energy_range(i, :)))
            
            do j = 1, max_num_states(i) - min_state + 1
                energy_state_dir = export_create_dir(trim(variation_dir)//"Energy States/", &
                    "Energy State "//export_to_string(j))
                    
                ! Variation Energy State Density
                call export_hstack(trim(energy_state_dir)//"density.dat", sum(den(:, :, j), 2))
                call export_hstack(trim(energy_state_dir)//"density.dat", transpose(sort_normalise(den(:, :, j))))
                call graph_colour(data_folder    = trim(energy_state_dir)//"density.dat", &
                                  output_file    = trim(energy_state_dir)//"density", &
                                  x_label        = "Layer", &
                                  y_label        = "Carrier Density [nm^{-2}]", &
                                  x_triangle     = 0.75d0, &
                                  y_triangle     = 0.75d0, &
                                  triangle_size  = 0.15d0, &
                                  column_label_1 = "d_{xy}", &
                                  column_label_2 = "d_{xz}", &
                                  column_label_3 = "d_{yz}")
                                  
                ! Variation Energy State Heatmap
                call export_hstack(trim(energy_state_dir)//"band structure"//".dat", log(sum(heatmap(:, :, j, :), 2)))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(energy_state_dir))
                
                ! Variation Energy State Energymap
                call graph_colour_MGX(data_folder    = trim(energymap_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                                      output_file    = trim(energy_state_dir)//"Energymap", &
                                      x_label        = "K_{M Γ X} [π/a]", &
                                      y_label        = "E - E_{cbm} [ev]", &
                                      x_triangle     = 0.2d0, &
                                      y_triangle     = 0.2d0, &
                                      triangle_size  = 0.15d0, &
                                      column_label_1 = "d_{xy}", &
                                      column_label_2 = "d_{xz}", &
                                      column_label_3 = "d_{yz}", &
                                      k_scale        = length_scale, &
                                      min_y          = minval(energy_range(i, :)), &
                                      max_y          = maxval(energy_range(i, :)))
                              
                if (cpu_is_master()) call main_status(j, max_state - min_state + 1, start_time)
            end do
            if (cpu_is_master()) then
                call main_clear_status()
                write(*, fmt = "(A)") "100%"
            end if
            
            ! Total Density
            call export_vstack(trim(path)//"Total_Density.dat", sum(den))
            call export_vstack(trim(path)//"Quantum_Well_Density.dat", sum(den(well_start:well_stop, :, :)))
            if (i > 2) then
                call graph_basic(data_folder = trim(path)//"Total_Density.dat", &
                                 output_file = trim(path)//"Total_Density", &
                                 x_label = x_label, &
                                 y_label = "Carrier Density [nm^{-2}]")
                call graph_basic(data_folder = trim(path)//"Quantum_Well_Density.dat", &
                                 output_file = trim(path)//"Quantum_Well_Density", &
                                 x_label = x_label, &
                                 y_label = "Carrier Density [nm^{-2}]")
            end if
            
            ! Symmetry Point Energymap
            do j = 1, max_state - min_state + 1
                ! Gamma
                call export_vstack(trim(gamma_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                    (/ variation_parameter_width(i), variation_parameter_depth(i), energymap((size(kp)+1) / 2, j), &
                    contribution((size(kp)+1) / 2, 1, j), & ! Light Band
                    dble(j) / dble(max_state - min_state + 1), & ! Seperate same Band Colours
                    (contribution((size(kp)+1) / 2, 2, j) + contribution((size(kp)+1) / 2, 3, j)) &
                    /)) ! Heavy Band
                ! M
                call export_vstack(trim(m_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                    (/ variation_parameter_width(i), variation_parameter_depth(i), energymap(1, j), &
                    contribution(1, 1, j), & ! Light Band
                    dble(j) / dble(max_state - min_state + 1), & ! Seperate same Band Colours
                    (contribution(1, 2, j) + contribution(1, 3, j)) /)) ! Heavy Band
                ! X
                call export_vstack(trim(x_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                    (/ variation_parameter_width(i), variation_parameter_depth(i), energymap(size(kp), j), &
                    contribution(size(kp), 1, j), & ! Light Band
                    dble(j) / dble(max_state - min_state + 1), & ! Seperate same Band Colours
                    (contribution(size(kp), 2, j) + contribution(size(kp), 3, j)) /)) ! Heavy Band
                ! GX
                call export_vstack(trim(gx_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                    (/ variation_parameter_width(i), variation_parameter_depth(i), energymap((size(kp)+1) * 9/16, j), &
                    contribution((size(kp)+1) * 9/16, 1, j), & ! Light Band
                    dble(j) / dble(max_state - min_state + 1), & ! Seperate same Band Colours
                    (contribution((size(kp)+1) * 9/16, 2, j) + contribution((size(kp)+1) * 9/16, 3, j)) /)) ! Heavy Band
                ! GM
                call export_vstack(trim(gm_dir)//"Energy_level_"//export_to_string(j, "0")//".dat", &
                    (/ variation_parameter_width(i), variation_parameter_depth(i), energymap((size(kp)+1) * 7/16, j), &
                    contribution((size(kp)+1) * 7/16, 1, j), & ! Light Band
                    dble(j) / dble(max_state - min_state + 1), & ! Seperate same Band Colours
                    (contribution((size(kp)+1) * 7/16, 2, j) + contribution((size(kp)+1) * 7/16, 3, j)) /)) ! Heavy Band
            end do
            if (i == num_variation**2) then
                ! Gamma
                call graph_colour_3d(trim(gamma_dir), trim(path)//"Gamma_Fence", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "fence", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gamma_dir), trim(path)//"Gamma_Surface", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "surface", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gamma_dir), trim(path)//"Gamma_Seperate", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "none", &
                                     plot_type          = "seperate", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0, &
                                     has_axis           = .false.)
                ! M
                call graph_colour_3d(trim(m_dir), trim(path)//"M_Fence", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "fence", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(m_dir), trim(path)//"M_Surface", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "surface", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(m_dir), trim(path)//"M_Seperate", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "none", &
                                     plot_type          = "seperate", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0, &
                                     has_axis           = .false.)
               ! X
               call graph_colour_3d(trim(x_dir), trim(path)//"X_Fence", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "fence", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(x_dir), trim(path)//"X_Surface", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "surface", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(x_dir), trim(path)//"X_Seperate", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "none", &
                                     plot_type          = "seperate", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0, &
                                     has_axis           = .false.)
                ! GX
                call graph_colour_3d(trim(gx_dir), trim(path)//"GX_Fence", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "fence", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gx_dir), trim(path)//"GX_Surface", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "surface", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gx_dir), trim(path)//"GX_Seperate", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "none", &
                                     plot_type          = "seperate", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0, &
                                     has_axis           = .false.)
                ! GM
                call graph_colour_3d(trim(gm_dir), trim(path)//"GM_Fence", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "fence", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gm_dir), trim(path)//"GM_Surface", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "box", &
                                     plot_type          = "surface", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0)
                call graph_colour_3d(trim(gm_dir), trim(path)//"GM_Seperate", &
                                     x_label            = "Width [Layer]", &
                                     y_label            = "Depth [eV]", &
                                     z_label            = "E - E_{cbm} [ev]", &
                                     z_max              = maxval(energy_range), &
                                     z_min              = minval(energy_range), &
                                     colour_box_x       = 0.9d0, &
                                     colour_box_y       = 0.5d0, &
                                     colour_box_size    = 0.1d0, &
                                     colour_box_type    = "none", &
                                     plot_type          = "seperate", &
                                     colour_box_label_1 = "Light", &
                                     colour_box_label_2 = "Heavy", &
                                     transparency       = 0.5d0, &
                                     pitch              = 75d0, &
                                     has_axis           = .false.)
            end if
            
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
