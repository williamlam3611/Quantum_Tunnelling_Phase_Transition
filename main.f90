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
    use spectral
    use density
    use export
    use graph
    use mpi
    implicit none
    
    integer,        parameter   :: num_variation      = 1 !25
    integer,        parameter   :: num_layers         = 20 !150
    
    integer                     :: well_1_start       = 1 
    integer                     :: well_1_stop        = 1 
    real*8                      :: well_1_start_depth = 0d0
    real*8                      :: well_1_stop_depth  = 0d0
    
    integer                     :: well_2_start       = 1!6
    integer                     :: well_2_stop        = 1!45!20
    real*8                      :: well_2_start_depth = 0d0 !-0.5d0
    real*8                      :: well_2_stop_depth  = 0d0

    real*8                      :: c_parameter = 0.09 !0.01 +0.5
    
    ! for varying depth
    real*8                      :: b_parameter = 0.645!0.6
    
    ! for varying width 
    real*8                      :: a_parameter = 2.5
        
    
    integer,        parameter   :: num_k_length       = 256
    integer,        parameter   :: num_energy         = 256
    integer,        parameter   :: output_bands(*)    = (/ 1, 2, 3 /)
    integer,        parameter   :: output_states(*)   = (/ 1, 2, 3 /)
    real*8                      :: energy_min         = -1.3d0
    real*8                      :: energy_max         = 0d0
    real*8,         parameter   :: length_scale       = 1d0
    real*8,         parameter   :: broadening         = 0.0025d0
    
    integer                     :: x_offset           = 0
    character(*),   parameter   :: x_label            = "Surface Layer Doping [-0.1eV]"
    
    logical,        parameter   :: auto_energy_min    = .true.
    logical,        parameter   :: auto_energy_max    = .true.
    
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
    
    real*8,         allocatable :: temp_heatmap_1(:, :, :, :), temp_heatmap_2(:, :, :, :), &
                                   heatmap_cpu(:, :, :, :), heatmap_total_cpu(:, :), &
                                   heatmap(:, :, :, :), heatmap_total(:, :)
    real*8,         allocatable :: temp_den_1(:, :, :), temp_den_2(:, :, :), &
                                   den_cpu(:, :, :), den_total_cpu(:), &
                                   den(:, :, :), den_total(:)
    
    character(256)              :: path, variation_dir, band_dir, energy_state_dir
    
    integer                     :: i, j, k, l, m, e, n
    
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
        !well_2_stop = i
        !variation_parameter_list(i) = well_2_stop - 1
        
        
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
        !call potential_add_well(pot(i, :), well_2_start, well_2_stop, well_2_start_depth, well_2_stop_depth)
        
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
    num_states = num_bands * num_layers
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
    
    ! Determine Maximum number of bands
    energy_min = energy_min + cbm
    energy_max = energy_max + cbm
    do i = 1, num_variation
        if (auto_energy_min) then
            energy_min = hqtlib_find_energy_min(hr, hrw, pot(i, :), broadening)
        end if
        if (auto_energy_max) then
            energy_max = hqtlib_find_energy_max(hr, hrw, pot(i, :), length_scale, broadening)
        end if
        max_num_states(i) = hqtlib_find_max_num_states(hr, hrw, pot(i, :), energy_min, energy_max)
        if (max_num_states(i) == 0) then
            max_num_states(i) = 1
        end if
        energy_range(i, :) = route_range_double(energy_min, energy_max, num_energy) - cbm
    end do
    
    ! Allocate Dynamic Arrays
    allocate(energy(num_states))
    allocate(weight(num_states, num_states))
    allocate(heatmap_cpu(      size(kp),   size(output_bands), size(output_states), num_energy))
    allocate(heatmap_total_cpu(size(kp),   num_energy))
    allocate(den_cpu(          num_layers, size(output_bands), size(output_states)))
    allocate(den_total_cpu(    num_layers))
    if (cpu_is_master()) then
        allocate(heatmap(      size(kp),   size(output_bands), size(output_states), num_energy))
        allocate(heatmap_total(size(kp),   num_energy))
        allocate(den(          num_layers, size(output_bands), size(output_states)))
        allocate(den_total(    num_layers))
        
    end if
    
    ! Set up Output Folder
    if (cpu_is_master()) then
        path = export_create_dir(out_dir, "Run ")
        call export_data(trim(path)//"density.dat", &
            "#          Total  Quantum Well 1  Quantum Well 2             Sum"// &
            "       Avg Total      Avg Well 1      Avg Well 2    Avg Well Sum")
        call export_data(trim(path)//"variation_parameter_list.dat", variation_parameter_list )
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
        allocate(gamma_energy_list(num_variation, num_layers* num_bands))
    end if 
    
    allocate(gamma_energy(num_layers* num_bands))
    
    do i = 1, num_variation
        den_cpu            = 0d0
        den_total_cpu      = 0d0
        heatmap_cpu        = 0d0
        heatmap_total_cpu  = 0d0
        if (cpu_is_master()) then
            den            = 0d0
            den_total      = 0d0
            heatmap        = 0d0
            heatmap_total  = 0d0
        end if
        
        ! determine gamma point 
        gamma_k_pos  = size(kw)
        gamma_energy = 0d0
        ! Spectral into Density and Heatmap
        do k = k_start, k_stop
            call hqtlib_find_energy_and_weight(energy, weight, num_found, kx(k), ky(k), hr, hrw, pot(i, :), 1, max_num_states(i))
            energy(:num_found) = energy(:num_found) - cbm
            
            if (k == size(kw)) then            
                gamma_energy = energy               
                call cpu_send_double(gamma_energy, size(gamma_energy), 0, i)
                cpu_from = cpu_get_id()
                
            end if 
            
            call hqtlib_find_density(temp_den_1, energy(:num_found), weight(:, :num_found), energy_range(i, :), broadening, &
                                     num_bands, crystal_length, sum(kw)) 
            den_total_cpu = den_total_cpu + sum(sum(temp_den_1, 2), 2) * kw(k)
            call hqtlib_find_density(temp_den_2, energy(:num_found), weight(:, :num_found), energy_range(i, :), broadening, &
                                     num_bands, crystal_length, sum(kw), bands = output_bands, energy_levels = output_states) 
            den_cpu       = den_cpu + temp_den_2 * kw(k)
            if (kl(k) .ne. -1) then
                call hqtlib_find_spectral_distribution(temp_heatmap_1, energy(:num_found), weight(:, :num_found), &
                                                       energy_range(i, :), broadening, num_bands)
                heatmap_total_cpu(kl(k), :) = heatmap_total_cpu(kl(k), :) + sum(sum(sum(temp_heatmap_1, 1), 1), 1)
                call hqtlib_find_spectral_distribution(temp_heatmap_2, energy(:num_found), weight(:, :num_found), &
                                                       energy_range(i, :), broadening, num_bands, &
                                                       bands = output_bands, energy_levels = output_states)
                heatmap_cpu(kl(k), :, :, :) = heatmap_cpu(kl(k), :, :, :) + sum(temp_heatmap_2, 1)
            end if
        end do
        
        if (cpu_is_master()) then
            call cpu_recv_double(gamma_energy, size(gamma_energy),  MPI_ANY_SOURCE , i)
            gamma_energy_list(i,:) = gamma_energy
        end if 
        
        den_cpu        = den_cpu * 1d-18
        den_total_cpu  = den_total_cpu * 1d-18
        call cpu_sum(den_cpu, den)
        call cpu_sum(den_total_cpu, den_total)
        call cpu_sum(heatmap_cpu, heatmap)
        call cpu_sum(heatmap_total_cpu, heatmap_total)
        ! No more multithreading for rest of cycle
        if (cpu_is_master()) then
            ! Export Meta Data
            variation_dir = export_create_dir(trim(path), "Variation "//export_to_string(i))
            call export_meta_data(trim(variation_dir)//"meta.dat", &
                num_k_length, num_layers, size(kp), num_energy, &
                broadening, crystal_length, length_scale, &
                minval(pot(i, :)), maxval(pot(i, :)), minval(den_total), maxval(den_total), &
                minval(energy_range(i, :)), maxval(energy_range(i, :)), minval(log(heatmap_total)), maxval(log(heatmap_total)))
                
            ! Export and Plot Variation Data
            call export_data(trim(variation_dir)//"potential.dat",      pot(i, :))
            call export_data(trim(variation_dir)//"gamma_energy.dat",   gamma_energy)
            call export_data(trim(variation_dir)//"density.dat",        den_total)
            call export_data(trim(variation_dir)//"band structure.dat", transpose(log(heatmap_total)))
            call graph_basic_plot("potential", "potential", &
                1, "Layer", "Potential [eV]", 1, trim(variation_dir))
            call graph_basic_plot("density", "density", &
                1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(variation_dir))
            call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(variation_dir))
            
            ! Export and Plot Variation Band Data
            do j = 1, size(output_bands)
                band_dir = export_create_dir(trim(variation_dir)//"Bands/", "Band "//export_to_string(output_bands(j)))
                call export_data(trim(band_dir)//"density"//".dat", &
                    sum(den(:, j, :), 2))
                call export_data(trim(band_dir)//"band structure"//".dat", &
                    transpose(log(sum(heatmap(:, j, :, :), 2))))
                call graph_basic_plot("density", "density", &
                    1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(band_dir))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(band_dir))
            end do
            
            ! Export and Plot Variation Energy State Data
            do j = 1, size(output_states)
                energy_state_dir = export_create_dir(trim(variation_dir)//"Energy States/", &
                    "Energy State "//export_to_string(output_states(j)))
                call export_data(trim(energy_state_dir)//"density"//".dat", &
                    sum(den(:, :, j), 2))
                call export_data(trim(energy_state_dir)//"band structure"//".dat", &
                    transpose(log(sum(heatmap(:, :, j, :), 2))))
                call graph_basic_plot("density", "density", &
                    1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(energy_state_dir))
                call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(energy_state_dir))
                do n = 1, size(output_bands)
                    band_dir = export_create_dir(trim(energy_state_dir)//"Bands/", "Band " &
                        //export_to_string(output_bands(n)))
                    call export_data(trim(band_dir)//"density"//".dat", &
                        den(:, n, j))
                    call export_data(trim(band_dir)//"band structure"//".dat", &
                        transpose(log(heatmap(:, n, j, :))))
                    call graph_basic_plot("density", "density", &
                        1, "Layer", "Carrier Density [nm^{-2}]", 1, trim(band_dir))
                    call graph_heatmap_plot("band structure", "band structure", trim(variation_dir)//"meta", trim(band_dir))
                end do
            end do
            
            ! Export and Plot Data
            call export_data(trim(path)//"density.dat", reshape((/ &
                sum(den_total), &
                sum(den_total(well_1_start:well_1_stop)), &
                sum(den_total(well_2_start:well_2_stop)), &
                sum(den_total(well_1_start:well_1_stop)) + &
                    sum(den_total(well_2_start:well_2_stop)), &
                sum(den_total) / size(den_total), &
                sum(den_total(well_1_start:well_1_stop)) / (abs(well_1_stop - well_1_start) + 1), &
                sum(den_total(well_2_start:well_2_stop)) / (abs(well_2_stop - well_2_start) + 1), &
                sum(den_total(well_1_start:well_1_stop)) + &
                    sum(den_total(well_2_start:well_2_stop)) / (abs(well_1_stop + well_2_stop - well_1_start - well_2_start) + 2) &
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
                call export_data(trim(path)//"gamma_energy_list.dat", gamma_energy_list(1:i,1:11))

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
                    
                    call export_data(trim(path)//"gamma_energy_list.dat", gamma_energy_list(1:i,1:11))              
                
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
    

end program main
