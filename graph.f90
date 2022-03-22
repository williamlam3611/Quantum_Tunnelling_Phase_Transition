module graph
    use export
    implicit none

    public :: graph_basic, graph_heatmap_plot, graph_multiple_plot, graph_multiple_coloured, &
              graph_colour, graph_colour_MGX, graph_colour_3D
              
contains
    subroutine graph_basic(data_folder, output_file, output_file_type, y_resolution, x_resolution, font, font_size, &
                           background_colour, x_label, y_label)
        character(*), optional :: data_folder
        character(*), optional :: output_file
        character(*), optional :: output_file_type
        integer,      optional :: y_resolution
        integer,      optional :: x_resolution
        character(*), optional :: font
        integer,      optional :: font_size
        character(*), optional :: background_colour
        character(*), optional :: x_label
        character(*), optional :: y_label
        character, allocatable :: buffer(:)
        
        call graph_buffer(buffer, "gnuplot -e """)
        if (present(data_folder)) call graph_buffer(buffer, "data_folder = '"//data_folder//"';")
        if (present(output_file)) call graph_buffer(buffer, "output_file = '"//output_file//"';")
        if (present(output_file_type)) call graph_buffer(buffer, "output_file_type = '"//output_file_type//"';")
        if (present(y_resolution)) call graph_buffer(buffer, "y_resolution = '"//export_to_string(y_resolution)//"';")
        if (present(x_resolution)) call graph_buffer(buffer, "x_resolution = '"//export_to_string(x_resolution)//"';")
        if (present(font)) call graph_buffer(buffer, "font = '"//font//"';")
        if (present(font_size)) call graph_buffer(buffer, "font_size = '"//export_to_string(font_size)//"';")
        if (present(background_colour)) call graph_buffer(buffer, "background_colour = '"//background_colour//"';")
        if (present(x_label)) call graph_buffer(buffer, "x_label = '"//x_label//"';")
        if (present(y_label)) call graph_buffer(buffer, "y_label = '"//y_label//"';")
        call graph_buffer(buffer, """ basic_plot.p")
        call execute_command_line(graph_to_string(buffer))
        
    end subroutine graph_basic

    subroutine graph_multiple_plot(data_file, output_file, x_label, y_label, x_offset, path)
        character(*), intent(in) :: data_file, output_file, x_label, y_label, path
        integer,      intent(in) :: x_offset ! data_column, 
        call execute_command_line("gnuplot -e """ &
            //"data_file='"//data_file//"';" &
            //"output_file='"//output_file//"';" &
            !//"data_column="//export_to_string(data_column)//";" &
            //"x_label='"//x_label//"';" &
            //"y_label='"//y_label//"';" &
            //"x_offset="//export_to_string(x_offset)//";" &
            //"path='"//trim(path)//"'" &
            //""" gamma_varation_plot.p")
    end subroutine graph_multiple_plot
    
    subroutine graph_heatmap_plot(data_file, output_file, meta_file, path)
        character(*), intent(in) :: data_file, output_file, meta_file, path
        call execute_command_line("gnuplot -e """ &
            //"data_file='"//data_file//"';" &
            //"output_file='"//output_file//"';" &
            //"meta_file='"//meta_file//"';" &
            //"path='"//path//"'" &
            //""" basic_heatmap.p")
    end subroutine graph_heatmap_plot
    
    subroutine graph_multiple_coloured(path, output_file, x_label, y_label)
        character(*), intent(in) :: path, output_file, x_label, y_label
        call execute_command_line("gnuplot -e """ &
            //"path='"//path//"';" &
            //"output_file='"//output_file//"';" &
            //"x_label='"//x_label//"';" &
            //"y_label='"//y_label//"'" &
            //""" colour_plot.p")
    end subroutine graph_multiple_coloured
    
    subroutine graph_colour(data_folder, output_file, output_file_type, y_resolution, x_resolution, font, font_size, &
                            background_colour, x_label, y_label, x_triangle, y_triangle, triangle_size, &
                            column_label_1, column_label_2, column_label_3)
        character(*), optional :: data_folder
        character(*), optional :: output_file
        character(*), optional :: output_file_type
        integer,      optional :: y_resolution
        integer,      optional :: x_resolution
        character(*), optional :: font
        integer,      optional :: font_size
        character(*), optional :: background_colour
        character(*), optional :: x_label
        character(*), optional :: y_label
        real*8,       optional :: x_triangle
        real*8,       optional :: y_triangle
        real*8,       optional :: triangle_size
        character(*), optional :: column_label_1
        character(*), optional :: column_label_2
        character(*), optional :: column_label_3
        character, allocatable :: buffer(:)
        
        call graph_buffer(buffer, "gnuplot -e """)
        if (present(data_folder)) call graph_buffer(buffer, "data_folder = '"//data_folder//"';")
        if (present(output_file)) call graph_buffer(buffer, "output_file = '"//output_file//"';")
        if (present(output_file_type)) call graph_buffer(buffer, "output_file_type = '"//output_file_type//"';")
        if (present(y_resolution)) call graph_buffer(buffer, "y_resolution = '"//export_to_string(y_resolution)//"';")
        if (present(x_resolution)) call graph_buffer(buffer, "x_resolution = '"//export_to_string(x_resolution)//"';")
        if (present(font)) call graph_buffer(buffer, "font = '"//font//"';")
        if (present(font_size)) call graph_buffer(buffer, "font_size = '"//export_to_string(font_size)//"';")
        if (present(background_colour)) call graph_buffer(buffer, "background_colour = '"//background_colour//"';")
        if (present(x_label)) call graph_buffer(buffer, "x_label = '"//x_label//"';")
        if (present(y_label)) call graph_buffer(buffer, "y_label = '"//y_label//"';")
        if (present(x_triangle)) call graph_buffer(buffer, "x_triangle = '"//export_to_string(x_triangle)//"';")
        if (present(y_triangle)) call graph_buffer(buffer, "y_triangle = '"//export_to_string(y_triangle)//"';")
        if (present(triangle_size)) call graph_buffer(buffer, "triangle_size = '"//export_to_string(triangle_size)//"';")
        if (present(column_label_1)) call graph_buffer(buffer, "column_label_1 = '"//column_label_1//"';")
        if (present(column_label_2)) call graph_buffer(buffer, "column_label_2 = '"//column_label_2//"';")
        if (present(column_label_3)) call graph_buffer(buffer, "column_label_3 = '"//column_label_3//"';")
        call graph_buffer(buffer, """ colour_plot.p")
        call execute_command_line(graph_to_string(buffer))
        
    end subroutine graph_colour

    subroutine graph_colour_3d(data_folder, output_file, output_file_type, y_resolution, x_resolution, font, font_size, &
                            background_colour, x_label, y_label, z_label, colour_box_x, colour_box_y, colour_box_size, &
                            colour_box_label_1, colour_box_label_2, colour_box_label_3, colour_box_type, plot_type, &
                            transparency, pitch, yaw, x_min, y_min, z_min, x_max, y_max, z_max, has_axis)
        character(*), optional :: data_folder
        character(*), optional :: output_file
        character(*), optional :: output_file_type
        integer,      optional :: y_resolution
        integer,      optional :: x_resolution
        character(*), optional :: font
        integer,      optional :: font_size
        character(*), optional :: background_colour
        character(*), optional :: x_label
        character(*), optional :: y_label
        character(*), optional :: z_label
        real*8,       optional :: colour_box_x
        real*8,       optional :: colour_box_y
        real*8,       optional :: colour_box_size
        character(*), optional :: colour_box_label_1
        character(*), optional :: colour_box_label_2
        character(*), optional :: colour_box_label_3
        character(*), optional :: colour_box_type
        character(*), optional :: plot_type
        real*8,       optional :: transparency
        real*8,       optional :: pitch
        real*8,       optional :: yaw
        real*8,       optional :: x_min
        real*8,       optional :: y_min
        real*8,       optional :: z_min
        real*8,       optional :: x_max
        real*8,       optional :: y_max
        real*8,       optional :: z_max
        logical,      optional :: has_axis
        character, allocatable :: buffer(:)
        
        call graph_buffer(buffer, "gnuplot -e """)
        if (present(data_folder)) call graph_buffer(buffer, "data_folder = '"//data_folder//"';")
        if (present(output_file)) call graph_buffer(buffer, "output_file = '"//output_file//"';")
        if (present(output_file_type)) call graph_buffer(buffer, "output_file_type = '"//output_file_type//"';")
        if (present(y_resolution)) call graph_buffer(buffer, "y_resolution = '"//export_to_string(y_resolution)//"';")
        if (present(x_resolution)) call graph_buffer(buffer, "x_resolution = '"//export_to_string(x_resolution)//"';")
        if (present(font)) call graph_buffer(buffer, "font = '"//font//"';")
        if (present(font_size)) call graph_buffer(buffer, "font_size = '"//export_to_string(font_size)//"';")
        if (present(background_colour)) call graph_buffer(buffer, "background_colour = '"//background_colour//"';")
        if (present(x_label)) call graph_buffer(buffer, "x_label = '"//x_label//"';")
        if (present(y_label)) call graph_buffer(buffer, "y_label = '"//y_label//"';")
        if (present(z_label)) call graph_buffer(buffer, "z_label = '"//z_label//"';")
        if (present(colour_box_x)) call graph_buffer(buffer, "colour_box_x = '"//export_to_string(colour_box_x)//"';")
        if (present(colour_box_y)) call graph_buffer(buffer, "colour_box_y = '"//export_to_string(colour_box_y)//"';")
        if (present(colour_box_size)) call graph_buffer(buffer, "colour_box_size = '"//export_to_string(colour_box_size)//"';")
        if (present(colour_box_label_1)) call graph_buffer(buffer, "colour_box_label_1 = '"//colour_box_label_1//"';")
        if (present(colour_box_label_2)) call graph_buffer(buffer, "colour_box_label_2 = '"//colour_box_label_2//"';")
        if (present(colour_box_label_3)) call graph_buffer(buffer, "colour_box_label_3 = '"//colour_box_label_3//"';")
        if (present(colour_box_type)) call graph_buffer(buffer, "colour_box_type = '"//colour_box_type//"';")
        if (present(plot_type)) call graph_buffer(buffer, "plot_type = '"//plot_type//"';")
        if (present(transparency)) call graph_buffer(buffer, "transparency = '"//export_to_string(transparency)//"';")
        if (present(pitch)) call graph_buffer(buffer, "pitch = '"//export_to_string(pitch)//"';")
        if (present(yaw)) call graph_buffer(buffer, "yaw = '"//export_to_string(yaw)//"';")
        if (present(x_min)) call graph_buffer(buffer, "x_min = '"//export_to_string(x_min)//"';")
        if (present(y_min)) call graph_buffer(buffer, "y_min = '"//export_to_string(y_min)//"';")
        if (present(z_min)) call graph_buffer(buffer, "z_min = '"//export_to_string(z_min)//"';")
        if (present(x_max)) call graph_buffer(buffer, "x_max = '"//export_to_string(x_max)//"';")
        if (present(y_max)) call graph_buffer(buffer, "y_max = '"//export_to_string(y_max)//"';")
        if (present(z_max)) call graph_buffer(buffer, "z_max = '"//export_to_string(z_max)//"';")
        if (present(has_axis)) call graph_buffer(buffer, "has_axis = '"//export_to_string(has_axis)//"';")
        call graph_buffer(buffer, """ colour_plot_3d.p")
        
        call execute_command_line(graph_to_string(buffer))
        
    end subroutine graph_colour_3D

    
    subroutine graph_colour_MGX(data_folder, output_file, output_file_type, y_resolution, x_resolution, font, font_size, &
                            background_colour, x_label, y_label, x_triangle, y_triangle, triangle_size, &
                            column_label_1, column_label_2, column_label_3, k_scale, min_y, max_y)
        character(*), optional :: data_folder
        character(*), optional :: output_file
        character(*), optional :: output_file_type
        integer,      optional :: y_resolution
        integer,      optional :: x_resolution
        character(*), optional :: font
        integer,      optional :: font_size
        character(*), optional :: background_colour
        character(*), optional :: x_label
        character(*), optional :: y_label
        real*8,       optional :: x_triangle
        real*8,       optional :: y_triangle
        real*8,       optional :: triangle_size
        character(*), optional :: column_label_1
        character(*), optional :: column_label_2
        character(*), optional :: column_label_3
        real*8,       optional :: k_scale
        real*8,       optional :: min_y
        real*8,       optional :: max_y
        character, allocatable :: buffer(:)
        call graph_buffer(buffer, "gnuplot -e """)
        if (present(data_folder)) call graph_buffer(buffer, "data_folder = '"//data_folder//"';")
        if (present(output_file)) call graph_buffer(buffer, "output_file = '"//output_file//"';")
        if (present(output_file_type)) call graph_buffer(buffer, "output_file_type = '"//output_file_type//"';")
        if (present(y_resolution)) call graph_buffer(buffer, "y_resolution = '"//export_to_string(y_resolution)//"';")
        if (present(x_resolution)) call graph_buffer(buffer, "x_resolution = '"//export_to_string(x_resolution)//"';")
        if (present(font)) call graph_buffer(buffer, "font = '"//font//"';")
        if (present(font_size)) call graph_buffer(buffer, "font_size = '"//export_to_string(font_size)//"';")
        if (present(background_colour)) call graph_buffer(buffer, "background_colour = '"//background_colour//"';")
        if (present(x_label)) call graph_buffer(buffer, "x_label = '"//x_label//"';")
        if (present(y_label)) call graph_buffer(buffer, "y_label = '"//y_label//"';")
        if (present(x_triangle)) call graph_buffer(buffer, "x_triangle = '"//export_to_string(x_triangle)//"';")
        if (present(y_triangle)) call graph_buffer(buffer, "y_triangle = '"//export_to_string(y_triangle)//"';")
        if (present(triangle_size)) call graph_buffer(buffer, "triangle_size = '"//export_to_string(triangle_size)//"';")
        if (present(column_label_1)) call graph_buffer(buffer, "column_label_1 = '"//column_label_1//"';")
        if (present(column_label_2)) call graph_buffer(buffer, "column_label_2 = '"//column_label_2//"';")
        if (present(column_label_3)) call graph_buffer(buffer, "column_label_3 = '"//column_label_3//"';")
        if (present(k_scale)) call graph_buffer(buffer, "k_scale = '"//export_to_string(k_scale)//"';")
        if (present(min_y)) call graph_buffer(buffer, "min_y = '"//export_to_string(min_y)//"';")
        if (present(max_y)) call graph_buffer(buffer, "max_y = '"//export_to_string(max_y)//"';")
        call graph_buffer(buffer, """ colour_plot_MGX.p")
        call execute_command_line(graph_to_string(buffer))
        
    end subroutine graph_colour_MGX
    
    subroutine graph_buffer(buffer, string)
        character, allocatable, intent(inout) :: buffer(:)
        character(*), intent(in)              :: string
        character, allocatable                :: temp_buffer(:)
        integer                               :: i
        if (allocated(buffer)) then
            temp_buffer = buffer
            deallocate(buffer)
            allocate(buffer(size(temp_buffer) + len(string)))
            buffer(1:size(temp_buffer)) = temp_buffer
            do i = 1, len(string)
                buffer(size(temp_buffer) + i) = string(i:i)     
            end do
        else
            allocate(buffer(len(string)))
            do i = 1, len(string)
                buffer(i) = string(i:i)     
            end do           
        end if
    
    end subroutine graph_buffer
    
    function graph_to_string(value) result(string)
        character, allocatable, intent(in) :: value(:)
        character(size(value))             :: string
        integer                            :: i
        do i = 1, size(value)
            string(i:i) = value(i)
        end do
    
    end function graph_to_string
    
end module graph
