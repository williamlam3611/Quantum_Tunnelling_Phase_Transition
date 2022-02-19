module graph
    use export
    implicit none

    public :: graph_basic_plot, graph_heatmap_plot, graph_multiple_plot, graph_multiple_coloured, graph_multiple_coloured_MGX, &
              graph_colour

contains
    subroutine graph_basic_plot(data_file, output_file, data_column, x_label, y_label, x_offset, path)
        character(*), intent(in) :: data_file, output_file, x_label, y_label, path
        integer,      intent(in) :: data_column, x_offset
        call execute_command_line("gnuplot -e """ &
            //"data_file='"//data_file//"';" &
            //"output_file='"//output_file//"';" &
            //"data_column="//export_to_string(data_column)//";" &
            //"x_label='"//x_label//"';" &
            //"y_label='"//y_label//"';" &
            //"x_offset="//export_to_string(x_offset)//";" &
            //"path='"//trim(path)//"'" &
            //""" basic_plot.p")
    end subroutine graph_basic_plot

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
    
    subroutine graph_multiple_coloured_MGX(path, output_file, meta_file, x_label, y_label)
        character(*), intent(in) :: path, output_file, meta_file, x_label, y_label
        call execute_command_line("gnuplot -e """ &
            //"path='"//path//"';" &
            //"output_file='"//output_file//"';" &
            //"meta_file='"//meta_file//"';" &
            //"x_label='"//x_label//"';" &
            //"y_label='"//y_label//"'" &
            //""" colour_plot_MGX.p")
    end subroutine graph_multiple_coloured_MGX
    
    subroutine graph_colour(data_folder, output_file, output_file_type, y_resolution, x_resolution, font, font_size, &
                            background_colour, x_label, y_label, x_triangle, y_triangle, triangle_size)
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
        call graph_buffer(buffer, """ colour_plot.p")
        call execute_command_line(graph_to_string(buffer))
        
    end subroutine graph_colour
    
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
