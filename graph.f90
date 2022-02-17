module graph
    use export
    implicit none

    public :: graph_basic_plot, graph_heatmap_plot, graph_multiple_plot, graph_multiple_coloured, graph_multiple_coloured_MGX

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
    
end module graph
