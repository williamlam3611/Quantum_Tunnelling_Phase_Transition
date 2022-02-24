set print "-";

is_folder = 0;
if (!exists("data_folder"))               data_folder       = GPVAL_PWD;
if (strstrt(substr(data_folder, 2, strlen(data_folder)), ".") == 0) is_folder = 1;
if (!exists("output_file") && is_folder)  output_file       = data_folder;
if (!exists("output_file") && !is_folder) output_file       = data_folder[1:strstrt(data_folder, ".") - 1];
if (!exists("output_file_type"))          output_file_type  = "svg";
if (!exists("y_resolution"))              y_resolution      = 1680;
if (!exists("x_resolution"))              x_resolution      = 1050;
if (!exists("font"))                      font              = "Helvetica";
if (!exists("font_size"))                 font_size         = "16";
if (!exists("background_colour"))         background_colour = "white";
if (!exists("x_label"))                   x_label           = "X";
if (!exists("y_label"))                   y_label           = "Y";
if (!exists("x_triangle"))                x_triangle        = 0.1;
if (!exists("y_triangle"))                y_triangle        = 0.2;
if (!exists("triangle_size"))             triangle_size     = 0.2;
if (!exists("column_label_1"))            column_label_1    = "red";
if (!exists("column_label_2"))            column_label_2    = "blue";
if (!exists("column_label_3"))            column_label_3    = "green";
if (!exists("k_scale"))                   k_scale           = 1.0;
if (!exists("min_y"))                     min_y             = -1.0;
if (!exists("max_y"))                     max_y             = 1.0;

max_col = 0;
max_row = 0;
probe_file = data_folder;
if (is_folder) {
    pwd = GPVAL_PWD;
    cd data_folder;
    data_files = system("ls -1 *.dat");
    cd pwd;
    
    probe_file = data_folder.data_files[1:strstrt(data_files, "\n") - 1];
};
stats probe_file nooutput;
max_col = STATS_columns;
max_row = STATS_records;

if (substr(output_file, strlen(output_file), strlen(output_file)) eq "/") {
    output_file = substr(output_file, 1, strlen(output_file) - 1);
}

set terminal output_file_type size y_resolution, x_resolution font font.", ".font_size background rgb background_colour;
set output output_file.".".output_file_type;
set key off;
set tics scale 0;

rgb(r,g,b) = 65536 * int(r * 256) + 256 * int(g * 256) + int(b * 256)

set multiplot layout 1, 3;

unset xlabel;
unset ylabel;
unset border;
unset xtics;
unset ytics;
set size triangle_size, triangle_size;
set origin x_triangle, y_triangle;
set label 1 column_label_3 at screen x_triangle, y_triangle centre textcolor rgb "black";
set label 2 column_label_2 at screen x_triangle + triangle_size, y_triangle centre textcolor rgb "black";
set label 3 column_label_1 at screen x_triangle + triangle_size / 2, y_triangle + triangle_size centre textcolor rgb "black";
plot 'colour_triangle.png' binary filetype=png with rgbalpha;
set xlabel;
set ylabel y_label;
set border;

set label 3 x_label at screen (2 - sqrt(2)), 0.02 centre;
set xtics (sprintf("%3.2f", sqrt(2) * k_scale) 0, \
    sprintf("%3.2f", 1.00 * k_scale) (max_row / 2) - ((sqrt(2) / 2) * (max_row / 2)), \
    sprintf("%3.2f", 0.50 * k_scale) (max_row / 2) - ((sqrt(2) / 4) * (max_row / 2)), \
    "0.00" (max_row / 2), \
    sprintf("%3.2f", 0.50 * k_scale) ((3 * max_row) / 4), \
    sprintf("%3.2f", 1.00 * k_scale) (max_row - 2));

set ytics
set size (2 - sqrt(2)), (1 - 0.03);
set origin 0, 0.03
set rmargin 0
set xrange [0:(max_row / 2)]
set yrange [min_y:max_y]
if (is_folder) {
    pwd = GPVAL_PWD;
    cd data_folder;
    plot for [data in data_files] data using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
    cd pwd;
} else {
    plot data_folder using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
};

set size (sqrt(2) - 1), (1 - 0.03);
set lmargin 0
unset rmargin
set origin (2 - sqrt(2)), 0.03;
unset ylabel
set ytics font ",0"
set xrange [((max_row / 2)):max_row]
if (is_folder) {
    pwd = GPVAL_PWD;
    cd data_folder;
    plot for [data in data_files] data using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
    cd pwd;
} else {
    plot data_folder using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
};

unset multiplot
