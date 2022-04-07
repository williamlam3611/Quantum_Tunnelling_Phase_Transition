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

if (is_folder) {
    pwd = GPVAL_PWD;
    cd data_folder;
    data_files = system("ls -1 *.dat");
    cd pwd;
};

if (substr(output_file, strlen(output_file), strlen(output_file)) eq "/") {
    output_file = substr(output_file, 1, strlen(output_file) - 1);
}

set terminal output_file_type size y_resolution, x_resolution font font.", ".font_size background rgb background_colour;
set output output_file.".".output_file_type;
set key off;
set xlabel x_label;
set ylabel y_label;
set tics scale 0;

if (is_folder) {
    pwd = GPVAL_PWD;
    cd data_folder;
    plot for [data in data_files] data using 0:1 with lines linewidth 2 linecolor rgbcolor 'blue';
    cd pwd;
} else {
    plot data_folder using 0:1 with lines linewidth 2 linecolor rgbcolor 'blue';
};
