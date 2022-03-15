set print "-";

is_folder = 0;
if (!exists("data_folder"))               data_folder        = GPVAL_PWD;
if (strstrt(substr(data_folder, 2, strlen(data_folder)), ".") == 0) is_folder = 1;
if (!exists("output_file") && is_folder)  output_file        = data_folder;
if (!exists("output_file") && !is_folder) output_file        = data_folder[1:strstrt(data_folder, ".") - 1];
if (!exists("output_file_type"))          output_file_type   = "svg";
if (!exists("y_resolution"))              y_resolution       = 1680;
if (!exists("x_resolution"))              x_resolution       = 1050;
if (!exists("font"))                      font               = "Helvetica";
if (!exists("font_size"))                 font_size          = "16";
if (!exists("background_colour"))         background_colour  = "white";
if (!exists("x_label"))                   x_label            = "X";
if (!exists("y_label"))                   y_label            = "Y";
if (!exists("z_label"))                   z_label            = "Z";
if (!exists("has_axis"))                  has_axis           = 1;
if (!exists("colour_box_x"))              colour_box_x       = 0.1;
if (!exists("colour_box_y"))              colour_box_y       = 0.2;
if (!exists("colour_box_size"))           colour_box_size    = 0.2;
if (!exists("colour_box_label_1"))        colour_box_label_1 = "red";
if (!exists("colour_box_label_2"))        colour_box_label_2 = "blue";
if (!exists("colour_box_label_3"))        colour_box_label_3 = "green";
if (!exists("colour_box_type"))           colour_box_type    = "triangle";
if (colour_box_type eq "box")             colour_box_type    = "box_rgb";
if (!exists("plot_type"))                 plot_type          = "surface";
if (!exists("transparency"))              transparency       = 1.0;
if (!exists("pitch"))                     pitch              = 45;
if (!exists("yaw"))                       yaw                = 45;


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

set xlabel x_label rotate parallel;
set ylabel y_label rotate parallel;
set zlabel z_label rotate parallel;

if (exists("x_min")) set xrange[x_min:];
if (exists("y_min")) set yrange[y_min:];
if (exists("z_min")) set zrange[z_min:];
if (exists("x_max")) set xrange[:x_max];
if (exists("y_max")) set yrange[:y_max];
if (exists("z_max")) set zrange[:z_max];

rgb(r,g,b) = 65536 * int(r * 255) + 256 * int(g * 255) + int(b * 255);


if (strstrt(colour_box_type, "none") == 0) {
    ratio = real(word(system("file ".colour_box_type.".png -b"), 4)) / real(word(system("file ".colour_box_type.".png -b"), 6));
    set pixmap 1 colour_box_type.".png" at screen colour_box_x, colour_box_y \
        size screen (colour_box_size * ratio), colour_box_size front center;
    if (strstrt(colour_box_type, "box") > 0) {
        set label 1 colour_box_label_1 \
            at screen colour_box_x, (colour_box_y - 0.01 - colour_box_size / 2.0) \
            centre textcolor rgb "black";
        set label 2 colour_box_label_2 \
            at screen colour_box_x, (colour_box_y + 0.01 + colour_box_size / 2.0) \
            centre textcolor rgb "black";
    };
    if (strstrt(colour_box_type, "triangle") > 0) {
        set label 1 colour_box_label_3 \
            at screen (colour_box_x - 0.01 -ratio * colour_box_size / 2.0), (colour_box_y - 0.01 - colour_box_size / 2.0) \
            centre textcolor rgb "black";
        set label 2 colour_box_label_2 \
            at screen (colour_box_x + 0.01 + ratio * colour_box_size / 2.0), (colour_box_y - 0.01 -colour_box_size / 2.0) \
            centre textcolor rgb "black";
        set label 3 colour_box_label_1 \
            at screen colour_box_x, (colour_box_y + 0.01 + colour_box_size / 2.0) \
            centre textcolor rgb "black";
    };
};

unset colorbox;
if (has_axis + 0 == 0) {
    unset border;
    unset xtics;
    unset ytics;
    unset xlabel;
    unset ylabel;
};
set size 1, 1;
set view (pitch + 0.0), (yaw + 0.0);
set xyplane 0.0;
if (strstrt(plot_type, "surface") > 0) {
    set dgrid3d sqrt(max_row), sqrt(max_row);
    set style fill transparent solid (transparency + 0.0);
    if (is_folder) {
        pwd = GPVAL_PWD;
        cd data_folder;
        splot for [i = 1:words(data_files)] word(data_files, i) \
                  using 1:2:3:(rgb(column(4), column(5), column(6))) \
                  with pm3d linecolor rgb variable;
        cd pwd;
    } else {
        splot data_folder using 1:2:3:(rgb(column(4), column(5), column(6))) with pm3d linecolor rgb variable;
    };
};
if (strstrt(plot_type, "fence") > 0) {
    if (is_folder) {
        pwd = GPVAL_PWD;
        cd data_folder;
        splot for [i = 1:words(data_files)] for [j = 1:sqrt(max_row)] word(data_files, i) every sqrt(max_row)::j \
            using 1:2:3:(rgb(column(4), column(5), column(6))) \
            with lines linecolor rgb variable;
        cd pwd;
    } else {
        splot data_folder for [j = 1:sqrt(max_row)] every sqrt(max_row)::j \
            using 1:2:3:(rgb(column(4), column(5), column(6))) \
            with lines linecolor rgb variable;
    };
};
if (strstrt(plot_type, "seperate") > 0) {
    system("mkdir -p '".output_file."'")
    if (has_axis + 0 == 1) {set xlabel y_label;};
    if (has_axis + 0 == 1) {set ylabel z_label;};
    if (exists("y_min")) {set xrange[y_min:];};
    if (exists("z_min")) {set yrange[z_min:];};
    if (exists("y_max")) {set xrange[:y_max];};
    if (exists("z_max")) {set yrange[:z_max];};
    do for [j=1:sqrt(max_row)] {
        unset output;
        set output output_file."/".word(x_label, 1)."_".j.".".output_file_type;
        if (is_folder) {
            pwd = GPVAL_PWD;
            cd data_folder;
            plot for [i = 1:words(data_files)] word(data_files, i) every sqrt(max_row)::j \
                using 2:3:(rgb(column(4), column(5), column(6))) \
                with lines linecolor rgb variable;
            cd pwd;
        } else {
            plot data_folder every sqrt(max_row)::i \
                using 2:3:(rgb(column(4), column(5), column(6))) \
                with lines linecolor rgb variable;
        };
    }; 
};
