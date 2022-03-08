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

if (is_folder) {
    pwd = GPVAL_PWD;
    cd folder;
    data_files = system("ls -1 *.dat");
    cd pwd;
};

if (substr(output_file, strlen(output_file), strlen(output_file)) eq "/") {
    output_file = substr(output_file, 1, strlen(output_file) - 1);
} 
#set terminal output_file_type size y_resolution, x_resolution font font.", ".font_size background rgb background_colour;


#set terminal povray
#set terminal x11; 
set terminal gif animate delay 5 loop 0 optimize
#set output "rot.gif"


set output output_file.".".output_file_type;
set key off;
set mouse mouseformat 3; 

set xlabel x_label;
set ylabel y_label;
#set yrange [:0]


rgb(r,g,b) = 65536 * int(r * 256) + 256 * int(g * 256) + int(b * 256);

#set multiplot layout 1, 2;

set size 1, 1;
pwd = GPVAL_PWD;
cd folder;
#splot data_folder.".dat" using (column(5)):(column(6)):1:(rgb($2, $3, $4)) with lines lc rgb variable; 


n = 100
do for [i=1:n] {
   set view 60, i*360/n
   splot data_folder.".dat" using (column(5)):(column(6)):1:(rgb($2, $3, $4)) with lines lc rgb variable;

}




cd pwd; 

set output;

#unset xlabel;
#unset ylabel;
#unset border;
#unset xtics;
#unset ytics;
#set size triangle_size, triangle_size;
#set origin x_triangle, y_triangle;
#set label 1 column_label_3 at screen x_triangle, y_triangle centre textcolor rgb "black";
#set label 2 column_label_2 at screen x_triangle + triangle_size, y_triangle centre textcolor rgb "black";
#set label 3 column_label_1 at screen x_triangle + triangle_size / 2, y_triangle + triangle_size centre textcolor rgb "black";
#plot 'colour_triangle.png' binary filetype=png with rgbalpha;

save folder.'3d_plot.gnu';

pause mouse keypress;  
#unset multiplot;

save folder.'animation.gnu'

