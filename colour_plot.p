cd path;
set print "-";

files = system("ls -1 *.dat")

set terminal svg size 1680,1050 name output_file font "Helvetica,16" background rgb 'white';
set output output_file.".svg";
set key off;
set tics scale 0;

set xlabel x_label;
set ylabel y_label;

rgb(r,g,b) = 65536 * int(r * 256) + 256 * int(g * 256) + int(b * 256)
plot for [data in files] data using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
