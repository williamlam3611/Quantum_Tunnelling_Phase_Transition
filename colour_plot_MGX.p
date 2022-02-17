set print "-";
getValue(row, col, filename) = system("awk "."'{if (NR == '".row."') print $'".col."'}'"." '".filename."'");

max_y   = getValue(8, 2, meta_file.".dat");
min_y   = getValue(8, 1, meta_file.".dat");
num_row = getValue(2, 3, meta_file.".dat");
num_col = getValue(2, 4, meta_file.".dat");

set terminal svg size 1680,1050 name output_file font "Helvetica,16" background rgb 'white';
set output path."/".output_file.".svg";
set key off;
set tics scale 0;

rgb(r,g,b) = 65536 * int(r * 256) + 256 * int(g * 256) + int(b * 256)

set multiplot layout 1, 3;
set size (0.2 * sqrt(2) * 1050 / 1680), 0.2
unset xlabel
unset ylabel
unset border
unset xtics
unset ytics
set origin 0.1, 0.2
set label 1 "dyz" at screen 0.1, 0.2 centre textcolor rgb "black";
set label 2 "dxz" at screen 0.277, 0.2 centre textcolor rgb "black";
set label 3 "dxy" at screen 0.188, 0.4 centre textcolor rgb "black";
plot 'colour_triangle.png' binary filetype=png with rgbalpha

set border
set label 3 "K_{M Γ X} [π/a]" at screen (2 - sqrt(2)), 0.02 centre;
set xtics (sprintf("%3.2f", sqrt(2) * getValue(4, 3, meta_file.".dat")) 0, \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file.".dat")) (num_col / 2) - ((sqrt(2) / 2) * (num_col / 2)), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file.".dat")) (num_col / 2) - ((sqrt(2) / 4) * (num_col / 2)), \
    "0.00" (num_col / 2), \
    sprintf("%3.2f", 0.50 * getValue(4, 3, meta_file.".dat")) ((3 * num_col) / 4), \
    sprintf("%3.2f", 1.00 * getValue(4, 3, meta_file.".dat")) (num_col - 2));

set ytics
set ylabel "E - E_{cbm} [ev]";

cd path;
files = system("ls -1 *.dat")

set size (2 - sqrt(2)), (1 - 0.03);
set origin 0, 0.03
set rmargin 0
set xrange [0:(num_row / 2)]
set yrange [min_y:max_y]
plot for [data in files] data using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;

set size (sqrt(2) - 1), (1 - 0.03);
set lmargin 0
unset rmargin
set origin (2 - sqrt(2)), 0.03;
unset ylabel
set ytics font ",0"
set xrange [((num_row / 2) + 1):num_row] 
plot for [data in files] data using 0:1:(rgb($2, $3, $4)) with lines lc rgb variable;
unset multiplot
