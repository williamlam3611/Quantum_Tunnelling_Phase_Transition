#cd path;
#pwd = GPVAL_PWD
#cd pwd
set print "-";

set terminal svg size 1680,1050 font "Helvetica,16" background rgb 'white';
set output output_file.".svg";
set key off;
set tics scale 0;

set xlabel x_label;

set ylabel y_label;

set multiplot layout 1,1
plot data_file_1 every ::0::15 using (column(0) + x_offset):data_column with lines linewidth 2 linecolor rgbcolor 'blue';

replot data_file_2 every ::0::15 using (column(0) + x_offset):data_column with lines linewidth 2 linecolor rgbcolor 'red';

unset multiplot;
