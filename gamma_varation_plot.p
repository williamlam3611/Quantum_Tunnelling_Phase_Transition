cd path;
set print "-";

set terminal svg size 1680,1050 name output_file font "Helvetica,16" background rgb 'white';
set output output_file.".svg";
set key off;
set tics scale 0;

set xlabel x_label;

set ylabel y_label;


set style line 2  lc rgb '#0025ad' lt 1 lw 2.5 # --- blue
set style line 3  lc rgb '#0042ad' lt 1 lw 2.5 #      .
set style line 4  lc rgb '#0060ad' lt 1 lw 2.5 #      .
set style line 5  lc rgb '#007cad' lt 1 lw 2.5 #      .
set style line 6  lc rgb '#0099ad' lt 1 lw 2.5 #      .
set style line 7  lc rgb '#00ada4' lt 1 lw 2.5 #      .
set style line 8  lc rgb '#00ad88' lt 1 lw 2.5 #      .
set style line 9  lc rgb '#00ad6b' lt 1 lw 2.5 #      .
set style line 10 lc rgb '#00ad4e' lt 1 lw 2.5 #      .
set style line 11 lc rgb '#00ad31' lt 1 lw 2.5 #      .
set style line 12 lc rgb '#00ad14' lt 1 lw 2.5#      .
set style line 13 lc rgb '#09ad00' lt 1 lw 2.5 # --- green

# using (column(0) + x_offset)
plot for [n=1:5] data_file.".dat" using (column(6)):n w lines ls n;



#plot data_file.".dat" using (column(0) + x_offset):1 with lines linewidth 2 linecolor rgbcolor 'blue';
#plot data_file.".dat" using (column(0) + x_offset):2 with lines linewidth 2 linecolor rgbcolor 'green';
#plot data_file.".dat" using (column(2) + x_offset):3 with lines linewidth 2 linecolor rgbcolor 'red';
#plot data_file.".dat" using (column(3) + x_offset):4 with lines linewidth 2 linecolor rgbcolor 'yellow';
#plot data_file.".dat" using (column(4) + x_offset):5  with lines linewidth 2 linecolor rgbcolor 'purple';
