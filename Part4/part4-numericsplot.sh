reset
# 1024,360

set border 3 lw 2
set xtics auto nomirror
set ytics auto nomirror
#set grid

set style line 1 lc rgbcolor "#FF0000" pt 9 ps 0.5
set style line 2 lc rgbcolor "#0000FF" pt 11 ps 0.5
set style line 3 lc rgbcolor "#FF00FF" pt 5 ps 0.5
set style line 4 lc rgbcolor "#000000" pt 7 ps 0.5

set multiplot
set xrange [0:39]
set size 0.5,1
set origin 0,0
set logscale y
set yrange [90:300]
set ytics (100,200,300)
set xlabel "t"
set ylabel "E(m)"
#set yrange [0:400]
set label 1 at graph 0.05, graph 0.95 "A" font "Arial Black, 20" front
plot "partitioning-regime-0.dat" u 1:2 ls 1 title "Binomial", "partitioning-regime-3.dat" u 1:2 ls 2 title "Subtractive", "partitioning-regime-1.dat" u 1:2 ls 3 title "Cluster", "partitioning-regime-2.dat" u 1:2 ls 4 title "Deterministic",  "analytic-regime-0.dat" u 1:2 w l ls 1 notitle, "analytic-regime-3.dat" u 1:2 w l ls 2 notitle, "analytic-regime-1.dat" u 1:2 w l ls 3 notitle, "analytic-regime-2.dat" u 1:2 w l ls 4 notitle

 #, "test.bak" u 1:2 title "Test"
set key bottom right
set origin 0.5,0
set logscale y
set ytics auto
set yrange [*:*]
set ylabel "V(m)"
set label 1 at graph 0.05, graph 0.95 "B" font "Arial Black, 20" front
plot "partitioning-regime-0.dat" u 1:3 ls 1 title "Binomial", "partitioning-regime-3.dat" u 1:3 ls 2 title "Subtractive", "partitioning-regime-1.dat" u 1:3 ls 3 title "Cluster", "partitioning-regime-2.dat" u 1:3 ls 4 title "Deterministic",  "analytic-regime-0.dat" u 1:3 w l ls 1 notitle, "analytic-regime-3.dat" u 1:3 w l ls 2 notitle, "analytic-regime-1.dat" u 1:3 w l ls 3 notitle, "analytic-regime-2.dat" u 1:3 w l ls 4 notitle

#, "test.bak" u 1:3 title "Test"
