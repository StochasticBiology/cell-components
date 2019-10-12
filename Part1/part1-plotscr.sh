reset
# suits 800,400
set multiplot
set size 1,0.5
set border 3 lw 2

set xlabel "t"
set ylabel "E(m)"

set xrange [0:200]
set yrange [0:*]

set origin 0,0.5
set key top right
set ytics 50
set xtics nomirror
set ytics nomirror
set label 1 at graph 0.05, graph 0.95 "A" font "Arial Black, 20" front
plot "testnuma.dat" u 1:2 every 2 w p pt 7 ps 0.5  lc rgbcolor "#000000" title "Mean (simulation)", "" u 1:($2+sqrt($3)) every 2 w p pt 7 ps 0.5 lc rgbcolor "#FF8888" title "SD (simulation)", "" u 1:($2-sqrt($3)) every 2  w p pt 7 ps 0.5  lc rgbcolor "#FF8888" notitle, "testtheorya.dat" u 1:2 w l ls 1  lc rgbcolor "#000000" title "Mean (theory)", "" u 1:($2+sqrt($3)) w l ls 1  lc rgbcolor "#FF8888" title "SD (theory)", "" u 1:($2-sqrt($3)) w l ls 1  lc rgbcolor "#FF8888" notitle
set origin 0,0
unset key #set key bottom right
set ytics 100
set label 1 at graph 0.05, graph 0.95 "B" font "Arial Black, 20" front
plot "testnumc.dat" u 1:2 every 2 w p pt 7 ps 0.5  lc rgbcolor "#000000" title "Mean (simulation)", "" u 1:($2+sqrt($3)) every 2 w p pt 7 ps 0.5 lc rgbcolor "#FF8888" title "SD (simulation)", "" u 1:($2-sqrt($3)) every 2 w p pt 7 ps 0.5  lc rgbcolor "#FF8888" notitle, "testtheoryc.dat" u 1:2 w l ls 1  lc rgbcolor "#000000" title "Mean (theory)", "" u 1:($2+sqrt($3)) w l ls 1  lc rgbcolor "#FF8888" title "SD (theory)", "" u 1:($2-sqrt($3)) w l ls 1  lc rgbcolor "#FF8888" notitle