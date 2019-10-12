#suits 1100,300
reset
set multiplot
set size 0.2,1
set border 3 lw 2

set xrange [0:30]
set yrange [400:1100]
set xtics nomirror
set ytics nomirror

set ylabel offset 2.5,0

set tmargin screen 0.95
set bmargin screen 0.15
x0 = 0.065
dx = 0.98/5
ddx = dx-0.07

unset key
set key at graph 1.2,1
set xlabel "t"
set ylabel "E(m)"
set lmargin screen x0+0*dx
set rmargin screen x0+0*dx+ddx
set label 1 at graph 0.05, graph 0.95 "A" font "Arial Black, 18" front
plot "tracenormal.txt" u 1:2 w l title "Theory", "c1normal.txt" u 1:2 lc rgbcolor "#AAAAAA" notitle w l, "trythisnormal.txt" u 1:2 title "Simulation" pt 7 ps 0.5 lc rgbcolor "#000000"

unset key
set xlabel "t"
set ylabel "V(m)"
set lmargin screen x0+1*dx
set rmargin screen x0+1*dx+ddx
set label 1 at graph 0.05, graph 0.95 "B" font "Arial Black, 18" front
plot "tracenormal.txt" u 1:3 w l,  "trythisnormal.txt" u 1:($3**2) pt 7 ps 0.5 lc rgbcolor "#000000"

#########

set xrange [0:30]
set yrange [81:107]
set xtics nomirror
set ytics nomirror

set xlabel "t"
set ylabel "E(m)"
set lmargin screen x0+2*dx
set rmargin screen x0+2*dx+ddx
set label 1 at graph 0.05, graph 0.95 "C" font "Arial Black, 18" front
plot "traceyeast.txt" u 1:2 w l title "Theory", "c1yeast.txt" u 1:2 lc rgbcolor "#AAAAAA" notitle w l, "trythisyeast.txt" u 1:2 every 2 title "Simulation" pt 7 ps 0.5 lc rgbcolor "#000000"

set xlabel "t"
set ylabel "V(m)"
set lmargin screen x0+3*dx
set rmargin screen x0+3*dx+ddx
set label 1 at graph 0.05, graph 0.95 "D" font "Arial Black, 18" front
plot "traceyeast.txt" u 1:3 w l title "Theory", "c1yeast.txt" u 1:3 lc rgbcolor "#AAAAAA" notitle w l, "c1yeast.txt" u 1:4 lc rgbcolor "#AAAAAA" notitle w l, "c1yeast.txt" u 1:5 lc rgbcolor "#AAAAAA" notitle w l,  "trythisyeast.txt" u 1:($3**2) every 2 pt 7 ps 0.5 lc rgbcolor "#000000" title "Simulation", "c2yeast.txt" u 1:2 lc rgbcolor "#0000FF" pt 2 ps 1 lw 3 notitle

unset key

set lmargin screen x0+4*dx
set rmargin screen x0+4*dx+ddx
set label 1 at graph 0.05, graph 0.95 "E" font "Arial Black, 18" front
set xlabel "m"
set ylabel "P(m)"
set yrange [*:0.016]
set xrange [500:1100]
set xtics 200
set ytics 0.005
set label 2 at 725, 0.015 "t = 1"
set label 3 at 850, 0.0142 "t = 3"
set label 4 at 1000, 0.012 "t = 5"
plot "pdfnormal.txt" u ($1 == 1 ? $2 : 1/0):3 w l lc rgbcolor "#FF0000", "" u ($1 == 3 ? $2 : 1/0):3 w l lc rgbcolor "#FF00FF",  "" u ($1 == 5 ? $2 : 1/0):3 w l lc rgbcolor "#0000FF", "pdfnormal2.txt" u ($1 == 1 && $3 > 1e-5 ? $2 : 1/0):3 every 3 w p pt 7 ps 0.5 lc rgbcolor "#000000", "" u ($1 == 3  && $3 > 1e-5 ? $2 : 1/0):3 every 3 w p pt 7 ps 0.5 lc rgbcolor "#000000",  "" u ($1 == 5  && $3 > 1e-5 ? $2 : 1/0):3 every 3 w p pt 7 ps 0.5 lc rgbcolor "#000000"
