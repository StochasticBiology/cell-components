reset
set multiplot
# svg 1024, 360
set border 3 lw 2
set xtics nomirror
set ytics nomirror

set style line 1 lc rgbcolor "#FF0000" pt 9 ps 0.5
set style line 2 lc rgbcolor "#FF00FF" pt 11 ps 0.5
set style line 3 lc rgbcolor "#0000FF" pt 5 ps 0.5
set style line 4 lc rgbcolor "#777777" pt 7 ps 0.5
set style line 5 lc rgbcolor "#000000" pt 13 ps 0.5

set tmargin screen 0.95
set bmargin screen 0.15
x1=0.1
x2=0.5
x3=0.57
x4=0.73
x5=0.76
x6=0.92
set size 0.5,1
set lmargin screen x1
set rmargin screen x2
set logscale y
set origin 0,0
set xlabel "Time" 
set ylabel "E(m)" offset 2.5,0
set key at graph 1.07,graph 0.4
set ytics ("10²" 100, "10³" 1000, "10⁴" 10000, "10⁵" 100000)
set label 1 at graph 0.05, graph 0.1 "A" font "Arial Black, 20" front
plot "./bneckbignuma.dat" w p ls 1 title "λ₁ = λ₁'", "./bneckbignumb.dat" w p ls 2 title "λ₁ = 0.02", "./bneckbignumc.dat" w p ls 3 title "λ₁ = 0.03", "./bneckbignumd.dat" w p ls 4 title "λ₁ = 0.04", "./bneckbignume.dat" w p ls 5 title "λ₁ = 0.05", "./bneckanala.dat" w l ls 1 notitle, "./bneckanalb.dat" w l ls 2 notitle, "./bneckanalc.dat" w l ls 3 notitle, "./bneckanald.dat" w l ls 4 notitle, "./bneckanale.dat" w l ls 5 notitle

set size 0.25,1
unset logscale y
set yrange [0:0.14]
set xtics 100
set xrange [0:250]
set ytics auto
set origin 0.5,0
set lmargin screen x3
set rmargin screen x4
m0 = 100000
n1 = 12
tau = 10
unset key
bneck(x) = 2**(-n1)*exp(n1*x*tau)*m0
pred(x) = sqrt(m0*(bneck(x)-m0)/bneck(x) * ((bneck(x)/m0)**(1./n1) + 1.)/((bneck(x)/m0)**(1./n1) - 1.))/m0
set ylabel "CV(m)"
set label 1 at graph 0.05, graph 0.95 "B" font "Arial Black, 20" front
plot pred(0.0693147) w l ls 1, pred(0.02) w l ls 2, pred(0.03) w l ls 3, pred(0.04) w l ls 4, pred(0.05) w l ls 5, "./bneckbignuma.dat" u 1:(sqrt($3)/$2) w p ls 1, "./bneckbignumb.dat" u 1:(sqrt($3)/$2) w p ls 2, "./bneckbignumc.dat" u 1:(sqrt($3)/$2) w p ls 3, "./bneckbignumd.dat" u 1:(sqrt($3)/$2) w p ls 4, "./bneckbignume.dat" u 1:(sqrt($3)/$2) w p ls 5, "./bneckanala.dat" u 1:(sqrt($3)/$2) w l ls 1, "./bneckanalb.dat" u 1:(sqrt($3)/$2) w l ls 2, "./bneckanalc.dat" u 1:(sqrt($3)/$2) w l ls 3, "./bneckanald.dat" u 1:(sqrt($3)/$2) w l ls 4, "./bneckanale.dat" u 1:(sqrt($3)/$2) w l ls 5

set origin 0.75,0
set lmargin screen x5 
set rmargin screen x6
set logscale x
set xlabel "b"
set border 9 lw 2
set y2tics auto nomirror
#set y2label "CV(m)"
unset ylabel
unset ytics
set xrange [1e2:5e5]
set y2range [0:0.14]
set xtics ("10²" 100, "10³" 1000, "10⁴" 10000, "10⁵" 100000)
set label 1 at graph 0.05, graph 0.1 "C" font "Arial Black, 20" front
plot  (pred( log(2**n1 * x / m0) / (n1 * tau))) axes x1y2 w l ls 1, "bneckcurve.dat" u (bneck($1)):(pred($1)) axes x1y2 w p ls 5 ps 1