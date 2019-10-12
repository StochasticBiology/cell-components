# simply grabs the last values from the simulated datasets and outputs them to a file with a reference depletion value
awk '{if($1 == 239) print 0.0693147, $1, $2, $3}' bneckbignuma.dat > bneckcurve.dat
awk '{if($1 == 239) print 0.02, $1, $2, $3}' bneckbignumb.dat >> bneckcurve.dat
awk '{if($1 == 239) print 0.03, $1, $2, $3}' bneckbignumb.dat >> bneckcurve.dat
awk '{if($1 == 239) print 0.04, $1, $2, $3}' bneckbignumb.dat >> bneckcurve.dat
awk '{if($1 == 239) print 0.05, $1, $2, $3}' bneckbignumb.dat >> bneckcurve.dat
