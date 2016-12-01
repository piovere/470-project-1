set xlabel "Position (dx)"
set ylabel "Flux (abritrary units)"
set grid
plot 'Flux1Group.dat' using 1 title '1st Group', \
     'Flux2Group.dat' using 1 title '2nd Group', \
     'Flux3Group.dat' using 1 title '3rd Group', \
     'Flux4Group.dat' using 1 title '4th Group', \

     
