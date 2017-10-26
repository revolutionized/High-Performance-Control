# Open a new window to plot the figure (similar to MATLABs figure(0))
set term windows 0
# Here we plot the results of the ODE
set xlabel "Time (seconds)"
set ylabel "Velocity"
plot 'ExactEulerResult.dat' using 1:2 with lines
# Replot is like MATLABs hold on
replot 'MarkovEulerResult.dat' using 1:2 with lines
# Open a new window to plot the figure (similar to MATLABs figure(1))
set term windows 1
# Here we plot the optimal control values for the ODE
set xlabel "Time (seconds)"
set ylabel "Control value"
plot 'ExactControl.dat' using 1:2 with lines
replot 'MarkovControl.dat' using 1:2 with lines