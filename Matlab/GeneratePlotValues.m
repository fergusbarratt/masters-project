% Plot ranges for iphnumvals
Evals = [1:1:8, 8.2:0.4:10, 10.2:0.2:12, 12.1:0.1:16, 17:25];
Dvals = [-20:1:-7,-6.8:0.1:-1.2, -1:0.05:1, 1.2:0.1:6.8, 7:1:20]; 

fprintf('Projected: %4.2fm. \r', length(Evals)*length(Dvals)/60)
PlotVals = iphnumvals(Evals, Dvals);
surfc(Evals, Dvals, PlotVals);