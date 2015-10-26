% Plot ranges for iphnumvals
Evals = [1:0.05:5];
Dvals = [-5:0.2:-3,-2.8:0.05:-1.2, -1:0.03:1, 1.2:0.05:2.8, 3:0.1:2]; 

fprintf('Projected: %4.2fm. \r', length(Evals)*length(Dvals)/60)
PlotVals = iphnumvals(Evals, Dvals, 60, 0.5, 0.01);
surfc(Evals, Dvals, PlotVals)