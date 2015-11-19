% Script for doing long runs of iphnumvals and saving to a file
% Plot ranges for iphnumvals
Evals = [4];
Dvals = [-1:1]; 

fprintf('Projected: %4.2fm. \r', length(Evals)*length(Dvals)/60)
PlotVals = iphnumvals(Evals, Dvals, 60, 0.5, 0.01);
save('PlotValues.mat', 'PlotVals')
