function [CriticalSlowingPlot, CriticalSlowingTimes] = Slowing(Evals, Dvals, varargin)
	% produces an array of slowing times for values of E and D
	try
		tlist = varargin{1};
	catch ME
		tlist = 0:1:600;
	end
	CriticalSlowingTimes = [];
	Ecounter = 1;
	Dcounter = 1;
	tic
	for E = Evals;
		Dcounter = 0;
		for D = Dvals;
			[a, CriticalSlowingTime, c] = CriticalSlowingComparison(E, D, 0, 0, tlist);
			CriticalSlowingTimes = horzcat(CriticalSlowingTimes, CriticalSlowingTime);
			if length(Dvals)>1
				t = toc;
				fprintf('%4.0fth D out of%4.0f after %4.2fs, E = %4.2f \n', Dcounter, length(Dvals), t, E)
				Dcounter = Dcounter+1;
			end

		end
		t = toc;
		fprintf('%4.0fth E out of%4.0f after %4.2fs \n', Ecounter, length(Evals), t)
		Ecounter = Ecounter+1;
	end
	t = toc;
	fprintf('\n %4.2fm for %d computations \n', t/60, length(Evals)*length(Dvals)*length(tlist))
	if length(Dvals == 1)
		CriticalSlowingPlot = plot(Evals, CriticalSlowingTimes);
	else
		CriticalSlowingPlot = surf(Evals, Dvals, CriticalSlowingTimes);
	end
