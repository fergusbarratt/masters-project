function [plo] = PlotQPSeries(Erange, det, functype, plottype, tracesys, varargin)
	% function to plot a series of Q functions along an E-line in figure from Carmichael paper
	% parameter values.
	N = 2;
	g = 25;
	kappa = 0.5;

	% default vals of cavity parameters
	cols = min(5, max(size(Erange))); % number of columns in subplots

	try
		Xrange = varargin{1};
		Yrange = varargin{2};
	catch ME
		Xrange = -5:0.3:5;
		Yrange = -5:0.3:5;
	end
	plo = figure;
	for E = Erange
		if functype == 'Q'
			h = real(qfunc(ptrace(rhoss(E, det, N, g, kappa), tracesys), Xrange, Yrange));
		elseif functype == 'W'
			h = real(wfunc(ptrace(rhoss(E, det, N, g, kappa), tracesys), Xrange, Yrange));
		end
		[P, Q] = find(Erange==E);
		subplot(ceil(max(size(Erange))/cols), cols, Q),
			switch plottype
				case 'cf'
					contourf(Xrange, Yrange, h);
				case 'c'
					contour(Xrange, Yrange, h);
				case 's'
					surfc(Xrange, Yrange, h);
				case 'm'
					meshc(Xrange, Yrange, h);

			end
		end
