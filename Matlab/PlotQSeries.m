function [plo] = plotQseries(Erange, det, plottype, varargin)
	% function to plot a series of Q functions along an E-line in figure from Carmichael paper
	% parameter values.
	N = 60;
	g = 10;
	kappa = 10;

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
		[P, Q] = find(Erange==E);
		subplot(ceil(max(size(Erange))/cols), cols, Q), 
			switch plottype
				case 'cf'
					contourf(Xrange, Yrange, qfunc(rhoss(E, det, N, g, kappa), Xrange, Yrange));
				case 'c'
					contour(Xrange, Yrange, qfunc(rhoss(E, det, N, g, kappa), Xrange, Yrange));
				case 's'
					surfc(Xrange, Yrange, qfunc(rhoss(E, det, N, g, kappa), Xrange, Yrange));
				case 'm'
					meshc(Xrange, Yrange, qfunc(rhoss(E, det, N, g, kappa), Xrange, Yrange));

			end
	end