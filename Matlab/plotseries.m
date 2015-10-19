function [plot] = plotseries(Erange, det, plottype, varargin)
	% function to plot a series of Q functions along an E-line in figure from Carmichael paper
	kappa = 25;
	g = 0.5;
	N = 60; % default vals of cavity parameters
	cols = 5; % number of columns in subplots
	xrangeexists = exist varargin{1};

	if (xrangeexists)
		Xrange = varargin{1};
	else
		Xrange = -5:0.5:5;	% default Q plotrange
	end
	if (yrangeexists)
		Yrange = varargin{2}
	else	
		Yrange = -5:0.5:5;
	end
	figure
	for E = Erange
		[P, Q] = find(Erange==E);
		subplot(ceil(max(size(Erange))/cols), cols, Q), 
			switch plottype
				case 'cf'
					contourf(Xrange, Yrange, Qss(E, kappa, g, det, N, Xrange, Yrange));
				case 'c'
					contour(Xrange, Yrange, Qss(E, kappa, g, det, N, Xrange, Yrange));
				case 's'
					surfc(Xrange, Yrange, Qss(E, kappa, g, det, N, Xrange, Yrange));
				case 'm'
					meshc(Xrange, Yrange, Qss(E, kappa, g, det, N, Xrange, Yrange));

			end
	end
