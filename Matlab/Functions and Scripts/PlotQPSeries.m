function [plo] = PlotQPSeries(Erange, Detrange, functype, varargin)
	% function to plot a series of Q functions along an E-line in figure from Carmichael paper
	% parameter values. inclusive of both ends of E range. Deletes old figure, specifies the figure size and renormalizes q at each stef
	% TODO: autowidth
	% TODO: matrices
	% to add to paper PlotQPSeries(4.5:0.2:5.5, 0, 'W'), N=60, g=10, kappa=1, gamma=1 - just greater than 1
	N = 60;
	g = 10;
	kappa = 1;
	gamma = 1.5;

	step_size = 0.1;
	cols = min(3, max(size(Erange))); % number of columns in E subplots
	rows = min(3, max(size(Detrange))); % number of rows in Det subplots
	
	try 
		plottype = varargin{3};
	catch ME
		plottype = 'c';
	end
	try
		tracesys = varargin{4};
	catch ME
		tracesys = 1;
	end
	try 
		Xrange = varargin{1};
		Yrange = varargin{2};
	catch ME
		Xrange = -8:step_size:8;
		Yrange = -8:step_size:8;
	end

	% DELETES CURRENT FIGURE BEFORE DRAWING NEW ONE.
	delete(gcf);
	figure('position', [0, 0, 720, 800]);
	for E = Erange
		for det = Detrange
			if functype == 'Q'
				h = real(qfunc(ptrace(rhoss(E/kappa, det/kappa, N, g, kappa, gamma), tracesys), Xrange, Yrange));
				h = h/volintegral(h, Xrange, Yrange); % Renormalise Q function at each step
				volintegral(h, Xrange, Yrange); %remove semicolon to output volume under Q
				chk = volintegral(h, Xrange, Yrange); 
				 	if (chk>1.05 || chk<0.95) %print if normalisation is further from 1 than some value
				 		'normalisation error'
				 		chk
				 	end
			elseif functype == 'W'
				h = real(wfunc(ptrace(rhoss(E/kappa, det/kappa, N, g, kappa, gamma), tracesys), Xrange, Yrange));
				h = h/volintegral(h, Xrange, Yrange); % Renormalise Q function at each step
				chk = volintegral(h, Xrange, Yrange);
			 	if (chk>1.05 || chk<0.95)
			 		'normalisation error'
			 		chk
			 	end
			end
			if isvector(Erange) && isscalar(Detrange)
				[P, Q] = find(Erange==E);
				subplot(ceil(max(size(Erange))/cols), cols, Q)
			elseif isvector(Detrange) && isscalar(Erange)
				[P, Q] = find(Detrange==det);
				subplot(rows, ceil(max(size(Detrange))/rows), Q)
			else
				error('works only for lines (for now)')
			end
				switch plottype
					case 'cf'
						contourf(Xrange, Yrange, h);
					case 'c'
						fprintf('|');
						plo = surf(Xrange, Yrange, h);
						plo.LineStyle = 'none';
						view(2);
						colormap('jet');
						colorbar
						ylim=get(gca,'YLim');
						xlim=get(gca,'XLim');
						text((xlim(1)-1),(ylim(1)-2.5),[num2str(E) ',' num2str(det) '; g: ' num2str(g) ', kappa: ' num2str(kappa) ', gamma: ' num2str(gamma)], 'VerticalAlignment','top', 'HorizontalAlignment','left')
						xlabel('Re(Q)')
						ylabel('Im(Q)')
					case 's'
						surfc(Xrange, Yrange, h);
					case 'm'
						meshc(Xrange, Yrange, h);
				end
			end
		end
	fprintf('\n')
