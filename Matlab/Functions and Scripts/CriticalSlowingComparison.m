function [h, slowingtime1, slowingtime2] = CriticalSlowingComparison(E1, det1, E2, det2, tlist, varargin)
	
	TDSols1 = timedependentsoln(E1, det1, tlist);
	TDSols2 = timedependentsoln(E2, det2, tlist);

	[r1, ss1] = timedependentsoln(E1, det1, tlist);
	[r2, ss2] = timedependentsoln(E2, det2, tlist);
	
	sss1 = ones(size(r1))*ss1;
	sss2 = ones(size(r1))*ss2;

	try
		ifplot = varargin{1}
	catch 
		ifplot = 'n';
	end
	try
		concrit = varargin{2};
	catch 
		concrit = 0.0001;
	end
	try
		res = varargin{3};
	catch
		res = 100;
	end


	if ifplot == 'y'
		h = plot(tlist, TDSols1, tlist, TDSols2, tlist, sss1,  tlist, sss2);
	else
		h = 'no plot';
	end

	ind1 = [];
	ind2 = [];
	for q = r1
		if abs((q-ss1))<=abs(concrit)
			ind1 = horzcat(ind1, find(TDSols1==q));
			if length(ind1)>=res
				break
			end
		else
			ind1 = [];
		end
	end
	for q = r2
		if abs((q-ss2))<=abs(concrit)
			ind2 = horzcat(ind2, find(TDSols2==q));
			if length(ind2)>=res
				break
			end
		else
			ind2 = [];
		end
	end

	% fprintf('line 1 has %4.0f points out of %4.0f \r', length(ind1), length(tlist))
	% fprintf('line 2 has %4.0f points out of %4.0f \r', length(ind2), length(tlist))
	if length(ind1)>=res
		slowingtime1 = tlist(min(ind1));
	else
		slowingtime1 = 0;
	end
	if length(ind2)>=res
		slowingtime2 = tlist(min(ind2));
	else
		slowingtime2 = 0;
	end