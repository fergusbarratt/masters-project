function [iphnumvals] = iphnumvals(Eplotvals, Dplotvals, varargin)

	% Generates an array of some function applied at a given range of driving strengths. Optional: cavity Hilbert space size, coupling g, coupling to vacuum kappa

	% solves the problem of a coherent, cavity-qubit detuning zero driving of cavity with a two-level atom
	% Eplotvals = amplitudes of driving field to plot,
	% kappa = mirror coupling,
	% gamma = spontaneous emission rate,
	% g = atom-field coupling,
	% d = detuning between drive and cavity-qubit resonance
	% N = size of Hilbert space for intracavity field (zero to N-1 photons)
	% count1 = photocount rate of light leaking out of cavity count2 = spontaneous emission rate
	% iphnum = intracavity field
	% dispstat('', 'init')
	% HBAR=1

	% set default values
	try 
		expectertype = varargin{1};
	catch
		expectertype = 'iphnum';
	end
	try
		N = varargin{2};
	catch
		N = 30;
	end
	try
		g = varargin{3};
	catch
		g = 50;
	end
	try
		kappa = varargin{4};
	catch
		kappa = 1;
	end
	fprintf(['starting run, mode ' expectertype ', g= %d, kappa = %d, N = %d \n'], g, kappa, N)
	tic
	function [rhos] = rhos(Hint, HBare, DriveStrength, Detuning, Collapse, ConjCollapse)
		% takes forms of bare and interaction hamiltonians, drive strength and collapse operator, solves for the density matrix in steady state
		% then returns steady state density matrix

			H = - Detuning * HBare + Hint + DriveStrength * (a' + a);
			L = -1i * ( spre(H)- spost(H) ) + 2*spre(Collapse)*spost(Collapse')-spre(ConjCollapse)-spost(ConjCollapse);
			rhos = steady(L);
	end

	function [expecter] = expecter(rho)
		% function to evaluate at each rho
		% TODO: multiple functions evaluated
		switch expectertype
		 	case 'iphnum'
		 		expecter = real(expect(C1dC1, rho)/(kappa)); % Intracavity Photon Number
		 	case 'purity'
				expecter = real(trace(rho^2)); % Purity
		 	otherwise
		 		'no such expecter function'
		 end 
	end

	%generate identities and constants
	ida = identity(N);
	idatom = identity(2);

	% Define cavity field and atomic operators
	a = tensor(destroy(N),idatom);
	sm = tensor(ida,sigmam);

	% Generate Hamiltonian components outside for loop
	Hi =  g*(sm' * a + sm * a');
	Hb = (sm' * sm + a' * a);
	C1 =  sqrt ( kappa ) * a;
	C1dC1 = C1' * C1;

	% vars for counting progress
	n = 0;
	ds = length(Dplotvals);
	es = length(Eplotvals);
	numberofvalues = es*ds;
	projectedtime = numberofvalues;

	% for loop initialisation - empty arrays for building
	iphnumvals = [];

	for E = Eplotvals
		row = [];

		for D = Dplotvals
			row = vertcat(expecter(rhos(Hi, Hb, E/kappa, D/kappa, C1, C1dC1)), row);
			% dispstat(D);
			fprintf('|');
		end
		iphnumvals = horzcat(iphnumvals, row);
		n = n+1;
		t = toc;
		if (n>0)
		fprintf(' %d/%d: %4.2fs elapsed out of projected %4.1fs \r', n, es, t, projectedtime);
		end
		if n==1;
			projectedtime = t*es;
		end

	end

	% arrayfun -> should be slower, might be faster? weird bugs
	% generate matrices of vals for arrayfun
	% kappas = kappa*ones(length(Dplotvals), length(Eplotvals));
	% gs = g*ones(length(Dplotvals), length(Eplotvals));
	% %ds = d*ones(length(Eplotvals));
	% Ns = N*ones(length(Dplotvals), length(Eplotvals));
	% Eplotvalsd = repmat(Eplotvals, [length(Dplotvals), 1]);
	% Dplotvalsd = repmat(Dplotvals', [1, length(Eplotvals)]);

	% iphnumvals = arrayfun(@ssiphnum, kappas, gs, Dplotvalsd, Ns, Eplotvalsd);

	timeelapsed = toc;

	fprintf('%4.2fs, %d values\n', timeelapsed, numberofvalues)
end
