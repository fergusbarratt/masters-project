function [iphnumvals] = iphnumvals(Eplotvals, Dplotvals, varargin)
		
	% Generates an array of expected cavity photon numbers for a given range of driving strengths.  

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
	
	% set default values
	try 
		N = varargin{1};
	catch ME
		N = 60;
	end
	try 
		g = varargin{2};
	catch ME
		g = 25;
	end
	try 
		kappa = varargin{3};
	catch ME
		kappa = 0.5;
	end

	tic

	%generate identities and constants
	ida = identity(N); 
	idatom = identity(2);

	% Define cavity field and atomic operators
	a = tensor(destroy(N),idatom);
	sm = tensor(ida,sigmam);

	% Generate Hamiltonian components outside for loop
	Hint =  g*(sm' * a + sm * a');
	HBare = (sm' * sm + a' * a);
	Collapse =  sqrt ( 2 * kappa ) * a;
	ConjCollapse = Collapse' * Collapse;

	function [ssiphnum] = ssiphnum(Hint, HBare, DriveStrength, Detuning, Collapse, ConjCollapse)
		% takes forms of bare and interaction hamiltonians, drive strength and collapse operator, solves for the density matrix in steady state
		% then generates steady state photon number

			H = - Detuning * HBare + Hint + DriveStrength * (a' + a);
			L = -1i * ( spre(H)- spost(H) ) + spre(Collapse)*spost(Collapse')-0.5*spre(ConjCollapse)-0.5*spost(ConjCollapse);
			rhoss = steady(L);
			ssiphnum = real(expect(ConjCollapse, rhoss))/2;

	end


	% vars for counting progress
	n = 1; 
	ds = length(Dplotvals); 
	es = length(Eplotvals);
	numberofvalues = es*ds;
	projectedtime = numberofvalues;

	% for loop initialisation - empty arrays for building
	row = [];
	iphnumvals = [];

	for E = Eplotvals
		row = [];
		t = toc;
		fprintf(' %d/%d: %4.2fs elapsed out of projected %4.1fs \r', n, es, t, projectedtime);
		for D = Dplotvals
			row = vertcat(ssiphnum(Hint, HBare, E, D, Collapse, ConjCollapse), row);
			% dispstat(D);
			fprintf('|');
		end
		iphnumvals = horzcat(iphnumvals, row);
		n = n+1;
		if n==2;
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