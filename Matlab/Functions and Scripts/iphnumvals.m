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
		expectertype = '22';
	end
	try
		N = varargin{2};
	catch
		N = 40;
	end
	try
		g = varargin{3};
	catch
		g = 10;
	end
	try
		kappa = varargin{4};
	catch
		kappa = 1;
	end
	fprintf(['starting run, mode ' expectertype ', g = %d, kappa = %d, N = %d \n'], g, kappa, N)
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
		atom = ptrace(rho, 2);
		switch expectertype
		 	case 'iphnum'
		 		expecter = real(expect(C1dC1, rho)/(kappa)); % Intracavity Photon Number
		 	case 'purity'
				expecter = real(trace(rho^2)); % Purity
			case '11'
				expecter = abs(atom(1, 1));
			case '22'
				expecter = abs(atom(2, 2));
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
	t =0;

	% for loop initialisation - empty arrays for building
	iphnumvals = zeros(ds, es);
	for rowind = 1:ds %for each e/row
		 row = zeros(1, es); %create a row of the right length

		for colind = 1:es % for each d/col index
			row(colind) = expecter(rhos(Hi, Hb, Eplotvals(colind)/kappa, Dplotvals(rowind)/kappa, C1, C1dC1));
			ptrace(rhos(Hi, Hb, Eplotvals(colind)/kappa, Dplotvals(rowind)/kappa, C1, C1dC1), 2) % set the values at (rowind, colind) to the value of expecter
			% dispstat(D);
			fprintf('|');
		end
		iphnumvals(rowind, :) = row; % add the e-row to the matrix
		n = n+1; % Counting and timing
		roundtime = toc;
		tic
		t = t + roundtime;
		if n>=1;
			projectedtime = (projectedtime+numberofvalues*roundtime)/2;
		end
		if (n>0)
			fprintf(' %d | Es: %d Ds: %d | %4.2fs elapsed out of projected %4.1fs \r', n, es, ds, t, projectedtime);
		end
	end

	fprintf('%4.2fs, %d values\n', t, numberofvalues)
end
