function [iphnumvals] = iphnumvals(Eplotvals, kappa, g, Dplotvals, N)
		
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

	tic

	%generate identities and constants
	ida = identity(N); 
	idatom = identity(2);

	% Define cavity field and atomic operators
	a = tensor(destroy(N),idatom);
	sm = tensor(ida,sigmam);

	function [ssiphnum] = ssiphnum(kappa, g, d, N, E)
		%Calculates the Liouvillian, finds the density matrix in steady state, outputs the photon number

		% get the hamiltonian for the values given
		% Hamiltonian cavity and qubit in resonance for input vals, d- cavity-qubit freq minus drive freq
		H = -d*(sm' *sm+ a'*a) + g*(sm'*a + a'*sm) + E*(a'+a); 

		% get the collapse operator for the values given, spont emission now commented out
		C = sqrt(2*kappa)*a;
		%C2 = sqrt(gamma)*sm;
		CdC = C'*C;
		%C2dC2 = C2'*C2;

		% Calculate the Liouvillian, spont emission lindbladians not included
		LH = -1i * (spre(H) - spost(H));
		L1 = spre(C)*spost(C')-0.5*spre(CdC)-0.5*spost(CdC); 
		% L2 = spre(C2)*spost(C2'')-0.5*spre(C2dC2)-0.5*spost(C2dC2); 
		L = LH+L1; %+L2;
		
		% Find steady state density matrix
		rhoss = steady(L);

		% expectation value of collapse product over 2kappa is expectation of number in cavity (QOT)
		ssiphnum = real(expect(CdC,rhoss)/(2*kappa));
	end

	% generate matrices of vals for arrayfun
	kappas = kappa*ones(length(Dplotvals), length(Eplotvals));
	gs = g*ones(length(Dplotvals), length(Eplotvals));
	%ds = d*ones(length(Eplotvals));
	Ns = N*ones(length(Dplotvals), length(Eplotvals));
	Eplotvalsd = repmat(Eplotvals, [length(Dplotvals), 1]);
	Dplotvalsd = repmat(Dplotvals', [1, length(Eplotvals)]);

	iphnumvals = arrayfun(@ssiphnum, kappas, gs, Dplotvalsd, Ns, Eplotvalsd);

	toc
end