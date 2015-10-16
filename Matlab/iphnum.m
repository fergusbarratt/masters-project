function [rhoss] = rhoss(E,kappa,g,det,N)
	
	% denity matrix for JC hamiltonian in steady state with conditions.

	% [count1, count2, iphnum] = probss(E,kappa,gamma,g,wc,w0,wl)
	% solves the problem of a coherently driven cavity with a two-level atom
	% E = amplitude of driving field, kappa = mirror coupling,
	% gamma = spontaneous emission rate, 
	% g = atom-field coupling,
	% wc = cavity frequency, 
	% w0 = atomic frequency, 
	% wl = driving field frequency, 
	% (together are det: detuning from cavity-qubit resonance)
	% N = size of Hilbert space for intracavity field (zero to N-1 photons)
	% count1 = photocount rate of light leaking out of cavity count2 = spontaneous emission rate
	% iphnum = intracavity field


	%generate identities and constants
	ida = identity(N); 
	idatom = identity(2);

	% Define cavity field and atomic operators
	a = tensor(destroy(N),idatom);
	sm = tensor(ida,sigmam);

	% Hamiltonian cavity and qubit in resonance, det- cavity-qubit freq minus drive freq
	HJC = -det*(sm' *sm+ a'*a) + g*(sm'*a + a'*sm) + E*(a'+a); 
	
	% Collapse operators, spont emission now commented out
	C1 = sqrt(2*kappa)*a;
	%C2 = sqrt(gamma)*sm;
	C1dC1 = C1'*C1;
	%C2dC2 = C2'*C2;
	
	% Calculate the Liouvillian, spont emission lindbladians gone
	LH = -1i * (spre(HJC) - spost(HJC));
	L1 = spre(C1)*spost(C1')-0.5*spre(C1dC1)-0.5*spost(C1dC1); 
	% L2 = spre(C2)*spost(C2'')-0.5*spre(C2dC2)-0.5*spost(C2dC2); 
	L = LH+L1; %+L2;
	
	% Find steady state
	rhoss = steady(L);
	
	% Calculate expectation values expectation of collapse conj*collapse
	% Monitor output field of the cavity. count1: output photon count due to kappa
	% count1 = expect(C1dC1,rhoss);
	% count2 = expect(C2dC2,rhoss). infield is now intracavity photon number iphnum
	% iphnum = real(count1/(2*kappa))
