function [iphnum, iphnums] = timedependentsoln(E, detuning, tlist, varargin) 

	% set default values
	try 
		N = varargin{1};
	catch ME
		N = 10;
	end
	try 
		g = varargin{2};
	catch ME
		g = 0.5;
	end
	try 
		kappa = varargin{3}
	catch ME
		kappa = 0.01;
	end

	%generate identities and constants
	ida = identity(N); 
	idatom = identity(2);

	% Define cavity field and atomic operators
	a = tensor(destroy(N),idatom);
	sm = tensor(ida,sigmam);


	% Hamiltonian cavity and qubit in resonance, detuning- cavity-qubit freq minus E freq
	HJC = -detuning*(sm' *sm+ a'*a) + g*(sm'*a + a'*sm) + E*(a'+a); 
	
	% system comes to steady state at modified kappa
	% generate initial state
	% Collapse operators, spont emission now commented out
	C1 = sqrt(2*kappa)*a;
	%C2 = sqrt(gamma)*sm;
	C1dC1 = C1'*C1;
	%C2dC2 = C2'*C2;	
	
	% BuildL
	LH = -1i * (spre(HJC) - spost(HJC));
	L1 = spre(C1)*spost(C1')-0.5*spre(C1dC1)-0.5*spost(C1dC1); 
	% L2 = spre(C2)*spost(C2'')-0.5*spre(C2dC2)-0.5*spost(C2dC2); 
	L = LH+L1; %+L2;


	rhos = steady(L);


	% explicit initial state
    psi0 = tensor(basis(N,1),basis(2,2));
    rho0 = psi0 * psi0';


	% time evolution from steady state to state w/ default kappa
	ES = ode2es(L, rho0);
	
	iphnum = esval(expect(C1dC1, ES)/(2*kappa), tlist);
	iphnums = expect(C1dC1, rhos)/(2*kappa);
