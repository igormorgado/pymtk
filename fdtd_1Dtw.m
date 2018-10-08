function [L2errEz, nmax] = fdtd_1Dtw(method)

	%1D traveling wave

	dx=1/40;
	dt=1/1600;
	Tmax=4;

	% number of timeteps
	nmax=ceil(Tmax/dt);
	% number of cells in the x-direction
	I=ceil(1/dx);

	Hy=zeros(I,1);
	Ez=zeros(I+1,1);
	errEz=zeros(1,nmax);
	L2errEz=zeros(1,nmax);

	if (strcmp(method, 'cg3434'))
		% CG:4/343
		Deh=buildD4(I);
		Dhe=buildD343(I-1); 
	elseif (strcmp(method,'explicit24'))
		% exact(2,4)
		Deh=buildAE4(I);
		Dhe=buildAH4(I);
	elseif (strcmp(method,'cg4'))
		% CG:4/4
		Deh=buildD4(I);
		Dhe=buildD4(I-1);
	elseif (strcmp(method,'cg343'))
		% CG:343/343
		Deh=buildD343(I);
		Dhe=buildD343(I-1); 
	elseif (strcmp(method,'yee222'))
		%standard FDTD
		Deh=buildD2(I);
		Dhe=buildD2(I-1); 
		dt=dx;
	end


	exactEz = @(x,t) sin(3*pi*(x+t));
	exactHy = @(x,t) sin(3*pi*(x+t));

	% -------------- n=0

	% Ez @ t = 0
	EzIC=exactEz((1:I-1)'*dx,0);

	% Hy @ t = (1/2)dt
	HyIC=exactHy(((0:I-1)'+0.5)*dx,0.5*dt);

	% -------------- n=1

	% Ez @ t=1*dt
	Ez(2:I)=EzIC+(dt/dx)*Dhe*HyIC;

	% BCs
	Ez(1)=exactEz(0,dt);
	Ez(I+1)=exactEz(1,dt);

	errEz=abs(Ez-exactEz((0:I)'*dx,dt));
	L2errEz(1)=sqrt(dx)*norm(errEz,2);

	% Hy @ t=(3/2)dt
	Hy=HyIC+(dt/dx)*Deh*Ez;

	% -------------- n=2:nmax
	for n=2:nmax
		% Ez @ t = n*dt
		Ez(2:I)=Ez(2:I)+(dt/dx)*(Dhe*Hy);
		% BCs
		Ez(1)=exactEz(0,n*dt);

		Ez(I+1)=exactEz(1,n*dt);
		errEz=abs(Ez-exactEz((0:I)'*dx,n*dt));
		L2errEz(n)=sqrt(dx)*norm(errEz,2);
		% Hy @ t = (n+1/2)dt
		Hy=Hy+(dt/dx)*Deh*Ez;
	end
endfunction
