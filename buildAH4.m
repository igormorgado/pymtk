function D = buildAH4(M)
	% the matrix A_H from Yefet and Petropoulos, ’A Staggered Fourth-Order
	% Accurate Explicit Finiite Difference Scheme for the Time-Domain Maxwell’s Equations’,
	% Computational Physics, 168(2001), pp286-315

	N=M-1;
	d4=[-23/24, 21/24, 3/24 ,-1/24];
	s4=[1/24,-9/8,9/8,-1/24];

	% construct 4th order 1st derivative matrix : (N+1) -> N
	r4 = sparse(ones(1,4),1:4,s4,1,N+1);
	c4 = sparse(1,1,s4(1),N-2,1);
	lnd=length(d4);
	rd4 = sparse(ones(1,lnd),1:lnd,d4,1,N+1);
	D = [rd4;toeplitz(c4,r4);-rd4(end:-1:1)];
endfunction
