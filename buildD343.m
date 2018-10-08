function [D, Q] = buildD343(N)
	d343=[(-199)/208, 271/312, 7/52, (-5)/104, 1/624];
	q343=[13/12, 7/8, 25/24];
	s4=[1/24,-9/8,9/8,-1/24];
	% construct 4th order DIV matrix : (N+1) -> N
	r4 = sparse(ones(1,4),1:4,s4,1,N+1);
	c4 = sparse(1,1,s4(1),N-2,1);

	lnd=length(d343);
	rd343 = sparse(ones(1,lnd),1:lnd,d343,1,N+1);
	D = [rd343;toeplitz(c4,r4);-rd343(end:-1:1)];
	lnq=length(q343);
	qdiag=ones(1,N);
	qdiag(1:lnq)=q343; qdiag(N:-1:(N-(lnq-1)))=q343;
	Q=sparse(diag(qdiag));
endfunction
