function [D, Q] = buildD4(N)

	d4=[(-4751)/5192, 909/1298, 6091/15576, (-1165)/5192, 129/2596, (-25)/15576];
	q4=[649/576, 143/192, 75/64, 551/576];
	s4=[1/24,-9/8,9/8,-1/24];

	% construct 4th order DIV matrix : (N+1) -> N
	r4 = sparse(ones(1,4),1:4,s4,1,N+1);
	c4 = sparse(1,1,s4(1),N-2,1);
	D = [zeros(1,N+1);
	toeplitz(c4,r4);
	zeros(1,N+1)];
	[dm,dn]=size(d4);
	D(1:dm,1:dn)=d4;
	D(N-dm+1:N,N+1-dn+1:N+1)=-flipud(fliplr(d4));
	lnq=length(q4);
	qdiag=ones(1,N);
	qdiag(1:lnq)=q4; 
	qdiag(N:-1:(N-(lnq-1)))=q4;
	Q=sparse(diag(qdiag));
endfunction
