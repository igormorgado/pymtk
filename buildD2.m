% 1D 2nd order centered finite difference operator
% from n+1 nodes to n cell-centers
% D : n x (n+1)

function hD=buildD2(n)

	%--------initialize D
	c=zeros(n,1); 
	c(1)=-1;
	r=zeros(1,n+1); 
	r(1)=-1; 
	r(2)=1;
	hD=toeplitz(c,r);
endfunction
