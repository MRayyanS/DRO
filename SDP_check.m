cvx_begin sdp

N = 10;
e = (1/N).*ones(N,1);


variable eta
variable lam(N,1)

J = eta^2 - e'*lam;

minimize (J)

subject to
 eta >= 0 
 lam >= 0

 for i = 1:N
     eta.*eye(100) - lam(i).*eye(100) >= 0
 end

cvx_end


eta
lam



