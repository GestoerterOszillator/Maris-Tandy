Symbols m, k2;
Vectors p, q, k;
Indices mu, nu;

Local Tr = -(g_(1,p,mu,q,mu) - g_(1,p,k,q,k)/k.k)/4;

trace4,1;
.sort

contract;

.sort
Print Tr;
.end

* Results are: 
* Tr = p.q + 2*p.k*q.k*k.k^-1;