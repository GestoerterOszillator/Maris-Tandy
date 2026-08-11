Symbols m, k2;
Vectors p, q, k;
Indices mu, nu;

Local Tr = (g_(1,mu,mu) - g_(1,k,k)/k.k)/4;

trace4,1;
.sort

contract;

.sort
Print Tr;
.end

* Results are: 
* Tr = 3;