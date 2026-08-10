Symbols m, k2;
Vectors p, q, k;
Indices mu, nu;

Local Tr1 = g_(1,mu,mu);

Local Tr2 = g_(1,k,k)/k2;

trace4,1;
.sort

contract;
.sort

* k^2 als Symbol einsetzen (falls k.k noch nicht durch k2 ersetzt wurde)
id k.k = k2;

.sort
Print Tr1, Tr2;
.end

* Results are: 
* Tr1 = 16;
* Tr2 = 4;
* Fuck, Metrik