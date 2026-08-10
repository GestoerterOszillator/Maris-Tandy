Symbols m, k2;
Vectors p, q, k;
Indices mu, nu;

* Teil 1: p-slash gamma_mu q-slash gamma_nu * delta^{mu nu}
* -> delta^{mu nu} kontrahiert einfach mu und nu miteinander
Local Tr1 = g_(1,p,mu,q,mu);

* Teil 2: p-slash gamma_mu q-slash gamma_nu * k^mu k^nu / k^2
Local Tr2 = g_(1,p,k,q,k)/k2;

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
* Tr1 = - 8*p.q;
* Tr2 = - 4*p.q + 8*p.k*q.k*k2^-1;
