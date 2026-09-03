Symbols m, A_plus, A_minus, B_plus, B_minus, E;
Vectors p, P, q, q_plus, q_minus, k;
Indices mu, nu;

Local Tr = -(g_(1,p,mu,q,mu) - g_(1,p,k,q,k)/k.k)/4;

trace4,1;
.sort

contract;

.sort
Print Tr;
.end