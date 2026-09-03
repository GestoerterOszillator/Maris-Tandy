Symbols m, A_plus, A_minus, B_plus, B_minus, E, F, G, H;
Vectors p, P, q, q_plus, q_minus, k;
Indices mu, nu;

Local Tr1 = g5_(1)*g_(1,mu)*(-i_*A_plus*g_(1,q_plus)+B_plus*gi_(1))
           *( E*g5_(1) )
           *(-i_*A_minus*g_(1,q_minus)+B_minus*gi_(1))*g_(1,nu);

trace4,1;
.sort

contract;

.sort
Local Result1 = Tr1*(d_(mu,nu) - k(mu)*k(nu)/k.k)/4;

.sort
Print +s;
.end