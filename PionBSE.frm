Symbols m, Aplus, Aminus, Bplus, Bminus, E, F, G, H;
Vectors p, P, q, qplus, qminus, k;
Indices mu, nu;
Off statistics;

Local TrE = g5_(1)
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*q.P*g_(1,q)*G + (g_(1,q,P) - g_(1,P,q))*H )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrF = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*q.P*g_(1,q)*G + (g_(1,q,P) - g_(1,P,q))*H )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrG = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*q.P*g_(1,q)*G + (g_(1,q,P) - g_(1,P,q))*H )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrH = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*q.P*g_(1,q)*G + (g_(1,q,P) - g_(1,P,q))*H )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

trace4,1;
.sort

contract;

.sort
Local ResultE = TrE *    (d_(mu,nu) - k(mu)*k(nu)/k.k)/4;
Local ResultF = TrF * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*((p.P)^2 - p.p*P.P));
Local ResultG = TrG * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*((p.P)^2 - p.p*P.P));
Local ResultH = TrH *    (d_(mu,nu) - k(mu)*k(nu)/k.k)/(16*((p.P)^2 - p.p*P.P));

.sort
Print +s;
.end