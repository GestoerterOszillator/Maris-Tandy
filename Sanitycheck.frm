Symbols E, F, G, H, Delta;
Vectors p, P;
Off statistics;

Local TrE = g5_(1)
           *g5_(1)*( E - i_*g_(1,P)*F - i_*p.P*g_(1,p)*G + (g_(1,p,P) - g_(1,P,p))*H );

Local TrF = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*p.P*g_(1,p)*G + (g_(1,p,P) - g_(1,P,p))*H );

Local TrG = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*p.P*g_(1,p)*G + (g_(1,p,P) - g_(1,P,p))*H );

Local TrH = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g5_(1)*( E - i_*g_(1,P)*F - i_*p.P*g_(1,p)*G + (g_(1,p,P) - g_(1,P,p))*H );

trace4,1;
.sort
contract;
id p.p*P.P = p.P^2 - Delta;
.sort
Local ResultE = TrE/4;
Local ResultF = TrF * i_/(4*Delta);
Local ResultG = TrG * i_/(4*Delta);
Local ResultH = TrH /(16*Delta);
.sort
contract;
.sort
Print +s;
.end