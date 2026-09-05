Symbols m, Aplus, Aminus, Bplus, Bminus, E, F, G, H, Delta;
Vectors p, P, q, qplus, qminus, k;
Indices mu, nu;
Off statistics;

Local TrEE = g5_(1)
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrEF = g5_(1)
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*g_(1,P) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrEG = g5_(1)
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*q.P*g_(1,q) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrEH = g5_(1)
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*(g_(1,q,P) - g_(1,P,q))
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrFE = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrFF = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*g_(1,P) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrFG = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*q.P*g_(1,q) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrFH = g5_(1)*(p.p*g_(1,P)-p.P*g_(1,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*(g_(1,q,P) - g_(1,P,q))
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrGE = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrGF = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*g_(1,P) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrGG = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*q.P*g_(1,q) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrGH = g5_(1)*(P.P/(p.P)*g_(1,p)-g_(1,P))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*(g_(1,q,P) - g_(1,P,q))
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrHE = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrHF = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*g_(1,P) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrHG = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*( - i_*q.P*g_(1,q) )
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

Local TrHH = g5_(1)*(g_(1,p,P)-g_(1,P,p))
           *g_(1,mu)*(-i_*Aplus*g_(1,qplus) + Bplus*gi_(1))
           *g5_(1)*(g_(1,q,P) - g_(1,P,q))
           *(-i_*Aminus*g_(1,qminus) + Bminus*gi_(1))*g_(1,nu);

trace4,1;
.sort
contract;
id p.p*P.P = p.P^2 - Delta;

.sort
Local EE = TrEE * (d_(mu,nu) - k(mu)*k(nu)/k.k)/4;
Local EF = TrEF * (d_(mu,nu) - k(mu)*k(nu)/k.k)/4;
Local EG = TrEG * (d_(mu,nu) - k(mu)*k(nu)/k.k)/4;
Local EH = TrEH * (d_(mu,nu) - k(mu)*k(nu)/k.k)/4;

Local FE = TrFE * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local FF = TrFF * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local FG = TrFG * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local FH = TrFH * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);

Local GE = TrGE * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local GF = TrGF * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local GG = TrGG * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);
Local GH = TrGH * i_*(d_(mu,nu) - k(mu)*k(nu)/k.k)/(4*Delta);

Local HE = TrHE * (d_(mu,nu) - k(mu)*k(nu)/k.k)/(16*Delta);
Local HF = TrHF * (d_(mu,nu) - k(mu)*k(nu)/k.k)/(16*Delta);
Local HG = TrHG * (d_(mu,nu) - k(mu)*k(nu)/k.k)/(16*Delta);
Local HH = TrHH * (d_(mu,nu) - k(mu)*k(nu)/k.k)/(16*Delta);

.sort
contract;
.sort
#write <results.txt> "%e", EE
#write <results.txt> "%e", EF
#write <results.txt> "%e", EG
#write <results.txt> "%e", EH
#write <results.txt> "%e", FE
#write <results.txt> "%e", FF
#write <results.txt> "%e", FG
#write <results.txt> "%e", FH
#write <results.txt> "%e", GE
#write <results.txt> "%e", GF
#write <results.txt> "%e", GG
#write <results.txt> "%e", GH
#write <results.txt> "%e", HE
#write <results.txt> "%e", HF
#write <results.txt> "%e", HG
#write <results.txt> "%e", HH
Print +s;
.end