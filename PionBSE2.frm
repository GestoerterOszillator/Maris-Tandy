Symbols n, Aplus, Aminus, Bplus, Bminus, Delta, E, F, G, H, p2, P2, q2, k2, pP, pq, Pq, pk, qk, Pk;
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
id p.p*P.P = (p.P)^2 - Delta;
* Delta = pP^2 - p2*P2
.sort
id qplus = q + P/2;
id qminus = q - P/2;

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
id p.p = p2;
id P.P = P2;
id k.k^n? = k2^n;
id q.q = q2;
id p.P^n? = pP^n;
id p.q = pq;
id P.q = Pq;
id p.k = pk;
id q.k = qk;
id P.k = Pk;

.sort
#write <Pion.jl> "function EE(p2, q2, P2, pP, pq, Pq, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", EE
#write <Pion.jl> "end\n"

#write <Pion.jl> "function EF(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", EF
#write <Pion.jl> "end\n"

#write <Pion.jl> "function EG(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", EG
#write <Pion.jl> "end\n"

#write <Pion.jl> "function EH(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", EH
#write <Pion.jl> "end\n"

#write <Pion.jl> "function FE(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", FE
#write <Pion.jl> "end\n"

#write <Pion.jl> "function FF(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", FF
#write <Pion.jl> "end\n"

#write <Pion.jl> "function FG(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", FG
#write <Pion.jl> "end\n"

#write <Pion.jl> "function FH(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", FH
#write <Pion.jl> "end\n"

#write <Pion.jl> "function GE(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", GE
#write <Pion.jl> "end\n"

#write <Pion.jl> "function GF(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", GF
#write <Pion.jl> "end\n"

#write <Pion.jl> "function GG(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", GG
#write <Pion.jl> "end\n"

#write <Pion.jl> "function GH(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", GH
#write <Pion.jl> "end\n"

#write <Pion.jl> "function HE(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", HE
#write <Pion.jl> "end\n"

#write <Pion.jl> "function HF(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", HF
#write <Pion.jl> "end\n"

#write <Pion.jl> "function HG(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", HG
#write <Pion.jl> "end\n"

#write <Pion.jl> "function HH(p2, q2, P2, k2, pP, pq, Pq, pk, qk, Pk, A, B)"
#write <Pion.jl> "    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq - q2; Pk = pP - Pq; Delta = pP^2 - p2*P2"
#write <Pion.jl> "    return %E", HH
#write <Pion.jl> "end\n"

.end