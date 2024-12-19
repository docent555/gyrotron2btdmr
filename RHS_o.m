function res = RHS_o(p, Fref, Delta)

res = -Fref - 1i*p.*(Delta - 1 + abs(p).^2);

end