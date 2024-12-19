function u = rhs_td(z, p, apar, aperp, ag, m, s, F)

u = 1i*aperp.^(s-2)./(ag.*apar).*F(z).*conj(p).^(s-1).*exp(-1i*(m-s)) ...
    - 1i*aperp.^2./apar.*p.*(abs(p).^2 - 1);

end
