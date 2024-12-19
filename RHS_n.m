function rhs = RHS_n(p, Fref, aperp, ag, apar, sn, m, phi)

rhs = 1i*aperp^(sn-2)./(ag*apar)*Fref.*conj(p).^(sn-1).*exp(-1i*(m-sn)*phi) ...
    - 1i*aperp^2/apar*p.*(abs(p).^2 - 1);

end