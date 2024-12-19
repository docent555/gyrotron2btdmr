function uv = rhsv_td(z, pv, apar, aperp, ag, m, s, a, reidx, imidx)

p = pv(reidx) + 1i*pv(imidx);

u = rhs_td(z, p, apar, aperp, ag, m, s, a);

uv = [real(u); imag(u)];

end