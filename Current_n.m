function c = Current_n(p, ai, aperp, sn, apar, I, Ne)

c = 1i*ai*aperp^sn/apar*I * sum(p, 2)/Ne;

end