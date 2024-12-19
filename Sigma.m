function S = Sigma(f1, f2, j1, j2, Sold, R, KAPPAJ, coeff_1i_m_C0_d_3_d_dt, coeff_1i_d_6)

S = coeff_1i_m_C0_d_3_d_dt/KAPPAJ * (f1 - f2) ...
    -coeff_1i_d_6*(j1 + R * j2) - R*Sold;

% S = -(-coeff_1i_m_C0_d_3_d_dt)/KAPPAJ * f1 ...
%     + (-coeff_1i_m_C0_d_3_d_dt)/KAPPAJ * f2 ...
%     -coeff_1i_d_6/KAPPAJ*(KAPPAJ * j1 + KAPPAJm1 * j2) - KAPPAJm1/KAPPAJ*Sold;

end