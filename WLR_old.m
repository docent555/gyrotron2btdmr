function w = WLR(sgn, cuEND, cuENDm1, cu_pEND, cu_pENDm1, W_PART, KAPPAJm1, KAPPAJ, coeff_dz_m_coeff_1i_d_6)

w = sgn * coeff_dz_m_coeff_1i_d_6/KAPPAJ*(2.0D0 * KAPPAJ * cu_pEND + ...
    2.0D0 * KAPPAJm1 * cuEND + KAPPAJ * cu_pENDm1 + KAPPAJm1 * cuENDm1) + W_PART;

% WL(IDX(step)) =  coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_p(1)  + 2.0D0 * cu(1)  + cu_p(2)    + cu(2))    + WL_PART;
% WR(IDX(step)) = -coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_p(Nz) + 2.0D0 * cu(Nz) + cu_p(Nzm1) + cu(Nzm1)) + WR_PART;

end