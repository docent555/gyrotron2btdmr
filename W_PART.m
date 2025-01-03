function u = W_PART(sgn, fEND, fENDm1, cuEND, cuENDm1, SigmaEND, SigmaENDm1, R, KAPPAJ, dz, ...
    coeff_1i_m_C0_d_3_d_dt, coeff_dz_m_coeff_1i_d_6)

u = sgn * dz * ((-coeff_1i_m_C0_d_3_d_dt/KAPPAJ) * (2*fEND + fENDm1) ...
    - R * (2.0D0 * SigmaEND + SigmaENDm1))...
-sgn * coeff_dz_m_coeff_1i_d_6 * (2.0D0 * R * cuEND + R * cuENDm1);

% WL_PART = -dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(1)...
%     + (-coeff_1i_m_C0_d_3_d_dt) * field(2)...
%     - (2.0D0 * Sigma0(IDX(step-1)) + Sigma1(IDX(step-1))));
% 
% 
% 
% WR_PART =  dz * ((-coeff_1i_m_C0_m_2_d_3_d_dt) * field(Nz)...
%     + (-coeff_1i_m_C0_d_3_d_dt) * field(Nzm1)...
%     - (2.0D0 * SigmaNz(IDX(step-1)) + SigmaNzm1(IDX(step-1))));
end