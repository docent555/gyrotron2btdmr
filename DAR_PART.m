function dp = DAR_PART(IL, IR, field, cu, WL, WR, W0, W1, WNzm1, WNz, PL, PR, SQRKAPPAJm1, KAPPAJm1, ...
    coeff_2_d_3_m_SQRDT, coeff_1i_m_C0_m_SQRDZ_d_dt, coeff_1i_m_SQRDZ, C1, C2, CL, CR, Nz)

Nzm1 = Nz - 1;

D_0_PART =   C1 * (IL + coeff_2_d_3_m_SQRDT * SQRKAPPAJm1 * (W0 * field(1)...
    + W1 * field(2) + WL) * exp(CL * PL));

D_MIDDLE_PART = 2.0D0 * (KAPPAJm1 - coeff_1i_m_C0_m_SQRDZ_d_dt) .* field(2:Nzm1)...
    - KAPPAJm1*(field(1:Nz - 2) + field(3:Nz)) -coeff_1i_m_SQRDZ * (KAPPAJm1*cu);

D_END_PART = - C2 * (IR + coeff_2_d_3_m_SQRDT * SQRKAPPAJm1 * (WNzm1 * field(Nzm1)...
    + WNz * field(Nz) + WR) * exp(CR * PR));

dp = [D_0_PART; D_MIDDLE_PART; D_END_PART];

end