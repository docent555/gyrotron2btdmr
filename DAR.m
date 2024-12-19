function d = DAR(cu_p, WL, WR, DPART, SQRTKAPPA, KAPPAJ, Nz, C1, C2, coeff_4_d_3_m_SQRDT, coeff_1i_m_SQRDZ) %#codegen

Nzm1 = Nz - 1;

d1  = DD( 1, WL, DPART(1), SQRTKAPPA,  C1, coeff_4_d_3_m_SQRDT);

d2_Nzm1 = -coeff_1i_m_SQRDZ * KAPPAJ*cu_p + DPART(2:Nzm1);

dNz = DD(-1, WR, DPART(Nz), SQRTKAPPA, C2, coeff_4_d_3_m_SQRDT);

d = [d1; d2_Nzm1; dNz];

% D(1)  =   coeff_C1_m_coeff_4_d_3_m_SQRDT * WL(IDX(step)) + D_0_PART;
% D(Nz) = - coeff_C2_m_coeff_4_d_3_m_SQRDT * WR(IDX(step)) + D_END_PART;

end

function dd = DD(sgn, W, DPART, SQRTKAPPA, C, coeff_4_d_3_m_SQRDT)
dd  = sgn * C * coeff_4_d_3_m_SQRDT * SQRTKAPPA * W + DPART;
end