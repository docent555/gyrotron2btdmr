function w = WLR(sgn, cu_pEND, cu_pENDm1, W_PART, coeff_dz_m_coeff_1i_d_6)

w = sgn * coeff_dz_m_coeff_1i_d_6 * (2.0D0 * cu_pEND +  cu_pENDm1) + W_PART;

end