function [OUTF, OUTJ, p, Eff, Omega, ConLow, jout] = gyroscr(Nz, Nzi, Nt, Ne, ZAxis, ZAxisi, TAxis, TauAxis, Delta, ...
    Ib, ag, ai, apar, aperp, kappa, m, F, W, dt, dz, dzi, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField) %#codegen

if size(ZAxis, 2) ~= 1
    ZAxis = ZAxis.';
end
if size(ZAxisi, 2) ~= 1
    ZAxisi = ZAxisi.';
end
if size(kappa, 1) ~= 1
    kappa = kappa.';
end

s = 1;
% E = complex(zeros(1,Nt));
PHL = complex(zeros(1,Nt));
PHR = complex(zeros(1,Nt));
W0 = complex(zeros(1,Nt));
W1 = complex(zeros(1,Nt));
WNz = complex(zeros(1,Nt));
WNzm1 = complex(zeros(1,Nt));
WR = complex(zeros(1,Nt));
FNz = complex(zeros(1,Nt));
FNzm1 = complex(zeros(1,Nt));
BNz = complex(zeros(1,Nt));
BNzm1 = complex(zeros(1,Nt));
SigmaNz = complex(zeros(1,Nt));
SigmaNzm1 = complex(zeros(1,Nt));
WL = complex(zeros(1,Nt));
F0 = complex(zeros(1,Nt));
F1 = complex(zeros(1,Nt));
B0 = complex(zeros(1,Nt));
B1 = complex(zeros(1,Nt));
Sigma0 = complex(zeros(1,Nt));
Sigma1 = complex(zeros(1,Nt));
fmax = zeros(Nt, 1);
Jmax = zeros(Nt, 1);
field = complex(zeros(Nz,1));
testfield = complex(zeros(Nz,1));
field_p = complex(zeros(Nz,1));
rfield_p = complex(zeros(Nz,1));
lfield_p = complex(zeros(Nz,1));
B_p = complex(zeros(Nz,1));
B = complex(zeros(Nz,1));
p_p = complex(zeros(Nzi,Ne));
Jref = complex(zeros(Nzi,1));
J_p = complex(zeros(Nz,1));
J = complex(zeros(Nz,1));
OUTF = complex(zeros(OUTNz, OUTNt));
OUTJ = complex(zeros(OUTNz, OUTNt));
Eff = zeros(1,Nt);
Omega = zeros(1,Nt);
ConLow = zeros(1,Nt);
ZAxis_ipart = ZAxis(ZAxis <= ZAxisi(end));
ZAxis_ipart_end = length(ZAxis_ipart);
p = zeros(Nzi, Ne);
% pv = zeros(NzI, 2*Ne);
% p0 = zeros(Ne,1);
% p0v = zeros(2*Ne,1);
% reidx = zeros(1,Ne);
% imidx = zeros(1,Ne);

if INTZ > 1
    IZ = 0:INTZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTZ:length(ZAxis);
end

% kpar2 = zeros(length(ZAxis),1);
% N = length(ZAxis);
AA = complex(zeros(Nz,1));
BB = complex(zeros(Nz-1,1));
CC = complex(zeros(Nz-1,1));
DD = complex(zeros(Nz,1));

% SQR2 = sqrt(2);
SQRT2M2 = 2.828427124746190;
SQRT2D2 = 0.707106781186548;
SQRTDT = sqrt(dt);
SQRDZ = dz*dz;
SQRT2M2minus2p5 = SQRT2M2 - 2.5;

h = dzi;
Nzm1 = Nz - 1;
C0 = s;
% CL = -1i*kpar2(1);
% CR = -1i*kpar2(Nz);
CL = -1i;
CR = -1i;
C1 = 1/sqrt(1i*pi*s);
C2 = 1/sqrt(1i*pi*s);

E = 1i*C0/3*dz/dt./kappa;

W0 =     (-2*E - 1/dz);
W1 =     (  -E + 1/dz);
WNz =   -(-2*E - 1/dz);
WNzm1 = -(  -E + 1/dz);

% Initial values
jout = 1;
field(:,1) = InitialField;
% FforP = SF(ZAxisi(1:Nzi-1));
OUTF(:, jout) = field(IZ,1);
th0 = 2*pi*(0:Ne-1)/Ne;
p0 = exp(1i*th0);
p0v = [real(p0) imag(p0)].';
reidx = 1:Ne;
imidx = Ne+1:2*Ne;

% optsp = odeset('RelTol',1e-8,'AbsTol',1e-10);
% SF = griddedInterpolant(ZAxis, field, 'spline');
% p = oscill_ode_td(SF, Nzi, ZAxisi, apar(1), aperp(1), ag(1), m, s, p0v, reidx, imidx, optsp);
p = repmat(p0,Nzi,1);
% Jref(:) = Ib(1) * sum(p, 2)/Ne;
% Jref(:) = Current_o(p, Ib(1), Ne);
Jref(:) = Current_n(p_p, ai(IDX(0)), aperp(IDX(0)), s, apar(IDX(0)), Ib(IDX(0)), Ne);
SJ = griddedInterpolant(ZAxisi, Jref,'spline');
J(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);
B(:,1) = J(:) - 1i*kpar2(:).*field(:);
OUTJ(:,jout) = J(IZ,1);

fmax(IDX(0)) = max(abs(field(:,1)));
Jmax(IDX(0)) = max(abs(J(:,1)));
F0(IDX(0)) = field(1);
F1(IDX(0)) = field(2);
B0(IDX(0)) = B(1);
B1(IDX(0)) = B(2);
Sigma0(IDX(0)) = 0;
Sigma1(IDX(0)) = 0;
FNz(IDX(0)) = field(Nz);
FNzm1(IDX(0)) = field(Nzm1);
BNz(IDX(0)) = B(Nz);
BNzm1(IDX(0)) = B(Nzm1);
SigmaNz(IDX(0)) = 0;
SigmaNzm1(IDX(0)) = 0;
Eff(IDX(0)) = 1 - sum(abs(p(Nzi,:)).^2)/Ne;
Omega(IDX(0)) = 0;
PHL(IDX(0)) = 0;
PHR(IDX(0)) = 0;

WL(IDX(0)) = -dz * (-1i/6 * (2 * B0(IDX(0)) + B1(IDX(0))));
WR(IDX(0)) =  dz * (-1i/6 * (2 * BNz(IDX(0)) + BNzm1(IDX(0))));

SHOW = 1;
if SHOW == 1
    [lhfmax, lhfabs, lhjmax, lhjabs, hFig] = makeFig(ZAxis, ZAxis_ipart, TauAxis);
end

%Coefficients
% coeff_1i_m_C0_m_2_d_3_d_dt = 1i*C0*2/3/dt;
coeff_1i_m_C0_d_3_d_dt = 1i*C0/3/dt;
coeff_1i_d_6 = 1i/6;
coeff_4_d_3_m_SQRDT = 4/3*SQRTDT;
coeff_2_d_3_m_SQRDT = 2/3*SQRTDT;
coeff_1i_m_SQRDZ = 1i*SQRDZ;
coeff_1i_m_C0_m_SQRDZ_d_dt = 1i*C0*SQRDZ/dt;
% coeff_C1_m_coeff_4_d_3_m_SQRDT = C1*coeff_4_d_3_m_SQRDT;
% coeff_C2_m_coeff_4_d_3_m_SQRDT = C2*coeff_4_d_3_m_SQRDT;
% coeff_exp_CL_m_dt = exp(CL*dt);
% coeff_CL_m_dt = CL*dt;
% coeff_exp_CR_m_dt = exp(CR*dt);
% coeff_CR_m_dt = CR*dt;
coeff_dz_m_coeff_1i_d_6 = dz*coeff_1i_d_6;

SQRTKAPPA = sqrt(kappa);
R = kappa(2:Nt)./kappa(1:Nt - 1);

wp  = @(sgn, fEND, fENDm1, cuEND, cuENDm1, SigmaEND, SigmaENDm1, R, KAPPAJ) W_PART(sgn, fEND, fENDm1, cuEND, cuENDm1, ...
    SigmaEND, SigmaENDm1, R, KAPPAJ, dz, coeff_1i_m_C0_d_3_d_dt,coeff_dz_m_coeff_1i_d_6);
w   = @(sgn, cu_pEND, cu_pENDm1, W_PART) WLR(sgn, ...
    cu_pEND, cu_pENDm1, W_PART,coeff_dz_m_coeff_1i_d_6);
ilr = @(step, u) ILR(step, u, coeff_4_d_3_m_SQRDT, SQRT2D2, SQRT2M2minus2p5);
dp  = @(IL, IR, field, B,  WL, WR, W0, W1, WNzm1, WNz, PL, PR, SQRKAPPA, KAPPA) DAR_PART(IL, IR, field, B, WL, WR, W0, W1, WNzm1, WNz, PL, PR, SQRKAPPA, KAPPA, ...
    coeff_2_d_3_m_SQRDT, coeff_1i_m_C0_m_SQRDZ_d_dt,coeff_1i_m_SQRDZ, C1, C2, CL, CR, Nz);
d   = @(B_p, WL, WR, DPART, SQRKAPPA, KAPPAJ) DAR(B_p, WL, WR, DPART, SQRKAPPA, KAPPAJ, Nz, C1, C2, coeff_4_d_3_m_SQRDT, coeff_1i_m_SQRDZ);
sgm = @(f1, f2, j1, j2, Sold, R, KAPPAJ) Sigma(f1, f2, j1, j2, Sold, R, KAPPAJ, coeff_1i_m_C0_d_3_d_dt, coeff_1i_d_6);

num_st_test_iter = 0;
fmax_glob_old = max(abs(field(:,1)));

fprintf('\n');
% timerVal = tic;
for step=1:Nt-1

    if SHOW == 1
        lhfmax.YData(1:step) = Omega(1:step);
        lhfmax.XData(1:step) = TauAxis(1:step);
        lhfabs.YData = abs(field);

        lhjmax.YData(1:step) = Eff(1:step);
        lhjmax.XData(1:step) = TauAxis(1:step);
        lhjabs.YData = abs(J(1:ZAxis_ipart_end));

        drawnow
    end     

    PL = h*kpar2(1)*kappa(IDX(step));
    PR = h*kpar2(Nz)*kappa(IDX(step));
    PHL(IDX(step)) = PHL(IDX(step-1)) + PL;
    PHR(IDX(step)) = PHR(IDX(step-1)) + PR;

    AA(1) =  1 - 4/3*C1*W0(IDX(step))*SQRTDT*SQRTKAPPA(IDX(step));
    AA(2:Nzm1) = -2*(kappa(IDX(step)) + 1i * SQRDZ/dt*C0);
    AA(Nz) = 1 + 4/3*C2*WNz(IDX(step))*SQRTDT*SQRTKAPPA(IDX(step));
    BB(1) =   -4/3*C1*W1(IDX(step))*SQRTDT*SQRTKAPPA(IDX(step));
    BB(2:Nzm1) = kappa(IDX(step));
    CC(1:Nzm1 - 1) = kappa(IDX(step));
    CC(Nzm1) = 4/3*C2*WNzm1(IDX(step))*SQRTDT*SQRTKAPPA(IDX(step));  

    M = spdiags([[CC; 0] AA [0 ;BB]], -1:1, Nz, Nz);

    WL_PART = wp(-1, field(1),  field(2), B(1),  B(2),    Sigma0(IDX(step-1)),  Sigma1(IDX(step-1)), R(IDX(step-1)), kappa(IDX(step)));
    WR_PART = wp( 1, field(Nz), field(Nzm1),B(Nz), B(Nzm1),  SigmaNz(IDX(step-1)), SigmaNzm1(IDX(step-1)), R(IDX(step-1)), kappa(IDX(step)));   

    IL = ilr(step, @ul);
    IR = ilr(step, @ur);
    
    D_PART = dp(IL, IR, field,B(2:Nzm1), WL(IDX(step-1)), WR(IDX(step-1)), W0(IDX(step-1)), ...
        W1(IDX(step-1)), WNzm1(IDX(step-1)), WNz(IDX(step-1)), PL, PR, SQRTKAPPA(IDX(step-1)), kappa(IDX(step-1)));

    % signal
    D_PART = D_PART + 2*F*exp(1i*W*TAxis(IDX(step))*kappa(IDX(step)) + 1i*PHR(IDX(step)));
    
    SF = griddedInterpolant(ZAxis, field,'spline');
    Fref = SF(ZAxisi);    
    % RHS0 = -Fref - 1i*p.*(Delta - 1 + abs(p).^2);
    % RHS0 = RHS_o(p, Fref, Delta);
    phi = 0; % временно (выяснить что такое phi)
    RHS0  = RHS_n(p, Fref, aperp(IDX(step-1)), ag(IDX(step-1)), apar(IDX(step-1)), s, m, phi);
    % RHS0 = 1i*aperp(IDX(step-1))^(s-2)./(ag(IDX(step-1))*apar(IDX(step-1)))*Fref.*conj(p).^(s-1).*exp(-1i*(m-s)*phi) ...
    % - 1i*aperp(IDX(step-1))^2/apar(IDX(step-1))*p.*(abs(p).^2 - 1);
    p_p = [p0; p(1:Nzi-1,:) + h * RHS0(1:Nzi-1,:)]; % так наверное
    % правильее (производная в начальной точке на шаг - метод Эйлера)
    % p_p = [p0; p(2:Nzi,:) + h * RHS0(2:Nzi,:)]; % вообще не правильно
    % (хотя считает, все равно)
    % p_p = p; % в производная в следующей по пространству точке берется
    % начиная с первого шага (внутри цикла)
      
    % Jref(:) = Ib(IDX(step)) * sum(p_p, 2)/Ne;    
    % Jref(:) = Current_o(p_p, Ib(IDX(step)), Ne);
    % Jref(:) = 1i*ai(IDX(step))*aperp(IDX(step))^s/apar(IDX(step))*Ib(IDX(step)) * sum(p_p, 2)/Ne;
    Jref(:) = Current_n(p_p, ai(IDX(step)), aperp(IDX(step)), s, apar(IDX(step)), Ib(IDX(step)), Ne);

    SJ = griddedInterpolant(ZAxisi, Jref,'spline');
    J_p(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);

    B_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);

    WL(IDX(step)) = w( 1, B_p(1),  B_p(2),    WL_PART);
    WR(IDX(step)) = w(-1, B_p(Nz), B_p(Nzm1), WR_PART);
    DD = d(B_p(2:Nzm1), WL(IDX(step)), WR(IDX(step)), D_PART, SQRTKAPPA(IDX(step)), kappa(IDX(step)));    
    
    field_p(:,1) = M \ DD; % для 1-го приближени методом Эйлера
    % rfield_p(:,1) = rtridag(CC,AA,BB,DD);
    % lfield_p(:,1) = ltridag(CC,AA,BB,DD);
    % field_p = (rfield_p + lfield_p)/2;
        
    num_insteps = 0;
    maxfield = max(abs(field_p(:,1)));
    testfield = field_p;
    % testfield = 100*ones(size(field_p));
    while 1
        num_insteps = num_insteps + 1;
        
        SF = griddedInterpolant(ZAxis, field_p,'spline');
        Fref_p = SF(ZAxisi);

        % RHS1 = -Fref_p - 1i*p_p.*(Delta - 1 + abs(p_p).^2);    
        % RHS1 = RHS_o(p_p, Fref_p, Delta);
        % RHS1 = 1i*aperp(IDX(step))^(s-2)./(ag(IDX(step))*apar(IDX(step)))*Fref_p.*conj(p_p).^(s-1).*exp(-1i*(m-s)*phi) ...
        %     - 1i*aperp(IDX(step))^2/apar(IDX(step))*p_p.*(abs(p_p).^2 - 1);
        RHS1 = RHS_n(p_p, Fref_p, aperp(IDX(step)), ag(IDX(step)), apar(IDX(step)), s, m, phi);

        p_p = [p0; p(1:Nzi-1,:) + h/2 * (RHS0(1:Nzi-1,:) + RHS1(2:Nzi,:))];        
        % Jref(:) = Ib(IDX(step)) * sum(p_p, 2)/Ne;
        % Jref(:) = Current_o(p_p, Ib(IDX(step)), Ne);
        Jref(:) = Current_n(p_p, ai(IDX(step)), aperp(IDX(step)), s, apar(IDX(step)), Ib(IDX(step)), Ne);
        SJ = griddedInterpolant(ZAxisi, Jref,'spline');
        J_p(1:ZAxis_ipart_end,1) = SJ(ZAxis_ipart);

        B_p(:,1) = J_p(:) - 1i*kpar2(:).*field_p(:);

        WL(IDX(step)) = w( 1, B_p(1),  B_p(2),    WL_PART);
        WR(IDX(step)) = w(-1, B_p(Nz), B_p(Nzm1), WR_PART);        

        DD = d(B_p(2:Nzm1), WL(IDX(step)), WR(IDX(step)), D_PART, SQRTKAPPA(IDX(step)), kappa(IDX(step)));
        
        % samosoglasovannoe pole
        field_p(:,1) = M \ DD;
        % rfield_p(:,1) = rtridag(CC,AA,BB,DD);
        % lfield_p(:,1) = ltridag(CC,AA,BB,DD);
        % field_p = (rfield_p + lfield_p)/2;

        maxdiff = max(abs(testfield - field_p));
        err = maxdiff/maxfield;
        if err < tol
            field(:,1) = field_p(:,1);
            J(:,1) = J_p(:,1);
            p(:,:) = p_p(:,:);
            B(:,1) = B_p(:,1);
            break
        end
        testfield = field_p;
        if num_insteps > 1000
            fprintf('\nToo many inner steps!\n');
            pause;
        end
    end
  
    fmax(IDX(step)) = max(abs(field(:,1)));
    Jmax(IDX(step)) = max(abs(J(:,1)));

    F0(IDX(step)) =  field(1);
    F1(IDX(step)) = field(2);
    FNz(IDX(step)) =  field(Nz);
    FNzm1(IDX(step)) = field(Nzm1);
    B0(IDX(step)) = B(1);
    B1(IDX(step)) = B(2);
    BNz(IDX(step)) = B(Nz);
    BNzm1(IDX(step)) = B(Nzm1);

    %     Omega(IDX(step)) = (angle(field(Nz)) - angle(FNz(IDX(step-1))))/dt;
    Omega(IDX(step)) = imag(log(field(Nz)/FNz(IDX(step-1))))/dt;
    Eff(IDX(step)) = 1 -  sum(abs(p(Nzi,:)).^2)/Ne;

    if (mod(num_st_test_iter,1000))
        fmax_glob_new = max(abs(field(:,1)));
        if abs(fmax_glob_new - fmax_glob_old)/fmax_glob_old < tol
            jout = jout + 1;
            OUTF(:, jout) = field(IZ,1);
            OUTJ(:, jout) = J(IZ,1);
            fprintf('Emergency exit!\n');
            return;
        end
        num_st_test_iter = num_st_test_iter + 1;
        fmax_glob_old = fmax_glob_new;
    end

    if mod(step,INTT) == 0
        jout = jout + 1;
        OUTF(:, jout) = field(IZ,1);
        OUTJ(:, jout) = J(IZ,1);
    end

    Sigma0(IDX(step))    = sgm(field(1),    F0(IDX(step - 1)),    B(1),     B0(IDX(step - 1)),    Sigma0(IDX(step - 1)), R(IDX(step - 1)), kappa(IDX(step)));    
    Sigma1(IDX(step))    = sgm(field(2),    F1(IDX(step - 1)),    B(2),     B1(IDX(step - 1)),    Sigma1(IDX(step - 1)), R(IDX(step - 1)), kappa(IDX(step)));    
    SigmaNz(IDX(step))   = sgm(field(Nz),   FNz(IDX(step - 1)),   B(Nz),    BNz(IDX(step - 1)),   SigmaNz(IDX(step - 1)), R(IDX(step - 1)), kappa(IDX(step)));
    SigmaNzm1(IDX(step)) = sgm(field(Nzm1), FNzm1(IDX(step - 1)), B(Nzm1),  BNzm1(IDX(step - 1)), SigmaNzm1(IDX(step - 1)), R(IDX(step - 1)), kappa(IDX(step)));    

    k = step + 1;

    ConLow(IDX(step)) = (2*imag(field(Nz)*conj(dfdzNz(step)) - field(1)*conj(dfdz0(step))) - Ib(IDX(step))*Eff(IDX(step)))/(Ib(IDX(step))*Eff(IDX(step)))*100;

    zabey = repmat('\b', 1, 162);
    fprintf(zabey + ...
        "Step = %8i   Time = %10.4f   Fmax = % 15.10e   Jmax = % 15.10e   W = % 15.10e   E = % 15.10e   CL = % 15.5f %%",...
        int64(step), TAxis(k), fmax(k), Jmax(k), Omega(IDX(step)), Eff(IDX(step)), ConLow(IDX(step)));
end

OUTJ(:,jout) = J(IZ,1);

fprintf("\n\n\n");

% ExecutionTime = toc(timerVal);

% hours = fix(ExecutionTime/3600);
% minutes = fix((ExecutionTime - fix(hours*3600))/60);
% seconds = ExecutionTime - hours*3600 - minutes*60;
% 
% fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

fprintf(" \n\n");

    function  f = ur(J)
        coder.inline("always");
        I = IDX(J);       
        f = ((WNzm1(I) .* FNzm1(I) + WNz(I) .* FNz(I) + WR(I)) .* SQRTKAPPA(I)) .* exp(CR .* (PHR(IDX(step)) - PHR(I)));        
    end

    function  f = ul(J)
        coder.inline("always");
        I = IDX(J);                
        f = ((W1(I)   .* F1(I)     + W0(I) .* F0(I)   + WL(I)) .* SQRTKAPPA(I)) .* exp(CL .* (PHL(IDX(step)) - PHL(I)));       
    end

    function j = IDX(j)
        coder.inline("always");
        j = j + 1;
    end

    function  dfdz = dfdzNz(j)
        coder.inline("always");
        dfdz = (WNzm1(IDX(j)) .* FNzm1(IDX(j)) + WNz(IDX(j)) .* FNz(IDX(j)) + WR(IDX(j)));
    end

    function  dfdz = dfdz0(j)
        coder.inline("always");
        dfdz = (W1(IDX(j)) .* F1(IDX(j)) + W0(IDX(j)) .* F0(IDX(j)) + WL(IDX(j)));
    end
end




