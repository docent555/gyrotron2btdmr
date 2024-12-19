function [OUTFre, OUTFim, OUTJre, OUTJim, pre, pim, ConLow, OUTZAxis, OUTTAxis, Eff, Omega, jout] = gyrotron() %#codegen
% function [OUTFre, OUTFim, OUTJre, OUTJim, pre, pim, ConLow, OUTZAxis, OUTTAxis, Eff, Omega, jout] = ...
% gyrotron(Ne, z0, zex, zin, Tend, Delta, I0, dz, dt, tol, INTT, INTZ, a0, nexl, nexr, mode, akp2) %#codegen

input_param = read_namelist('input_fortran.in', 'param');

Ne = input_param.ne;
z0 = input_param.z0;
zex = input_param.zex;
zexi = input_param.zexi;
TauEnd = input_param.tauend;
Delta = input_param.delta;
Nz = input_param.nz;
Ntau = input_param.ntau;
Nzi = input_param.nzi;
tol = input_param.tol;
INTT = input_param.intt;
INTZ = input_param.intz;
a0 = input_param.a0;
nexl = input_param.nexl;
nexr = input_param.nexr;
e = input_param.e;
B = input_param.b;
m = input_param.m;
Ib_bar = input_param.ib_bar;
g_bar = input_param.g_bar;
U_bar = input_param.u_bar;
F0_bar = input_param.f0_bar;
W = input_param.w;

Lz = zex - z0;
Lzi = zexi - z0;

if Lzi > Lz
    fprintf('\nLzi must be less than Lz!\n');
    pause
end

dz = Lz/(Nz - 1);
dzi = Lzi/(Nzi - 1);
dtau = TauEnd/(Ntau - 1);

% for i=1:Nz
%     ZAxis(i) = (i-1) * dz + z0;
% end
ZAxis = (z0:dz:zex)';

% for i=1:Nzi
%     ZAxisi(i) = (i-1) * dzi + z0;
% end
ZAxisi = (z0:dzi:zexi)';

% for i=1:Nt
%     TAxis(i) = (i-1) * dt;
% end
TauAxis = (0:dtau:TauEnd)';

% BEGIN Вычисление стационарнцх величин   
    gamma_bar = 1 + U_bar/511;
    Beta_bar = sqrt(1 - 1/(gamma_bar*gamma_bar));
    Beta2_bar = 1 - 1/(gamma_bar*gamma_bar);
    Beta_par_bar = Beta_bar/sqrt(g_bar*g_bar + 1);
    Beta_par2_bar = Beta2_bar/(g_bar*g_bar + 1);
    Beta_perp2_bar = Beta2_bar - Beta_par2_bar;
    Beta_perp_bar = sqrt(Beta_perp2_bar);
% END Вычисление стационарных величин

% BIGIN Вычисление нестационарных величин в точках времени без "крышки"
    U_tau = U_bar*ones(Ntau,1); % temporarly
    % U_tau = U_bar*tanh(0.01*TauAxis + 0.1);    

    g_tau = g_bar*sqrt(U_tau./((1 + g_bar^2)*U_bar - g_bar^2*U_tau));

    gamma_tau = 1 + U_tau/511;
    Beta_tau = sqrt(1 - 1./(gamma_tau.*gamma_tau));
    Beta_par_tau = Beta_tau./sqrt(g_tau.*g_tau + 1);
    
    apar_tau = Beta_par_tau/Beta_par_bar;
% END Вычисление нестационарных величин в точках времени без "крышки"

% BEGIN Новое время        
    TAxis_aux = (4*apar_tau)./g_tau.^2.*TauAxis; % точки времени с "крышкой", соответствующие равномерной шкале "безкрышечного"
   
    St = griddedInterpolant(TAxis_aux, TauAxis);

    dt = dzi;
    TAxis = TAxis_aux(1):dt:TAxis_aux(Ntau); % новая шкала времени (dt = dzi)

    TauAxis_aux = St(TAxis); % точки "безкрышечного" времени, соответствующие времени с "крышкой"   
    TauAxis_aux = TauAxis_aux';  

    Nt = length(TAxis); % число точек времени в которых будет найдено решение
% END Новое время

% BEGIN Вычисление НЕстационарнцх величин в точках времени с "крышкой"
    U = U_bar*ones(Nt,1); % здесь нужно вычислить U в точках шкалы TauAxis_aux
    % U = U_bar*tanh(0.01*TauAxis_aux + 0.1);

    g = g_bar*sqrt(U./((1 + g_bar^2)*U_bar - g_bar^2*U));

    gamma = 1 + U/511;
    Beta = sqrt(1 - 1./(gamma.*gamma));
    Beta2 = 1 - 1./(gamma.*gamma);
    Beta_par = Beta./sqrt(g.*g + 1);
    Beta_par2 = Beta2./(g.*g + 1);
    Beta_perp2 = Beta2 - Beta_par2;
    Beta_perp = sqrt(Beta_perp2);

    Ib = Ib_bar*exp(B*(sqrt(U)-sqrt(U_bar)));

    F = F0_bar;

    ag = gamma/gamma_bar;
    ai = Ib/Ib_bar;
    apar = Beta_par/Beta_par_bar;
    aperp = Beta_perp/Beta_perp_bar;
    kappa = g.^2./(4*apar);
% END Вычисление НЕстационарнцх величин в точках времени с "крышкой"

% BEGIN Вычисление kpar2
    kpar2 = zeros(Nz,1);
    
    fileID = fopen("rr.dat");
    d = fscanf(fileID, '%f %f', [2 Inf]);
    fclose(fileID);
    d = d';
    
    Nz3 = Nz - nexl - nexr;
    
    S1 = griddedInterpolant(d(:,1), d(:,2),'linear');
    x = linspace(d(1,1), d(end,1), Nz3)';
    kpar2(nexl + 1:Nz - nexr) = S1(x);
    
    kpar2(1:nexl) = kpar2(nexl + 1);
    kpar2(Nz - nexr + 1:end) = kpar2(Nz - nexr);
    
    fileID = fopen('delta.dat','w');
    for i=1:Nz
        fprintf(fileID,'%17.8e\t%17.8e\n', (i-1)*dz, kpar2(i));
    end
    fclose(fileID);
    
    % kpar2 = kpar2 + Delta - 1i*e;
    kpar2(:) = -Delta/kappa(1);
    Delta = 0;
% END Вычисление kpar2 

InitialField = zeros(Nz,1);
ZBEG = z0;
ZEND = zex;
IND1 = (ZAxis > ZBEG & ZAxis < ZEND);
InitialField(IND1,1) = a0*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
% InitialField(IND1,1) = sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;
% InitialField = ones(length(ZAxis),1) + 1i*ones(length(ZAxis),1);
% InitialField(IND1,1) = a0*cos(pi * (ZAxis(IND1)) / (ZEND - ZBEG)).^2;

if INTT < 1
    error('Too small Tend');
end

if INTZ < 1
    error('Too small Zend');
end

if INTT > 1 && INTZ > 1
    OUTNt = fix((Nt-1)/INTT) + 1;
    OUTNz = fix((Nz-1)/INTZ) + 1;
elseif INTT == 1 && INTZ > 1
    OUTNt = Nt;
    OUTNz = fix((Nz-1)/INTZ) + 1;
elseif INTT > 1 && INTZ == 1
    OUTNt = fix((Nt-1)/INTT) + 1;
    OUTNz = Nz;
else
    OUTNt = Nt;
    OUTNz = Nz;
end

%OUTJ = zeros(OUTNz,OUTNt);
%OUTF = zeros(OUTNz,OUTNt);
OUTZAxis = zeros(OUTNz,1);
OUTTAxis = zeros(OUTNt,1);

OUTZAxis(1) = z0;
for i = 2:OUTNz
    OUTZAxis(i) = (i-1)*INTZ*dz + z0;
end

OUTTAxis(1) = 0;
for i = 2:OUTNt
    OUTTAxis(i) = (i-1)*INTT*dt;
end

fileID = fopen('input_m.txt','w');
% fprintf(fileID,'Nz1 = %i\n', int64(Nz1));
fprintf(fileID,'Nzi = %i\n', int64(Nzi));
fprintf(fileID,'Nt = %i\n', int64(Nt));
fprintf(fileID,'Ne = %i\n', int64(Ne));
fprintf(fileID,'Lz = %f\n', Lz);
fprintf(fileID,'Lzi = %f\n', Lzi);
fprintf(fileID,'Tend = %f\n', TauEnd);
% fprintf(fileID,'ZetaEx = %f\n', ZetaEx);
fprintf(fileID,'Delta = %f\n', Delta);
fprintf(fileID,'Ib = %f\n', Ib_bar);
fprintf(fileID,'dz = %f\n', dz);
fprintf(fileID,'dt = %f\n', dt);
fprintf(fileID,'tol = %g\n', tol);
fprintf(fileID,'INTT = %i\n', int64(INTT));
fprintf(fileID,'INTZ = %i\n', int64(INTZ));
fprintf(fileID,'a0 = %f\n', a0);
fclose(fileID);

timerVal = tic;
[OUTF, OUTJ, p, Eff, Omega, ConLow, jout] = gyroscr(Nz, Nzi, Nt, Ne, ZAxis, ZAxisi, TAxis, TauAxis_aux, Delta, ...
    Ib, ag, ai, apar, aperp, kappa, m, F, W, dt, dz, dzi, tol, kpar2, INTT, INTZ, OUTNz, OUTNt, InitialField);
ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

Folder = 'results/';

hash = datestr(now,30);
FolderName=sprintf('%s%s', Folder, hash);
mkdir(FolderName);

% Вывод в .mat файл
%     fileResults = sprintf('%s/%s', FolderName, 'results.mat');
%     save(fileResults,"RES","-v7.3");

% Вывод в .dat файл
OUTFre = real(OUTF);
OUTFim = imag(OUTF);
OUTJre = real(OUTJ);
OUTJim = imag(OUTJ);
pre = real(p);
pim = imag(p);

% fileID = fopen('fre_m.dat','w');
% for idx0 = 1:OUTNz
%     fprintf(fileID, "%e\t", OUTZAxis(idx0));
%     for idx1 = 1:OUTNt
%         fprintf(fileID, "%e\t%f\t", OUTFre(idx0, idx1));
%     end
%     fprintf(fileID, "\n");
% end
% fclose(fileID);
% 
% fileID = fopen('fim_m.dat','w');
% for idx0 = 1:OUTNz
%     fprintf(fileID, "%e\t", OUTZAxis(idx0));
%     for idx1 = 1:OUTNt
%         fprintf(fileID, "%e\t%f\t", OUTFim(idx0, idx1));
%     end
%     fprintf(fileID, "\n");
% end
% fclose(fileID);
% 
% fileID = fopen('ire_m.dat','w');
% for idx0 = 1:OUTNz
%     fprintf(fileID, "%e\t", OUTZAxis(idx0));
%     for idx1 = 1:OUTNt
%         fprintf(fileID, "%e\t%f\t", OUTJre(idx0, idx1));
%     end
%     fprintf(fileID, "\n");
% end
% fclose(fileID);
% 
% fileID = fopen('iim_m.dat','w');
% for idx0 = 1:OUTNz
%     fprintf(fileID, "%e\t", OUTZAxis(idx0));
%     for idx1 = 1:OUTNt
%         fprintf(fileID, "%e\t%f\t", OUTJim(idx0, idx1));
%     end
%     fprintf(fileID, "\n");
% end
% fclose(fileID);
% 
% fileID = fopen('e_m.dat','w');
% for idx1 = 1:Nt
%     fprintf(fileID, "%e\t%e\n", TAxis(idx1), Eff(idx1));
% end
% fclose(fileID);
% 
% fileID = fopen('w_m.dat','w');
% for idx1 = 1:Nt
%     fprintf(fileID, "%e\t%e\n", TAxis(idx1), Omega(idx1));
% end
% fclose(fileID);

fileParameters = sprintf('%s/%s', FolderName, 'parameters.txt');
fileResultsFre = sprintf('%s/%s', FolderName, 'fre.dat');
fileResultsFim = sprintf('%s/%s', FolderName, 'fim.dat');
fileResultsJre = sprintf('%s/%s', FolderName, 'jre.dat');
fileResultsJim = sprintf('%s/%s', FolderName, 'jim.dat');
fileResultsEff = sprintf('%s/%s', FolderName, 'e.dat');
fileResultsW = sprintf('%s/%s', FolderName, 'w.dat');
fileResultsConLow = sprintf('%s/%s', FolderName, 'cl.dat');
fileT = sprintf('%s/%s', FolderName, 'Time.dat');
fileZ = sprintf('%s/%s', FolderName, 'Z.dat');
fileResultsPre = sprintf('%s/%s', FolderName, 'pre.dat');
fileResultsPim = sprintf('%s/%s', FolderName, 'pim.dat');

fileID = fopen(fileParameters ,'w');
fprintf(fileID,'Nz = %f\n', Nz);
fprintf(fileID,'Nt = %f\n', Nt);
fprintf(fileID,'Ne = %f\n', Ne);
% fprintf(fileID,'ZetaEx = %f\n', ZetaEx);
% fprintf(fileID,'TauEnd = %f\n', TauEnd);
fprintf(fileID,'Delta = %f\n', Delta);
fprintf(fileID,'Ib = %f\n', Ib_bar);
fprintf(fileID,'dz = %f\n', dz);
fprintf(fileID,'dt = %f\n', dt);
fprintf(fileID,'tol = %g\n', tol);
fprintf(fileID,'Last tau index = %g\n', jout);
fprintf(fileID,'g = %f\n', g_bar);
fclose(fileID);
 
OUTFre = [OUTZAxis OUTFre];
OUTFim = [OUTZAxis OUTFim]; 
OUTJre = [OUTZAxis OUTJre];
OUTJim = [OUTZAxis OUTJim];
pre = [ZAxisi pre];
pim = [ZAxisi pim];

ConLow = [TauAxis_aux ConLow'];
OUTEff = [TauAxis_aux Eff'];
OUTOmega = [TauAxis_aux Omega'];

save(fileResultsFre, 'OUTFre', '-ascii');
save(fileResultsFim, 'OUTFim', '-ascii');
save(fileResultsJre, 'OUTJre', '-ascii');
save(fileResultsJim, 'OUTJim', '-ascii');
save(fileResultsEff, 'OUTEff', '-ascii');
save(fileResultsW, 'OUTOmega', '-ascii');
save(fileResultsConLow, 'ConLow', '-ascii');
save(fileResultsPre, 'pre', '-ascii');
save(fileResultsPim, 'pim', '-ascii');

save(fileT, 'TauAxis_aux', '-ascii');
save(fileZ, 'ZAxis', '-ascii');
end
