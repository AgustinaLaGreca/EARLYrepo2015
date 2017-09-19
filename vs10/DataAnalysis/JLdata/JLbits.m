% JLbits - bits & pieces from JLbeat analysis

% autocorrelation (todo: truncate spikes)
dt = D.AD(1).dt_ms;
% 40 dB contra
AC=[]; for icond = 1:13; [D DS]=readTKabF('RG10219','1-6-',icond); Y = D.AD(1).samples; Y = Y-mean(Y); AC = [AC, xcorr(Y,1000)]; icond, end; dplot(dt, AC);
% 40 dB ipsi
AC=[]; for icond = 1:13; [D DS]=readTKabF('RG10219','1-8-',icond); Y = D.AD(1).samples; Y = Y-mean(Y); AC = [AC, xcorr(Y,1000)]; icond, end; dplot(dt, AC);
% crosscorr
ExpID = 'RG10219';
Rec_ipsi = '1-4-'; Rec_contra = '1-5-';
Rec_ipsi = '1-8-'; Rec_contra = '1-6-';
NsamXC = 1000;
XC=[]; 
for icond = 1:13, 
    [D DS]=readTKabF(ExpID,Rec_contra,icond); Yc = D.AD(1).samples; Yc = Yc-mean(Yc); 
    dt = D.AD(1).dt_ms;
    freq_c = DS.fcar(icond, 2);
    [D DS]=readTKabF(ExpID,Rec_ipsi,icond); Yi = D.AD(1).samples; Yi = Yi-mean(Yi); 
    freq_i = DS.fcar(icond, 1);
    % time warp: stretch contra wave to equalize its stim freq w ipsi stim
    Nsam = numel(Yc);
    tt_c = dt*(0:Nsam-1)'; % true time axis of contra (& ipsi)
    Nsamnew = round(Nsam*freq_c/freq_i);
    tt_cnew = linspace(0,max(tt_c), Nsamnew)'; % denser time axis
    Yc = interp1(tt_c, Yc, tt_cnew);
    XC = [XC, xcorr(Yi, Yc, NsamXC)];
end; 
tau = dt*(-NsamXC:NsamXC)';
plot(tau, XC)
% Xcorr, but now for cyclic averages
XC=[]; 
for icond = 1:13, 
    Ji = JbB_50dB(icond);
    Yi = Ji.Y1;
    Yc = Ji.Y2;
    dt = Ji.dt;
    freq_i = Ji.Freq1;
    freq_c = Ji.Freq2;
    % time warp: stretch contra wave to equalize its stim freq w ipsi stim
    Nsam = numel(Yc);
    tt_c = dt*(0:Nsam-1)'; % true time axis of contra (& ipsi)
    Nsamnew = round(Nsam*freq_c/freq_i);
    tt_cnew = linspace(0,max(tt_c), Nsamnew)'; % denser time axis
    Yc = interp1(tt_c, Yc, tt_cnew);
    XC = [XC, xcorr(Yi, Yc, NsamXC)];
end; 
tau = dt*(-NsamXC:NsamXC)';
plot(tau, XC)

