function JL_wcartoon(Uidx_B, DT, T0);
% JL_wcartoon - cartoon comparing monaural & binaural responses
%    JL_wcartoon(Uidx_B, DT, T0);
%      DT is wnbdow length
%      T0 is array of start times [B I C]

global DBc
if isempty(DBc),
    DBc = structMerge('UniqueRecordingIndex', JLdbase, JL_rec_select);
    DBc = DBc([DBc.CanBeUsed]);
end

% find monaural partners
D = JLdatastruct(Uidx_B);
Dsim = structmatch(DBc, D, {'icell_run' 'SPL' 'freq'});
Uidx_I = Dsim([Dsim.chan]=='I').UniqueRecordingIndex
Uidx_C = Dsim([Dsim.chan]=='C').UniqueRecordingIndex

[W, Wb] = JLgetRec(Uidx_B);
[dum, Wi] = JLgetRec(Uidx_I);
[dum, Wc] = JLgetRec(Uidx_C);

dt = W.dt;
Nsm = 10;
Wb = smoothen(Wb, Nsm);
Wi = smoothen(Wi, Nsm);
Wc = smoothen(Wc, Nsm);


Ti = 1e3/D.freq; % ms ipsi cycle
Tc = 1e3/(D.freq+4); % ms contra cycle

T0_B = T0(1)
tshift_i = mod(T0_B, Ti)
tshift_c = mod(T0_B, Tc)
T0_I = Ti*floor(T0(2)/Ti)+tshift_i
T0_C = Tc*floor(T0(3)/Tc)+tshift_c

Twin = [0 DT]
[dum, dum, Wb] = cutFromWaveform(dt, Wb, T0_B, Twin);
[dum, dum, Wi] = cutFromWaveform(dt, Wi, T0_I, Twin);
[dum, dum, Wc] = cutFromWaveform(dt, Wc, T0_C, Twin);

Wb = Wb-mean(Wb);
Wi = Wi-mean(Wb);
Wc = Wc-mean(Wb);
Yshift = 3*std(Wb);
dplot(dt, Wb, 'k', 'linewidth', 2);
xdplot(dt, Wc+Yshift, 'color', [0.8 0 0 ], 'linewidth', 2);
xdplot(dt, Wi+2*Yshift, 'color', [0 0 0.8], 'linewidth', 2);






