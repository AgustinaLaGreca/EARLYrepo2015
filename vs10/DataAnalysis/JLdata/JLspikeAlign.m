function T = JLspikeAlign(Rec, DT, PreWin);
% JLspikeAlign - align APs of JL data to their max EPSP slopes
%   JLspikeAlign(Rec, DT, PreWin);
%     Rec: recording ID (see JLdatastruct)
%      DT: smoothing (ms) applied to dV/dt to find EPSP slopes [0.1 ms]
%  PreWin: pre-AP-peak window (ms) for steepest EPSPs slope to occur [-0.5 0.1] ms] 
% 

[DT, PreWin] = arginDefaults('DT/PreWin', 0.1, [-0.5 -0.1]); 

% retrieve JLbeat output
S=JLdatastruct(Rec);
%Jb=JLgetBeatStruct(S);
%dt = Jb.dt;

[SPTraw, Snip, APslope] = deal([]);
for ii=1:numel(S),
    s = S(ii);
    [D, DS, L]=readTKABF(s);
    rec = D.AD(1).samples; % entire recording, including pre- and post-stim parts
    NsamTot = D.NsamPerSweep/D.Nchan;
    cutWin = [s.APwindow_start s.APwindow_end];
    %
    dt = D.dt_ms;
    [rec, SPTraw_ii, Tsnip, Snip_ii, APslope_ii] = ...
        APtruncate2(rec, [s.APthrSlope, s.CutoutThrSlope], dt, cutWin, 1); % last 1: replace APs by nans
    SPTraw = [SPTraw; SPTraw_ii];
    Snip = [Snip, Snip_ii];
    APslope = [APslope, APslope_ii];
end
% ispike = 11:60;
% Snip = Snip(:,ispike);

T = local_Timing(Snip, Tsnip, DT, dt, PreWin, APslope);

set(figure,'units', 'normalized', 'position', [0.0266 0.173 0.852 0.767])
subplot(2,1,1);
plot(T.Tsnip, T.Snip);
ylabel('V_{rec} (mV)');
subplot(2, 1,2);
plot(T.TsnipE, T.Snip);
xlabel('Time (ms)');
ylabel('V_{rec} (mV)');

set(figure,'units', 'normalized', 'position', [0.598 0.236 0.355 0.659])
subplot(3,1,1);
plot(T.maxSlope , -T.TmaxSlope, '.')
ylim([0 inf])
xlabel('max EPSP'' (dV/dt)');
ylabel('EPSP-AP latency (ms)');
subplot(3,1,2);
plot(T.maxSlope , T.APpeak-T.APdip, '.')
xlabel('max EPSP'' (dV/dt)');
ylabel('APpeak-APdip (mV)')
%ylim([0 inf])
subplot(3,1,3);
plot(-T.TmaxSlope , T.APpeak-T.APdip, '.')
xlim([0 inf])
xlabel('EPSP-AP latency (ms)');
ylabel('APpeak-APdip (mV)')

%====================================================
%====================================================

function T = local_Timing(Snip, Tsnip, DT, dt, PreWin, APslope);
dSnip = smoothen(diff(Snip), DT, dt)/dt;
TdSnip = Tsnip(1:end-1)+dt/2;
qPre = betwixt(TdSnip, PreWin);
dPreAP = dSnip(qPre, :);
TdPreAP = TdSnip(qPre);
[maxSlope, imaxSlope] = max(dPreAP); % peak dV/dt in pre-AP window
TmaxSlope = TdPreAP(imaxSlope); TmaxSlope = TmaxSlope(:).';
TsnipE = SameSize(TmaxSlope, Tsnip);
TsnipE = SameSize(Tsnip, TsnipE) - TsnipE;
%dsize(Tsnip, Snip, dSnip, dPreAP, TmaxSlope, TsnipE)
APpeak = max(Snip);
APdip = min(Snip);
T = CollectInStruct(DT, dt, PreWin, Tsnip, Snip, '-', ...
    TsnipE, dSnip, TdSnip, dPreAP, TdPreAP, '-', ...
    APslope, TmaxSlope, maxSlope, APpeak, APdip);





