% Xfeature_sandbox - test ground for Xfeature functionality derived from
% ITD-2009 grapes.

% ======NEW deristats==========
[D, DS1, E1] = readTKABF('RG09110', '3-1-',8); % 625 Hz
[CH, RW, AP, EPSP]=deristats(D, 0.07, -8);
[APw, T] = getsnips({RW.sam RW.dt}, [AP.tEPSPmaxRate], [-0.7 1.7]);
%hist(denan([AP.tpeak]-[AP.tEPSPmaxRate]), 31);
Iap = categorize([AP.tpeak]-[AP.tEPSPmaxRate], 0.1:0.1:0.6);
catplot(T, APw, Iap)
[EPSPw, T] = getsnips({RW.sam RW.dt}, [EPSP.tmaxRate], [-0.7 1.7]);
Iep = categorize([EPSP.tpeak]-[AP.tEPSPmaxRate], 0.1:0.1:0.6);
XcycleHist(RW, AP, EPSP, [10 50]);


[D, DS1, E1] = readTKABF('RG09110', '3-1-',8); % 625 Hz
[CH, RW]=deristats(D, 0.07, -8, [], 3);
APd = CH([CH.isAPdownFlank]); [APdw, T] = getsnips({RW.sam RW.dt}, [APd.tmaxrate], [0 1]);
CatPlot(T, APdw, Iap, 'residue')
EPSP = CH([CH.isEPSP]); [EPSPw, T] = getsnips({RW.sam RW.dt}, [EPSP.tmaxrate], [-0.5 0.8]);
Iep = categorize([EPSP.maxRate], 10);

X=xspikes(D);
%===========NICE cell===========
%   - CF ~ 1400 Hz
%   - tau = 3.9 ms (BN 1 2:3 4:5)
%   - fairly monotonic rate fnc 
%   - large prepotential

% 3-1-FS   237::5000 Hz  60 dB       no mod    23 x 10 x 50/400 ms   L
[D, DS1, E1] = readTKABF('RG09110', '3-1-',1); % 237 Hz
X1(1) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.802],0, 0.5);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',2); % 272 Hz
X1(2) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',3); % 312 Hz
X1(3) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',4); % 312 Hz
X1(4) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',5); % 412 Hz
X1(5) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',6); % 474 Hz
X1(6) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',7); % 544 Hz
X1(7) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',8); % 625 Hz
X1(8) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',9); % 718 Hz
X1(9) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',10); % 825 Hz
X1(10) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);
[D, DS1, E1] = readTKABF('RG09110', '3-1-',11); % 947 Hz
X1(11) = xspikes(D,1:11,1,[0.4 0.7 1.5 2.4 2.6 2.8], 0, 0.4);

M = spikemetrics(X1(11),0,9, -[0.45 0.2], 1.4, 0.5, 0.4); % note need for explicot APthr

%===========RG09106 pen10unit13 (SGSR cell 14)===========
%   - CF 2020 Hz
%   - nonmonotonic rate fnc 
%   - no clear prespike, but clear EPSP and large AP amp

%-104 14-22-SPL 2020 Hz 0:5:70 dB no mod 15 x 10 x 50/400 ms L 
[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',2); % 5 dB (Poster)
Y(1) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M5dB = SpikeMetrics(Y(1),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',5); % 20 dB
Y(2) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M20dB = SpikeMetrics(Y(2),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',6); % 25 dB (Poster)
Y(3) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M25dB = SpikeMetrics(Y(3),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',10); % 45 dB
Y(4) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M45dB = SpikeMetrics(Y(4),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',11); % 50 dB
Y(5) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M50dB = SpikeMetrics(Y(5),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',12); % 55 dB (Poster)
Y(6) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M55dB = SpikeMetrics(Y(6),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',13); % 60 dB
Y(7) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M60dB = SpikeMetrics(Y(7),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

[D, DS1, E1] = readTKABF('RG09106', '14-22-SPL',14); % 65 dB
Y(8) = xspikes(D,1:11, 1, [1 2 3.8 5 6 7], 0, 1.2);
M65dB = SpikeMetrics(Y(8),'>0',9, -[0.45 0.2], 4, 0.5, 0.4);

save Xspikes_storage X1 Y M*

load Xspikes_storage
% templates based on EPSP-AP interval
[TM, T, NN, Ie] = Mtemplate(M5dB,'itvEPSP_AP', [0 0.25 0.27 0.31 0.35 0.4 0.6],-2, '\Deltat(EPSP-AP) ', 'ms');
[TM, T, NN, Ie] = Mtemplate(M55dB,'itvEPSP_AP', [0 0.25 0.27 0.31 0.35 0.4 0.6],-2, '\Deltat(EPSP-AP) ', 'ms');
% SPL has an effect on amplitude of APs, not so much the temporal
% characteristics

% hist of inter-AP intervals
subplot(2,2,1);hist(M5dB.itvAPPrevAP,linspace(1,50,30)); title('5 dB'); xlim([0 50]); 
subplot(2,2,2);hist(M25dB.itvAPPrevAP,linspace(1,50,30)); title('25 dB'); xlim([0 50]); 
subplot(2,2,3);hist(M45dB.itvAPPrevAP,linspace(1,50,30)); title('45 dB'); xlim([0 50]); 
subplot(2,2,4);hist(M55dB.itvAPPrevAP,linspace(1,50,30)); title('55 dB'); xlim([0 50]);

% effect of AP-event interval on EPSP size
MrkFail=struct('color', [0 0 0.7], 'marker', 'o', 'markerfacecolor', [0 0 0.7], 'linestyle', 'none');
MrkFail(2)=struct('color', [0.9 0 0], 'marker', 'o', 'markerfacecolor', [0.7 0 0], 'linestyle', 'none');
f1; qscatter(M5dB, 'itvPrevAP', 'MaxEPSPrate', 'hasAP', [], MrkFail); xlog125([0.5 50]);
f2; qscatter(M55dB, 'itvPrevAP', 'MaxEPSPrate', 'hasAP', [], MrkFail); xlog125([0.5 50]);

% template based on EPSP-AP interval
[TM, T, NN] = Mtemplate(M,'itvEPSP_AP', [0 0.24 0.25 0.27 0.31 0.4 0.6],1);
[TM, T, NN] = Mtemplate(M,'MaxEPSPrate', [15 17 20 25 30 inf],1);

% downward slope of APs & EPSPs have non-overlapping distributions
f1; qscatter(M25dB, 'itvPrevAP', 'MinEPSPrate'); xlog125([0.5 50])
f2; qscatter(M25dB, 'itvPrevAP', 'MinAPrate', 'hasAP'); xlog125([0.5 50])

% no clear relation rising and fallin EPSP slope
qscatter(M25dB, 'MaxEPSPrate', 'MinEPSPrate');
% tight connection AP size and AP down rate 
qscatter(M25dB, 'MaxAP', 'MinAPrate');

% align APs
[Ev, T, TM] = getevent(X1(1),find(X1(1).catIdx>0)); plot(T, Ev)
tmshift=[0 -0.03 -0.075 -0.12 -0.15 -0.15]; Scal = [1 1.1 1.2 1.3 1 1];
clf; for icat=1:6, xplot(T+tmshift(icat), Scal(icat)*TM(:,icat), lico(icat),'linewidth',2); end; grid on
% lin comb of TM & time-integral to reconstruct V(intra)
ir=60:193; 
TMi = 2*cumsum(TM(ir,:)-0.12)*X1(1).dt+TM(ir,:);
clf; for icat=1:6, xplot(T(ir)+tmshift(icat), TMi(:,icat), lico(icat),'linewidth',2); end; grid on
ir=70:193; 
TMi = 3*cumsum(TM(ir,:)-0.15)*X1(1).dt+TM(ir,:);
clf; for icat=1:6, xplot(T(ir)+tmshift(icat), TMi(:,icat), lico(icat),'linewidth',2); end; grid on

% principal component analysis
[Ev, Te, TM, Ienf] = getevent(X1(2),-(1:4)); % non-failures
[Coeff, Score] = princomp(Ev);
Nc = 3; REv = Score(:,1:Nc)*Coeff(:,1:Nc)'; % reconstructed events
for iev=1:1e3, ev=Ev(:,iev); rev = REv(:, iev); plot(Te, ev-mean(ev)); xplot(Te, rev,'r'); ff; pause; end


% scatter plots
qscatter(M,'MaxPrespike', 'MaxEPSP',  'hasAP');

% ==========AP size is determined by mix of itvPrevAP & EPSP size
M = SpikeMetrics(X1(6),0,9, -[0.45 0.2], [], 0.5, 0.4);
qscatter(M, 'MaxEPSP', 'MaxAP', 'itvPrevAP',[0.5 1 2 5 10]);
qscatter(M, 'MaxEPSP', 'MaxAP', 'itvPrevEvent',[0.5 1 2 5 10]); % itv to prev event, no AP needed
% subdivision according to itvPrevXXX does not help to make Prespikes
% correlate to anything (further checks needed)
qscatter(M, 'MaxPrespike', 'MaxAP', 'itvPrevEvent',[0.5 1 2 3.5 5 10]);
% delay reduction of PCA??

% deduce "EPSP threshold" that best predicts AP occurence. Evaluate how it
% varies with stimulus parameters.

% ======OLD deristats==========
[D, DS1, E1] = readTKABF('RG09110', '3-1-',8); % 625 Hz 60 dB
S = deristats(D, 0.2,  [0.068  1.203; 0.654  3.168], [0.019 0.109; 0.300 0.442], [0.084 2.038; 1.149 2.066]);
tap = S.Tfall(:,S.isAP); % falling times of APs
[Sap, T]=getsnips(D, mean(tap), [1 1]); % AP snippets
Iap = categorize(S.MagnFall(S.isAP), 6); % categorize according to fall size

[D, DS1, E1] = readTKABF('RG09110', '3-6-',8); % 625 Hz 20 dB
S = deristats(D, 0.2);
tap = S.Tfall(:,S.isAP); % falling times of APs
[Sap, T]=getsnips(D, mean(tap), [1 1]); % AP snippets
Iap = categorize(S.MagnFall(S.isAP), 6); % categorize according to fall size

tep = S.Tfall(:,S.isEPSP); 
[Sep, T]=getsnips(D, mean(tep), [1 1]); % EPSP snippets
Iep = categorize(S.MagnFall(S.isEPSP), 5);

f1; CatPlot(T, Sap, Iap);
f2; CatPlot(T, Sep, Iep);




% 3-2-FS   237::5000 Hz  50 dB       no mod    23 x 10 x 50/400 ms   L
% 3-3-FS   237::5000 Hz  40 dB       no mod    23 x 10 x 50/400 ms   L
% 3-4-FS   650 Hz        50 dB       no mod    1 x 100 x 400/1000 ms L
% 3-5-FS   237::5000 Hz  30 dB       no mod    23 x 10 x 50/400 ms   L
% 3-6-FS   237::5000 Hz  20 dB       no mod    23 x 10 x 50/400 ms   L
% 3-7-FS   237::5000 Hz  10 dB       no mod    23 x 10 x 50/400 ms   L
% 3-8-FS   237::5000 Hz  0 dB        no mod    23 x 10 x 50/400 ms   L
% 3-9-FS   237::5000 Hz  70 dB       no mod    23 x 10 x 50/400 ms   L
% 3-10-FS  2900 Hz       70 dB       no mod    1 x 10 x 50/400 ms    L
% 3-11-FS  2900 Hz       70 dB       no mod    1 x 10 x 400/400 ms   L
% 3-12-FS  2900 Hz       70 dB       no mod    1 x 10 x 400/400 ms   L
% 3-13-FS  650 Hz        50 dB       no mod    1 x 50 x 400/1000 ms  L
% 3-14-FS  650 Hz        50 dB       no mod    2 x 50 x 400/1000 ms  L
% 3-15-SPL 1250 Hz       0:5:70 dB   no mod    15 x 10 x 50/400 ms   L
% 3-16-SPL 650 Hz        0:5:70 dB   no mod    15 x 10 x 50/400 ms   L
% 3-33-SPL 1450 Hz       0:5:70 dB   no mod    15 x 10 x 50/400 ms   L
% 3-34-SPL 1450 Hz       -25:5:40 dB no mod    14 x 10 x 50/400 ms   L



[D2, DS2, E2] = readTKABF('RG09110', '3-2-',1);
X2 =xspikes(D2,1:11,1,[0.4 1.5 2.5 2.75]);

[D3, DS3, E3] = readTKABF('RG09110', '3-3-',1);
X3 =xspikes(D3,1:11,1,[0.4 1 2.1 2.3],1);

[D5, DS5, E5] = readTKABF('RG09110', '3-5-',1);
X5 = xspikes(D5,1:11,1,[0.4 1 2.1 2.3],1);













