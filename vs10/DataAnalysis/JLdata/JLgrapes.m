function JLgrapes;
% JLgrapes selected MSO analyses

% function JLbeat(ExpID, RecID, icond, Nxmean, Vclip);

function RG09178
% RG10178 >>>>cell 1 ("EPSPs plus Aps");  Weakly binaural. Spikes only at 100 Hz. 
%      70 dB 2xbin
% many FS at different levels (!)
% 1 x CFS; % NTD, rho=1,-1,0; % BN
JLreadBeats('RG09178', 1);
JLbeatPlot(JbB_70dB(1)); % 100 Hz  R = [0.31 0.5 0.39] 6 spikes; rho=0.55; multiple peaks per car cycle?! Not very binaural
JLbeatPlot(JbB2_70dB(1));% 100 Hz  R = [0 0 0] 1 spikes; rho=0.54; similar EPSPs, but spikes ar gone
% no spikes at higher freqs
JLbeatPlot(JbB_70dB(2)) % 250 Hz again: higher harmonics dominate;
JLbeatPlot(JbB_70dB(3)) % 400 Hz looks better - but why trust this one? ..
% ... spectral is peak always @ ~750 Hz, indep of Fcar!?
JLbeatPlot(JbB_70dB(4)); % 550 Hz
JLbeatPlot(JbB_70dB(5)); % 700 Hz
JLbeatPlot(JbB_70dB(6)); % 850 Hz
JLbeatPlot(JbB_70dB(7)); % 1000 Hz
JLbeatPlot(JbB_70dB(8)); % 1150 Hz  NOISY
%
% RG10178 >>>>cell 3 ("big spikes, not much EPSPs?");  decent binaural cell
%      70 dB bin
%      70 dB, 100 Hz: beat freq varied
%      CFS, FS, NTD, everything @ 70 dB
%  only spikes @ 100 Hz
JLreadBeats('RG09178', 3); 
JLbeatPlot(JbB_70dB(1)); % 100 Hz  R = [0.74 0.92 0.7] 27 spikes; rho=0.64; crappy STA though
JLbeatPlot(JbB_70dB(2))  % 250 Hz  R = [0 0 0] 0 spikes; rho=0.32;
JLbeatPlot(JbB_70dB(3)); % 400 Hz  R = [0 0 0] 0 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(4)); % 550 Hz  R = [0 0 0] 0 spikes; rho=0.1;
JLbeatPlot(JbB_70dB(5)); % 700 Hz  R = [0 0 0] 0 spikes; rho=0.097;
JLbeatPlot(JbB_70dB(6)); % 850 Hz  R = [0 0 0] 0 spikes; rho=0.041;
% FS sustained & onset 
ABFplot({'RG09178', '3-4-FS', 1}, 2:11, 1); % 63 Hz, 70 dB
ABFplot({'RG09178', '3-4-FS', 2}, 2:11, 1);  % 78 Hz
ABFplot({'RG09178', '3-4-FS', 3}, 2:11, 1);  % 96 Hz
ABFplot({'RG09178', '3-4-FS', 4}, 2:11, 1);  % 118 Hz
ABFplot({'RG09178', '3-4-FS', 5}, 2:11, 1);  % 146 Hz
ABFplot({'RG09178', '3-4-FS', 6}, 2:11, 1);  % 179 Hz
%==============100-Hz only======================
APslope = -15; Nav = 1;
JLbeat('RG09178', '3-7-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-8-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-9-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
% 30-s beats
JLbeat('RG09178', '3-10-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-14-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-15-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-20-',1, Nav, APslope); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-17-',1, Nav, APslope); % 10-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-18-',1, Nav, APslope); % 2-Hz beat, 100-Hz carrier

JLbeat('RG09178', '3-16-',1,5,[-15 15]); % 4-Hz beat, 100-Hz carrier
JLbeat('RG09178', '3-19-',1,5,[-15 15]); % 4-Hz beat, 100-Hz carrier

JLbeat('RG09178', '3-17-',1,5,[-15 15])
JLbeat('RG09178', '3-17-',1,5,[-15 15])
JLbeat('RG09178', '3-17-',1,5,[-15 15])
JLbeat('RG09178', '3-17-',1,5,[-15 15])

%========single tones===========================
Rec70  = {'1-5-' '1-3-' '1-4-' } ; % 70 dB
JLtone('RG09178', Rec70, 5,5); % 109 Hz
JLtone('RG09178', Rec70, 6,5); % 134 Hz
JLtone('RG09178', Rec70, 7,5); % 165 Hz
JLtone('RG09178', Rec70, 8,5); % 203 Hz
JLtone('RG09178', Rec70, 9,5); % 250 Hz
JLtone('RG09178', Rec70, 10,5); % 308 Hz
JLtone('RG09178', Rec70, 11,5); % 379 Hz
JLtone('RG09178', Rec70, 12,5); % 467 Hz
JLtone('RG09178', Rec70, 13,5); % 574 Hz
JLtone('RG09178', Rec70, 14,5); % 707 Hz
JLtone('RG09178', Rec70, 15,5); % 871 Hz
%
Rec60  = {'1-12-' '1-13-' '1-11-' } ; % 60 dB
JLtone('RG09178', Rec60, 5,5); % 109 Hz
JLtone('RG09178', Rec60, 6,5); % 134 Hz
JLtone('RG09178', Rec60, 7,5); % 165 Hz
JLtone('RG09178', Rec60, 8,5); % 203 Hz
JLtone('RG09178', Rec60, 9,5); % 250 Hz
JLtone('RG09178', Rec60, 10,5); % 308 Hz
JLtone('RG09178', Rec60, 11,5); % 379 Hz
JLtone('RG09178', Rec60, 12,5); % 467 Hz
JLtone('RG09178', Rec60, 13,5); % 574 Hz
JLtone('RG09178', Rec60, 14,5); % 707 Hz
JLtone('RG09178', Rec60, 15,5); % 871 Hz
Rec50  = {'1-17-' '1-16-' '1-14-' } ; % 60 dB
JLtone('RG09178', Rec50, 5,5); % 109 Hz
JLtone('RG09178', Rec50, 6,5); % 134 Hz
JLtone('RG09178', Rec50, 7,5); % 165 Hz
JLtone('RG09178', Rec50, 8,5); % 203 Hz
JLtone('RG09178', Rec50, 9,5); % 250 Hz
JLtone('RG09178', Rec50, 10,5); % 308 Hz
JLtone('RG09178', Rec50, 11,5); % 379 Hz
JLtone('RG09178', Rec50, 12,5); % 467 Hz
JLtone('RG09178', Rec50, 13,5); % 574 Hz
JLtone('RG09178', Rec50, 14,5); % 707 Hz
JLtone('RG09178', Rec50, 15,5); % 871 Hz

%=======check temporal alignment
Rec109Hz_ipsi = {'1-5-' '1-12' '1-17' '1-21' '1-24' '1-26'};
Rec109Hz_contra = {'1-3-' '1-13' '1-16' '1-20' '1-23' '1-17'}
JLtone('RG09178', Rec109Hz_ipsi, 5, 5); % 109 Hz 20..70 dB SPL, 109 Hz
JLtone('RG09178', Rec109Hz_ipsi, 8, 5); % 109 Hz 20..70 dB SPL, 203 Hz

JLtone('RG09178', Rec109Hz_contra, 5, 5); % 109 Hz 20..70 dB SPL, 109 Hz
JLtone('RG09178', Rec109Hz_contra, 8, 5); % 109 Hz 20..70 dB SPL, 203 Hz

f1; JLtone('RG09178', {'3-24' '3-25' '3-23'}, 4, 15); % 109 Hz 20..70 dB SPL, 109 Hz

%Rec70dBIpsi


% R10198
APslope = -8;
JLbeat('RG10198', '2-15-BFS',1,Nav, APslope);
JLbeat('RG10198', '2-15-BFS',7,Nav, APslope);

function RG10191
% RG10201 cell 1 ("big spikes, no clear EPSPs");  huge, atypical spikes. Only weakly binaural.
%      60 dB 2xbin/ipsi/contra bookkeeping problems
JLreadBeats('RG10191', 1);
% RG10191 cell 1: important to truncate the huge spikes
% == Note that ipsi ear EPSPs have a shallower rising part - this is
% == consistent across carrier freqs. Also, Ipsi is reduced @ 300 Hz, where
% == contra is still thriving.
% == Combination of ipsi & contra is somewhat nonlinear here. At higher
% carrier freqs, it is linearized.
% ---60 dB  many APs, mostly C-locked despite (frequent) equality of C & I waves
JLbeatPlot([JbB_60dB(1) JbB2_60dB(1)]);% 100 Hz  R = [0.21 0.88 0.16] 276 spikes; rho=0.26;C-drv but ~equal C/I wvs
JLbeatPlot([JbB_60dB(2) JbB2_60dB(2)]);% 150 Hz  R = [0.42 0.82 0.32] 459 spikes; rho=0.46;~C driven<->eq C/I wvs
JLbeatPlot([JbB_60dB(3) JbB2_60dB(3)]);% 200 Hz  R = [0.45 0.93 0.37] 533 spikes; rho=0.64;C driven, C wv larger
JLbeatPlot([JbB_60dB(4) JbB2_60dB(4)]);% 250 Hz  R = [0.37 0.9 0.31] 447 spikes; rho=0.57; same
JLbeatPlot([JbB_60dB(5) JbB2_60dB(5)]);% 300 Hz  R = [0.26 0.91 0.18] 377 spikes; rho=0.58;
JLbeatPlot(JbB2_60dB(6));  % 350 Hz  R = [0.23 0.89 0.16] 182 spikes; rho=0.52;
JLbeatPlot(JbB2_60dB(7));  % 400 Hz  R = [0.4 0.85 0.27] 185 spikes; rho=0.46;
JLbeatPlot(JbB2_60dB(8));  % 450 Hz  R = [0.51 0.71 0.21] 172 spikes; rho=0.4;
% higher freqs not available (lost cell) see below
% ===linearity
ifr = 1; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % noisy but linear
ifr = 2; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % linear
ifr = 3; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % linear
ifr = 4; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % linear
ifr = 5; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % linear
ifr = 6; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % small deviations
ifr = 7; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 8; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 9; JLbeatPlot({JbB2_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % clearly lost the cell in bin
% stability. 1-2-BFS vs 1-6 BFS. Shapes are generally well preserved. Ampl has grown a little.
ifr = 2; JLbeatPlot({JbB_60dB(ifr)  JbB2_60dB(ifr)  JbB2_60dB(ifr)}, 3, 'x3');
ifr = 3; JLbeatPlot({JbB_60dB(ifr)  JbB2_60dB(ifr)  JbB2_60dB(ifr)}, 3, 'x3');%shapes preserved. slight Ampl diff
ifr = 4; JLbeatPlot({JbB_60dB(ifr)  JbB2_60dB(ifr)  JbB2_60dB(ifr)}, 3, 'x3'); % same
ifr = 5; JLbeatPlot({JbB_60dB(ifr)  JbB2_60dB(ifr)  JbB2_60dB(ifr)}, 3, 'x3'); % same
% RG10201 >>>>cell 2 ("EPSP like data?")
%      60 dB bin/ipsi/contra
%      CFS
% -9   2-3-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
% -10  2-4-BFS 100*:50*:1000 Hz -60|60 dB 4 Hz beat 19 x 1 x 9000/11000 ms B
% -11  2-5-BFS 100*:50*:1000 Hz 60|-60 dB 4 Hz beat 19 x 1 x 9000/11000 ms B
% XXX -12  2-6-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
JLreadBeats('RG10191', 2);
% ---60 dB  
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0 0 0] 1 spikes; rho=0.34; I mlt peaked
JLbeatPlot(JbB_60dB(2));  % 150 Hz  R = [0.58 0.75 0.34] 9 spikes; rho=0.2;
JLbeatPlot(JbB_60dB(3));  % 200 Hz  R = [0.39 1 0.33] 2 spikes; rho=0.58; I 3-peaked; C 2-peaked
JLbeatPlot(JbB_60dB(4));  % 250 Hz  R = [0.98 0.76 0.62] 3 spikes; rho=0.27; same
JLbeatPlot(JbB_60dB(5));  % 300 Hz  R = [0.93 0.98 0.86] 6 spikes; rho=0.6;
JLbeatPlot(JbB_60dB(6));  % 350 Hz  R = [0.95 0.99 0.91] 3 spikes; rho=0.53;
JLbeatPlot(JbB_60dB(7));  % 400 Hz  R = [0 0 0] 0 spikes; rho=0.63;
JLbeatPlot(JbB_60dB(8));  % 450 Hz  R = [0 0 0] 0 spikes; rho=0.52;
JLbeatPlot(JbB_60dB(9));  % 500 Hz  R = [1 0.99 0.97] 2 spikes; rho=0.53;
JLbeatPlot(JbB_60dB(10));  % 550 Hz  R = [0 0 0] 0 spikes; rho=0.5;
JLbeatPlot(JbB_60dB(11));  % 600 Hz  R = [0 0 0] 1 spikes; rho=0.53;
JLbeatPlot(JbB_60dB(12));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.49;
% RG10191 >>> cell 3 ("1.5 mV spikes, not sure about EPSPs")
%      60 dB bin/ipsi/contra
%      CFS, FS, 
% -14  3-2-BFS 100*:50*:1000 Hz 60 dB     4 Hz beat 19 x 1 x 9000/11000 ms B
% -15  3-3-BFS 100*:50*:700 Hz  60|-60 dB 4 Hz beat 13 x 1 x 9000/11000 ms B
% -16  3-4-BFS 100*:50*:700 Hz  -60|60 dB 4 Hz beat 13 x 1 x 9000/11000 ms B
JLreadBeats('RG10191', 3);
% ---60 dB  % many spikes, mostly driven; wv not terribly predictive
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0.76 0.87 0.58] 449 spikes; rho=0.45; STA: APs occur very late
JLbeatPlot(JbB_60dB(2));  % 150 Hz  R = [0.85 0.61 0.49] 250 spikes; rho=0.28; STA: APs occur very late
JLbeatPlot(JbB_60dB(3));  % 200 Hz  R = [0.9 0.87 0.74] 78 spikes; rho=0.24;
JLbeatPlot(JbB_60dB(4));  % 250 Hz  R = [0.87 0.75 0.55] 28 spikes; rho=0.099;
JLbeatPlot(JbB_60dB(5));  % 300 Hz  R = [0.87 0.8 0.63] 16 spikes; rho=0.14;
JLbeatPlot(JbB_60dB(6));  % 350 Hz  R = [0.81 0.69 0.48] 18 spikes; rho=0.11;
JLbeatPlot(JbB_60dB(7));  % 400 Hz  R = [0.88 0.92 0.7] 4 spikes; rho=0.11;
JLbeatPlot(JbB_60dB(8));  % 450 Hz  R = [0 0 0] 1 spikes; rho=0.084;
% ===linearity
% ----60 dB
ifr = 1; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3');  % linear
ifr = 1; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 6, 'x3'); % linear; EI cell?
% 100 Hz  R = [0.76 0.87 0.58] 449 spikes; rho=0.45;
% 100 Hz  R = [0.88 0.031 0.037] 144 spikes; rho=0.24;
% 100 Hz  R = [0.076 0.88 0.09] 236 spikes; rho=0.21;
ifr = 2; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % pretty linear
ifr = 2; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 6, 'x3');
% 150 Hz  R = [0.85 0.61 0.49] 250 spikes; rho=0.28;
% 150 Hz  R = [0.91 0.061 0.086] 54 spikes; rho=0.28;
% 150 Hz  R = [0.089 0.56 0.056] 182 spikes; rho=0.082;
ifr = 3; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % pretty linear
ifr = 3; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 6, 'x3');
% 200 Hz  R = [0.9 0.87 0.74] 78 spikes; rho=0.24;
% 200 Hz  R = [0.79 0.11 0.089] 36 spikes; rho=0.19;
% 200 Hz  R = [0.15 0.73 0.11] 72 spikes; rho=0.046;
ifr = 4; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 4; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 6, 'x3');
% 250 Hz  R = [0.87 0.75 0.55] 28 spikes; rho=0.099;
% 250 Hz  R = [0.71 0.28 0.15] 24 spikes; rho=0.054;
% 250 Hz  R = [0.13 0.63 0.056] 90 spikes; rho=0.047;
ifr = 5; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 5; JLbeatPlot({JbB_60dB(ifr)  JbI_60dB(ifr)  JbC_60dB(ifr)}, 6, 'x3');
% 300 Hz  R = [0.87 0.8 0.63] 16 spikes; rho=0.14;
% 300 Hz  R = [0.7 0.23 0.31] 26 spikes; rho=0.035;
% 300 Hz  R = [0.25 0.66 0.13] 86 spikes; rho=0.065;




function RG10195
% RG10195 cell 3 ("EPSPs plus Aps or maybe multi unit?")  onset-like reponses
%      65 dB bin/ipsi/contra
%      75 dB bin/ipsi/contra
% FS, FM, BN
%=========FS data clearly show multipeaked character 
JLtone2('RG10195','3-2-F', 4:6, [0 0 1.3], 'r'); xlim([40 60]) % contra
JLtone2('RG10195','3-3-F', 4:6, [0 0.4 1.5], 'r'); xlim([40 60]); % ipsi
% BN: [D,ds,E]=readTKabF('RG10195','3-18-BN', 1, 1000)
% FM: clampfit('RG10195', '3-11-FM',1)
JLreadBeats('RG10195', 3);
% ==65 dB
JLbeatPlot(JbB_65dB(1));  % 100 Hz  R = [0 0 0] 0 spikes; rho=0.47; C very multi-peaked
JLbeatPlot(JbB_65dB(2));  % 200 Hz  R = [0 0 0] 0 spikes; rho=0.62; note high rho. I&C somewhat multi-peaked
JLbeatPlot(JbB_65dB(3));  % 300 Hz  R = [0 0 0] 0 spikes; rho=0.58; ditto
JLbeatPlot(JbB_65dB(4));  % 400 Hz  R = [0 0 0] 0 spikes; rho=0.52;
JLbeatPlot(JbB_65dB(5));  % 500 Hz  R = [0 0 0] 0 spikes; rho=0.54; C dominates
JLbeatPlot(JbB_65dB(6));  % 600 Hz  R = [0 0 0] 0 spikes; rho=0.57; C dominates
JLbeatPlot(JbB_65dB(7));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.52;
JLbeatPlot(JbB_65dB(8));  % 800 Hz  R = [0 0 0] 1 spikes; rho=0.38;
JLbeatPlot(JbB_65dB(9));  % 900 Hz  R = [0 0 0] 0 spikes; rho=0.25;
JLbeatPlot(JbB_65dB(10));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.14;
JLbeatPlot(JbB_65dB(11));  % 1100 Hz  R = [0 0 0] 0 spikes; rho=0.11;
JLbeatPlot(JbB_65dB(12));  % 1200 Hz  R = [0 0 0] 0 spikes; rho=0.084;
JLbeatPlot(JbB_65dB(15));  % 1500 Hz  R = [0 0 0] 0 spikes; rho=0.015;
JLbeatPlot(JbB_65dB(20));  % 2000 Hz  R = [0 0 0] 0 spikes; rho=-0.0062; I contrib gone
% ==75 dB
JLbeatPlot(JbB_75dB(1));  % 100 Hz  R = [0 0 0] 0 spikes; rho=0.36; C is double peaked and dominates
JLbeatPlot(JbB_75dB(2));  % 200 Hz  R = [0 0 0] 0 spikes; rho=0.28; marginal C double-peakedness; C dominates
JLbeatPlot(JbB_75dB(3));  % 300 Hz  R = [0 0 0] 0 spikes; rho=0.37; C&I double peaks; C dominates
JLbeatPlot(JbB_75dB(4));  % 400 Hz  R = [0 0 0] 0 spikes; rho=0.34; peaks merging
JLbeatPlot(JbB_75dB(5));  % 500 Hz  R = [0 0 0] 0 spikes; rho=0.46;
JLbeatPlot(JbB_75dB(6));  % 600 Hz  R = [0 0 0] 0 spikes; rho=0.34;
JLbeatPlot(JbB_75dB(7));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.27; everything turned ~sinusoidal; C dominates
JLbeatPlot(JbB_75dB(10));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.096;
% =====linearity====
% --65 dB 
% to check: is rho consistently higher in Bin re Mon? Also note relatively
% large interval btw Bin & Mon recordings.
ifr = 1; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % ~consistent multipeaks
ifr = 2; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % stable shapes; mutual suppr
ifr = 3; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % ditto
ifr = 4; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % some variation in shape
ifr = 5; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % ditto
ifr = 6; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % I strongly suppressed. C rel. unaffected
ifr = 7; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % I strongly suppressed. C rel. unaffected
ifr = 8; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); % I strongly suppressed. C rel. unaffected
ifr = 9; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3); 
ifr = 10; JLbeatPlot({JbB_65dB(ifr) JbI_65dB(ifr) JbC_65dB(ifr) }, 3);
% --75 dB 
ifr = 1; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % consistent multipeaks. C dominates
ifr = 2; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % ~linear
ifr = 3; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % ~linear
ifr = 4; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % large deviations from linearity
ifr = 5; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % larg deviations; strong I suppr
ifr = 6; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % stable shapes; strong I suppr; mild C suppr
ifr = 7; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % ditto
ifr = 8; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % strong mutual suppr
ifr = 10; JLbeatPlot({JbB_75dB(ifr) JbI_75dB(ifr) JbC_75dB(ifr) }, 3); % same

function RG10196;
% RG10196 >>>> cell 2 ("EPSPs and Aps")  beauty
% Many spikes. Small freq steps -> maybe helps the study of merging peaks.
%      80 dB bin/ipsi/contra [50-Hz steps]
% FS
JLreadBeats('RG10196', 2);
% ===80 dB
% at 50 Hz, spikes are not phase-locked, but waveform is clearly driven (see also the raw data)
% on the other hand, the rho value is low. 
JLbeatPlot(JbB_80dB(1));  % 50 Hz  R = [0.055 0.15 0.15] 53 spikes; rho=0.019; raw data 
% @ 100 Hz, I is unusual (dipped not peaked). C is more regular, double-peaked.
JLbeatPlot(JbB_80dB(2));  % 100 Hz  R = [0.26 0.8 0.15] 67 spikes; rho=0.53;
% 150 Hz: unusual peaks. Locked t oI&C, but no good binaurality
JLbeatPlot(JbB_80dB(3));  % 150 Hz  R = [0.47 0.51 0.067] 212 spikes; rho=0.32;
% 200 Hz: same story: cell phaselocks to I & C but does no coincidence
% detection. "The refractory" dips in the CHISTs, however, show that this
% is not two units.
JLbeatPlot(JbB_80dB(4));  % 200 Hz  R = [0.54 0.48 0.022] 372 spikes; rho=0.53;
% 250-300 Hz; I & C are double peaked, which also shapes the CHISTs. Cell is no
% binaural. B inherits the double peaks from I and C.
JLbeatPlot(JbB_80dB(5));  % 250 Hz  R = [0.3 0.73 0.2] 292 spikes; rho=0.35;
JLbeatPlot(JbB_80dB(6));  % 300 Hz  R = [0.71 0.55 0.35] 157 spikes; rho=0.51;
% 350-450 Hz: cell now does a decent binaural job. Double-peakedness has become less prominent.
JLbeatPlot(JbB_80dB(7));  % 350 Hz  R = [0.77 0.78 0.56] 150 spikes; rho=0.42;
JLbeatPlot(JbB_80dB(8));  % 400 Hz  R = [0.8 0.75 0.6] 155 spikes; rho=0.5;
JLbeatPlot(JbB_80dB(9));  % 450 Hz  R = [0.79 0.8 0.62] 128 spikes; rho=0.49;
% 500+ Hz: mono peaked, good binaurality
JLbeatPlot(JbB_80dB(10));  % 500 Hz  R = [0.78 0.86 0.64] 202 spikes; rho=0.47;
JLbeatPlot(JbB_80dB(11));  % 550 Hz  R = [0.8 0.89 0.69] 295 spikes; rho=0.43; C dominates
JLbeatPlot(JbB_80dB(12));  % 600 Hz  R = [0.82 0.89 0.71] 381 spikes; rho=0.44;
JLbeatPlot(JbB_80dB(13));  % 650 Hz  R = [0.79 0.85 0.68] 345 spikes; rho=0.43;
JLbeatPlot(JbB_80dB(14));  % 700 Hz  R = [0.78 0.85 0.67] 327 spikes; rho=0.4;
JLbeatPlot(JbB_80dB(15));  % 750 Hz  R = [0.8 0.83 0.68] 239 spikes; rho=0.42;
JLbeatPlot(JbB_80dB(16));  % 800 Hz  R = [0.81 0.84 0.7] 163 spikes; rho=0.45;
JLbeatPlot(JbB_80dB(17));  % 850 Hz  R = [0.78 0.79 0.6] 132 spikes; rho=0.46;
% 900+Hz: spike counts drop, but binaurality is still decent
JLbeatPlot(JbB_80dB(18));  % 900 Hz  R = [0.78 0.8 0.66] 82 spikes; rho=0.44;
JLbeatPlot(JbB_80dB(19));  % 950 Hz  R = [0.73 0.8 0.68] 70 spikes; rho=0.4;
JLbeatPlot(JbB_80dB(20));  % 1000 Hz  R = [0.73 0.71 0.58] 60 spikes; rho=0.4;
JLbeatPlot(JbB_80dB(21));  % 1050 Hz  R = [0.71 0.81 0.61] 44 spikes; rho=0.34;
JLbeatPlot(JbB_80dB(22));  % 1100 Hz  R = [0.6 0.75 0.54] 33 spikes; rho=0.32;
JLbeatPlot(JbB_80dB(23));  % 1150 Hz  R = [0.68 0.65 0.52] 29 spikes; rho=0.26;
JLbeatPlot(JbB_80dB(24));  % 1200 Hz  R = [0.66 0.74 0.53] 34 spikes; rho=0.25;
JLbeatPlot(JbB_80dB(25));  % 1250 Hz  R = [0.78 0.67 0.55] 20 spikes; rho=0.23;
JLbeatPlot(JbB_80dB(26));  % 1300 Hz  R = [0.55 0.59 0.47] 30 spikes; rho=0.19;
JLbeatPlot(JbB_80dB(27));  % 1350 Hz  R = [0.5 0.7 0.35] 24 spikes; rho=0.18;
JLbeatPlot(JbB_80dB(28));  % 1400 Hz  R = [0.81 0.7 0.8] 17 spikes; rho=0.16;
JLbeatPlot(JbB_80dB(29));  % 1450 Hz  R = [0.31 0.65 0.15] 16 spikes; rho=0.13;
JLbeatPlot(JbB_80dB(30));  % 1500 Hz  R = [0.27 0.65 0.6] 12 spikes; rho=0.13;
JLbeatPlot(JbB_80dB(31));  % 1550 Hz  R = [0.57 0.62 0.47] 17 spikes; rho=0.12;
JLbeatPlot(JbB_80dB(32));  % 1600 Hz  R = [0.64 0.49 0.22] 14 spikes; rho=0.1;
JLbeatPlot(JbB_80dB(33));  % 1650 Hz  R = [0.37 0.48 0.25] 16 spikes; rho=0.086;
JLbeatPlot(JbB_80dB(34));  % 1700 Hz  R = [0.38 0.39 0.22] 9 spikes; rho=0.077;
JLbeatPlot(JbB_80dB(35));  % 1750 Hz  R = [0.4 0.47 0.52] 9 spikes; rho=0.063;
JLbeatPlot(JbB_80dB(36));  % 1800 Hz  R = [0.23 0.61 0.22] 14 spikes; rho=0.058;
JLbeatPlot(JbB_80dB(37));  % 1850 Hz  R = [0.31 0.44 0.42] 10 spikes; rho=0.046;
JLbeatPlot(JbB_80dB(38));  % 1900 Hz  R = [0.53 0.61 0.2] 7 spikes; rho=0.04;
JLbeatPlot(JbB_80dB(39));  % 1950 Hz  R = [0.094 0.22 0.11] 14 spikes; rho=0.021;
JLbeatPlot(JbB_80dB(40));  % 2000 Hz  R = [0.43 0.37 0.55] 9 spikes; rho=0.02;
% ===linearity
ifr = 1; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3);  % C too noisy, I noisy but ~consistent
% 100 Hz: strange I shape. Consistent B-I.
% analog waveformss? But note unclear filter settings during recording (JL's notes).
clampfit(JbC_80dB(2)); % ~analog  waveform in C ..
clampfit(JbI_80dB(1)); % .. but not in I
JLbeatPlot(JbI_80dB(2));
ifr = 2; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % C ~ok; I linear.
ifr = 3; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % ~linear
ifr = 3; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 6);  %spk addition test (see above remark on failingbinaurility)
% 150 Hz  R = [0.47 0.51 0.067] 212 spikes; rho=0.32;
% 150 Hz  R = [0.98 0.13 0.13] 80 spikes; rho=0.23;
% 150 Hz  R = [0.1 0.44 0.13] 46 spikes; rho=0.35;
ifr = 4; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % I better than C
ifr = 5; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % ditto
ifr = 6; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % ditto
ifr = 7; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % ditto
ifr = 8; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % ditto. 
% from here linearity is lost, but note JL's notes about losing the cell
ifr = 9; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % not very linear
ifr = 10; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % not very linear, either
ifr = 11; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % C is facilitated by C !?
% ===peak merging
JLbeatPlot({JbB_80dB(4) JbB_80dB(5) JbB_80dB(6) }, 3);
%
% RG10196 >>>> cell 4 ("EPSPs and Aps") 
%      60 dB bin
%      70 dB bin: SINKS
%      80 dB bin/ipsi/contra
JLreadBeats('RG10196', 4)
% many spikes
% ===80 dB
% 50 Hz: I noisy but consistent; C below thr
JLbeatPlot(JbB_80dB(1));  % 50 Hz  R = [0.49 0.027 0.056] 97 spikes; rho=0.099;
% 100 Hz: I&C consistent, and characterized by steep *downward* slopes. 
% Is this an artifact due to missed APs?
JLbeatPlot(JbB_80dB(2));  % 100 Hz  R = [0.74 0.31 0.089] 179 spikes; rho=0.36;
% 15-+ Hz. in previous cell, I&C phase locking seems to add @ low freqs: no
% coincidence detection, but simple pooling of "I & C spikes". 
% Note similarity btw I,C wavefeorms and corresponding CHISt.
JLbeatPlot(JbB_80dB(3));  % 150 Hz  R = [0.54 0.37 0.063] 355 spikes; rho=0.29; I strongly double-peaked
JLbeatPlot(JbB_80dB(4));  % 200 Hz  R = [0.37 0.6 0.089] 398 spikes; rho=0.41;
JLbeatPlot(JbB_80dB(5));  % 250 Hz  R = [0.3 0.74 0.16] 339 spikes; rho=0.3; % I doubl peaked
JLbeatPlot(JbB_80dB(6));  % 300 Hz  R = [0.51 0.6 0.25] 236 spikes; rho=0.38; % I & B double-peaked.
JLbeatPlot(JbB_80dB(7));  % 350 Hz  R = [0.71 0.49 0.37] 169 spikes; rho=0.34; multi peaks merge
JLbeatPlot(JbB_80dB(8));  % 400 Hz  R = [0.61 0.68 0.39] 168 spikes; rho=0.42; ditto
JLbeatPlot(JbB_80dB(9));  % 450 Hz  R = [0.68 0.71 0.45] 184 spikes; rho=0.39;
% 500+ Hz: peaks merge; many spikes; decent binaurality.
JLbeatPlot(JbB_80dB(10));  % 500 Hz  R = [0.73 0.78 0.54] 245 spikes; rho=0.39;
JLbeatPlot(JbB_80dB(11));  % 550 Hz  R = [0.76 0.78 0.56] 320 spikes; rho=0.38;
JLbeatPlot(JbB_80dB(12));  % 600 Hz  R = [0.73 0.82 0.57] 322 spikes; rho=0.38;
JLbeatPlot(JbB_80dB(13));  % 650 Hz  R = [0.71 0.78 0.51] 262 spikes; rho=0.39;
JLbeatPlot(JbB_80dB(14));  % 700 Hz  R = [0.74 0.78 0.53] 273 spikes; rho=0.39;
JLbeatPlot(JbB_80dB(15));  % 750 Hz  R = [0.71 0.77 0.52] 221 spikes; rho=0.4;
JLbeatPlot(JbB_80dB(16));  % 800 Hz  R = [0.68 0.74 0.51] 206 spikes; rho=0.38;
JLbeatPlot(JbB_80dB(17));  % 850 Hz  R = [0.67 0.69 0.44] 159 spikes; rho=0.37;
JLbeatPlot(JbB_80dB(18));  % 900 Hz  R = [0.68 0.71 0.52] 144 spikes; rho=0.33;
JLbeatPlot(JbB_80dB(19));  % 950 Hz  R = [0.67 0.59 0.39] 153 spikes; rho=0.27;
JLbeatPlot(JbB_80dB(20));  % 1000 Hz  R = [0.57 0.56 0.39] 143 spikes; rho=0.25;
% ===70 dB SINKS - skip
% ===60 dB: lower best freq
JLbeatPlot(JbB_60dB(1));  % 50 Hz  R = [0.76 0.58 0.007] 4 spikes; rho=0.046; too NOISY
JLbeatPlot(JbB_60dB(2));  % 100 Hz  R = [0.96 0.41 0.6] 3 spikes; rho=0.23; noisy, but consistent shapes
JLbeatPlot(JbB_60dB(3));  % 150 Hz  R = [0.97 0.58 0.63] 15 spikes; rho=0.41; nice I shape
JLbeatPlot(JbB_60dB(4));  % 200 Hz  R = [0.81 0.73 0.65] 10 spikes; rho=0.59; few, well-timed spks; nice waveforms
JLbeatPlot(JbB_60dB(5));  % 250 Hz  R = [0.64 0.89 0.58] 7 spikes; rho=0.49;
JLbeatPlot(JbB_60dB(6));  % 300 Hz  R = [0.85 0.8 0.86] 5 spikes; rho=0.56;
JLbeatPlot(JbB_60dB(7));  % 350 Hz  R = [0 0 0] 0 spikes; rho=0.52;
% (almost) no spikes >400 Hz, despite hig hrho
JLbeatPlot(JbB_60dB(8));  % 400 Hz  R = [0 0 0] 1 spikes; rho=0.6; high rho, but no spikes
JLbeatPlot(JbB_60dB(9));  % 450 Hz  R = [0 0 0] 0 spikes; rho=0.58;
JLbeatPlot(JbB_60dB(10));  % 500 Hz  R = [0 0 0] 0 spikes; rho=0.56;
JLbeatPlot(JbB_60dB(11));  % 550 Hz  R = [0 0 0] 0 spikes; rho=0.53;
JLbeatPlot(JbB_60dB(12));  % 600 Hz  R = [0 0 0] 0 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(13));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.35;
JLbeatPlot(JbB_60dB(14));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.3;
JLbeatPlot(JbB_60dB(15));  % 750 Hz  R = [0 0 0] 0 spikes; rho=0.29; % I remains 2-peaked
JLbeatPlot(JbB_60dB(16));  % 800 Hz  R = [0 0 0] 1 spikes; rho=0.27;
JLbeatPlot(JbB_60dB(17));  % 850 Hz  R = [0 0 0] 0 spikes; rho=0.16;
JLbeatPlot(JbB_60dB(18));  % 900 Hz  R = [0 0 0] 1 spikes; rho=0.13;
% finally spikes become sinks..see JL's note on losing the cell
JLbeatPlot(JbB_60dB(19));  % 950 Hz  R = [0.5 0.66 0.38] 3 spikes; rho=0.093; 
JLbeatPlot(JbB_60dB(20));  % 1000 Hz  R = [0.35 0.52 0.56] 11 spikes; rho=0.033;
% linearity 
% ---80 dB
ifr = 1; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % 50 Hz noise but I is consistent
ifr = 2; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % 100 Hz I&C lin
ifr = 3; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % 150 Hz I&C lin
ifr = 4; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % linear; I multipeaked
JLbeatPlot(JbI_80dB(4)); % 200 Hz  R = [0.78 0.06 0.089] 273 spikes; rho=0.27; mpeak also in I-stim, wv+CHIST
ifr = 5; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % linear; I multipeaked
JLbeatPlot(JbI_80dB(5)); % 250 Hz  R = [0.59 0.066 0.059] 167 spikes; rho=0.25;  mpeak also in I-stim, wv+CHIST
ifr = 5; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 6);
% 250 Hz  R = [0.3 0.74 0.16] 339 spikes; rho=0.3;
% 250 Hz  R = [0.59 0.066 0.059] 167 spikes; rho=0.25;
% 250 Hz  R = [0.048 0.84 0.027] 320 spikes; rho=0.2;
ifr = 6; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % linear
JLbeatPlot(JbI_80dB(6)); % 300 Hz  R = [0.68 0.042 0.051] 135 spikes; rho=0.22; mpeak also in I-stim, wv+CHIST
ifr = 7; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3);  % linear
JLbeatPlot(JbI_80dB(7));  % 350 Hz  R = [0.71 0.068 0.083] 126 spikes; rho=0.26; mpeaks merging
ifr = 8; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % C inconsistent
JLbeatPlot(JbI_80dB(8));  % 400 Hz  R = [0.77 0.081 0.061] 121 spikes; rho=0.32;
ifr = 10; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 3); % mono peaks; reduced linearity 
ifr=10; JLbeatPlot({JbI_80dB(ifr) JbC_80dB(ifr) JbB_80dB(ifr) },6);
% 500 Hz  R = [0.87 0.049 0.039] 239 spikes; rho=0.4;
% 500 Hz  R = [0.011 0.63 0.055] 206 spikes; rho=0.35;
% 500 Hz  R = [0.73 0.78 0.54] 245 spikes; rho=0.39;
ifr = 11; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 6);
% 550 Hz  R = [0.76 0.78 0.56] 320 spikes; rho=0.38;
% 550 Hz  R = [0.92 0.011 0.032] 368 spikes; rho=0.45;
% 550 Hz  R = [0.056 0.71 0.052] 304 spikes; rho=0.36;
ifr = 15; JLbeatPlot({JbB_80dB(ifr) JbI_80dB(ifr) JbC_80dB(ifr) }, 6);
% 750 Hz  R = [0.71 0.77 0.52] 221 spikes; rho=0.4;
% 750 Hz  R = [0.86 0.043 0.04] 370 spikes; rho=0.46;
% 750 Hz  R = [0.026 0.67 0.062] 249 spikes; rho=0.33;
qq=[]

function RG10197;
% RG10197 >>>> cell 5 ("big spikes")
%      70 dB 2xbin/ipsi/contra
% -5   5-1-BFS 50*:50*:1000 Hz 70 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -6   5-2-BFS 50*:50*:1000 Hz -70|70 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -7   5-3-BFS 50*:50*:1000 Hz 70|-70 dB 4 Hz beat 20 x 1 x 6000/8500 ms  B
% -8   5-4-BFS 50*:50*:200 Hz  70 dB     4 Hz beat 4 x 1 x 6000/8500 ms   B
% FS
JLreadBeats('RG10197', 5);
% ===70 dB  
% Onset-like responses. Spikes occur early in trace.
% 50 Hz: spikes are I-locked, but I wv is noise?!
JLbeatPlot(JbB_70dB(1));  % 50 Hz  R = [0.78 0.16 0.14] 38 spikes; rho=0.017;
JLbeatPlot({JbB_70dB(1) JbB_70dB(1) JbB2_70dB(1) }, 3); % .. it's less noisy than it looks ..
JLbeatPlot(JbB_70dB(2));  % 100 Hz  R = [0.57 0.062 0.86] 2 spikes; rho=0.14;
JLbeatPlot(JbB_70dB(3));  % 150 Hz  R = [0 0 0] 1 spikes; rho=0.051;
% 200+ Hz wave start to look better - no spikes though
JLbeatPlot(JbB_70dB(4));  % 200 Hz  R = [0 0 0] 0 spikes; rho=0.13;
JLbeatPlot(JbB_70dB(5));  % 250 Hz  R = [0 0 0] 0 spikes; rho=0.083;
clampfit(JbB_70dB(5)); % onset-like response
JLbeatPlot(JbB_70dB(6));  % 300 Hz  R = [0 0 0] 0 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(7));  % 350 Hz  R = [0 0 0] 0 spikes; rho=0.21;
JLbeatPlot(JbB_70dB(8));  % 400 Hz  R = [0 0 0] 1 spikes; rho=0.28;
JLbeatPlot(JbB_70dB(9));  % 450 Hz  R = [0 0 0] 0 spikes; rho=0.25;
JLbeatPlot(JbB_70dB(10)); % 500 Hz  R = [0 0 0] 0 spikes; rho=0.27;
JLbeatPlot(JbB_70dB(11)); % 550 Hz  R = [0 0 0] 1 spikes; rho=0.25;
JLbeatPlot(JbB_70dB(12)); % 600 Hz  R = [0 0 0] 1 spikes; rho=0.24;
% higher freqs: still hardly any spikes; multi peaks in C
JLbeatPlot(JbB_70dB(13)); % 650 Hz  R = [0.95 0.6 0.5] 3 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(14)); % 700 Hz  R = [0.63 1 0.58] 2 spikes; rho=0.14;
JLbeatPlot(JbB_70dB(15)); % 750 Hz  R = [0.74 0.89 0.96] 3 spikes; rho=0.13;
JLbeatPlot(JbB_70dB(16));  % 800 Hz  R = [0.5 0.53 0.76] 4 spikes; rho=0.12;
JLbeatPlot(JbB_70dB(17));  % 850 Hz  R = [0.55 0.42 0.31] 7 spikes; rho=0.13;
JLbeatPlot(JbB_70dB(18));  % 900 Hz  R = [0.25 0.41 0.37] 6 spikes; rho=0.13;
JLbeatPlot(JbB_70dB(19));  % 950 Hz  R = [0.13 0.53 0.56] 6 spikes; rho=0.16;
JLbeatPlot(JbB_70dB(20));  % 1000 Hz  R = [0.32 0.22 0.47] 9 spikes; rho=0.19;
% consistency across repeats
ifr = 1; JLbeatPlot({JbB_70dB(ifr) JbB_70dB(ifr) JbB2_70dB(ifr)}, 3); % noisy waveforms, but I look alike
ifr = 2; JLbeatPlot({JbB_70dB(ifr) JbB_70dB(ifr) JbB2_70dB(ifr)}, 3); % very consistent
ifr = 3; JLbeatPlot({JbB_70dB(ifr) JbB_70dB(ifr) JbB2_70dB(ifr)}, 3); % ditto
ifr = 4; JLbeatPlot({JbB_70dB(ifr) JbB_70dB(ifr) JbB2_70dB(ifr)}, 3); % almost exact copy
% linearity
% ---70 dB
ifr = 1; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 50 Hz; noisy, but I wv is consistent
ifr = 2; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 100 Hz; consistent shapes; C suppr
ifr = 3; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % small variations & reduction
ifr = 4; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % amazingly identical
ifr = 5; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % C identical; I reduced and ~changed
ifr = 6; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % linear, maye slight reduction
ifr = 7; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % linear
ifr = 8; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3);  %linear
ifr = 9; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % linear; C somewhat facilitated
ifr = 10; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % linear
ifr = 11; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % linear
ifr = 12; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3);  %linear
ifr = 13; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin
ifr = 14; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin
ifr = 15; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin
ifr = 16; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin; maybe slight I reduction
ifr = 17; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin, but I reduced
ifr = 18; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % C lin, I affected
ifr = 19; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % ditto
ifr = 20; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % lin,, but some ampl changes
% RG10197 >>>> cell 8 ("big spikes")
%      70 dB 2xbin/ipsi/contra
%      FS
% -19  8-1-BFS 50*:50*:1000 Hz 70 dB     4 Hz beat 20 x 1 x 6000/8500 ms  B
% -20  8-2-BFS 100*:50*:600 Hz 70|-70 dB 4 Hz beat 11 x 1 x 6000/8500 ms  B
% -21  8-3-BFS 100*:50*:600 Hz -70|70 dB 4 Hz beat 11 x 1 x 6000/8500 ms  B
% -22  8-4-BFS 100*:50*:600 Hz 70 dB     4 Hz beat 11 x 1 x 6000/8500 ms  B
% FS & some additional 70-dB Bin beats (SGSR data missing)
JLreadBeats('RG10197', 8);
% ===FS; 70 dB
ifr=10; % 102 Hz; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f1; ifr=15; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f2; ifr=16; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f3; ifr=17; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f4; ifr=18; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f5; ifr=19; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
f6; ifr=20; 
JLtone2('RG10197','8-5-FS', ifr, 0, 'k'); % bin
JLtone2('RG10197','8-6-FS', ifr, 0, 'r');  % contra
JLtone2('RG10197','8-7-FS', ifr, 0, 'b'); % ipsi
delete(findobj('type','axes', 'tag', 'legend'))
subsaf; linkaxes(findobj('type','axes', 'tag', ''));
% ===70 dB  
JLbeatPlot(JbB_70dB(1));  % 50 Hz  R = [0 0 0] 0 spikes; rho=-0.0011;
JLbeatPlot(JbB_70dB(2));  % 100 Hz  R = [0 0 0] 0 spikes; rho=0.31; C inhib? I looks "analog"
JLbeatPlot(JbB_70dB(3));  % 150 Hz  R = [0 0 0] 0 spikes; rho=0.24; C looks inhib again
JLbeatPlot(JbB_70dB(4));  % 200 Hz  R = [0.98 0.63 0.53] 8 spikes; rho=0.77; I analog, but spks well I-locked..
JLbeatPlot(JbB_70dB(5));  % 250 Hz  R = [0.96 0.83 0.76] 55 spikes; rho=0.49; STAs?
JLbeatPlot(JbB_70dB(6));  % 300 Hz  R = [0.95 0.8 0.73] 110 spikes; rho=0.67;
% For all of these recs, it is hard to reconcile the small magnitude of the
% ipsi waveform with its contrib to evoke well-timed spikes.
% At these freqs, the C contrib stems from negative dips in the raw traces.
% the I (which evokes APs!) doesn't show much:
clampfit(JbI_70dB(5)); 
Db = readTKABF(JbB_70dB(6)); 
Dc = readTKABF(JbC_70dB(5));
subplot(3,1,3); delete(gca); subplot(3,1,2); delete(gca); 
set(gca, 'posi', [0.1300    0.1093    0.7750    0.7155]);
ylim([-1 1.5]); set(gcf,'units', 'normalized', 'position', [0.0227 0.364 0.944 0.543]);
xdplot(Db.dt_ms, Db.AD(1).samples, 'k')
xdplot(Dc.dt_ms, Dc.AD(1).samples, 'r') % note the similarity between B & C.
xlim(500+[-15 15]);
%
JLbeatPlot(JbB_70dB(7));  % 350 Hz  R = [0.96 0.82 0.77] 77 spikes; rho=0.58;
JLbeatPlot(JbB_70dB(8));  % 400 Hz  R = [0.89 0.89 0.8] 39 spikes; rho=0.63;
JLbeatPlot(JbB_70dB(9));  % 450 Hz  R = [0.93 0.83 0.8] 15 spikes; rho=0.54;
JLbeatPlot(JbB_70dB(10));  % 500 Hz  R = [0.94 0.83 0.79] 9 spikes; rho=0.48;
JLbeatPlot(JbB_70dB(11)); % 550 Hz  R = [0 0 0] 1 spikes; rho=0.4; I wv stays small
JLbeatPlot(JbB_70dB(12));  % 600 Hz  R = [0 0 0] 0 spikes; rho=0.37;
JLbeatPlot(JbB_70dB(13));  % 650 Hz  R = [0 0 0] 1 spikes; rho=0.32;
JLbeatPlot(JbB_70dB(14));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.3;
JLbeatPlot(JbB_70dB(15));  % 750 Hz  R = [0 0 0] 0 spikes; rho=0.24;
JLbeatPlot(JbB_70dB(16));  % 800 Hz  R = [0 0 0] 0 spikes; rho=0.24;
% I wv grows somewhat, but remains markedly different from C. I is
% more "rectified" than C)
JLbeatPlot(JbB_70dB(17));  % 850 Hz  R = [0 0 0] 0 spikes; rho=0.21;
JLbeatPlot(JbB_70dB(18));  % 900 Hz  R = [0 0 0] 0 spikes; rho=0.22;
JLbeatPlot(JbB_70dB(19));  % 950 Hz  R = [0 0 0] 0 spikes; rho=0.2;
JLbeatPlot(JbB_70dB(20));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.16;
% linearity
% --70 dB
ifr = 1; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 100 Hz; smooth I wv. C has single neg pk
ifr = 2; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 150 Hz; very linear
ifr = 3; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 200 Hz; very smooth I wv.
ifr = 4; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); % 250 Hz; I wv gets "dented"; C has neg pk
% --300 Hz; linear. But strange waveforms & spike timing re waveforms.
ifr = 5; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % 300 Hz; I split; small + C pk precedes large - pk
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 5); % spks timed on 2nd, larger, I pk
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 6); % I hist not split
% B 300 Hz  R = [0.95 0.8 0.73] 110 spikes; rho=0.67; facilitation
% I 300 Hz  R = [0.96 0.15 0.1] 24 spikes; rho=0.25;
% C 300 Hz  R = [0 0 0] 0 spikes; rho=0.71;
% --350 Hz: % comparable to 300 Hz
ifr = 6; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); 
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 5); % ditto
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 6); % ditto
% 350 Hz  R = [0.96 0.82 0.77] 77 spikes; rho=0.57; facilitation again
% 350 Hz  R = [0.97 0.14 0.13] 25 spikes; rho=0.18;
% 350 Hz  R = [0 0 0] 0 spikes; rho=0.56;
% --400 Hz
ifr = 7; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % comparable to 300 Hz
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 5); % ditto
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 6); % ditto
% 400 Hz  R = [0.89 0.89 0.8] 39 spikes; rho=0.63;
% 400 Hz  R = [0.93 0.02 0.0086] 35 spikes; rho=0.17;
% 400 Hz  R = [0 0 0] 0 spikes; rho=0.56;
% --450 Hz .. less spikes but still the same overall picture
ifr = 8; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); 
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 5); % ditto
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 6); % ditto
% 450 Hz  R = [0.93 0.83 0.8] 15 spikes; rho=0.54;
% 450 Hz  R = [0.95 0.068 0.049] 20 spikes; rho=0.11;
% 450 Hz  R = [0 0 0] 0 spikes; rho=0.41;
% --500 Hz .. less spikes but still the same overall picture
ifr = 9; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); 
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 5); % ditto
JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 6); % ditto
% 500 Hz  R = [0.94 0.83 0.79] 9 spikes; rho=0.48;
% 500 Hz  R = [0.95 0.14 0.088] 16 spikes; rho=0.066;
% 500 Hz  R = [0 0 0] 0 spikes; rho=0.45;
ifr = 10; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % linear; I wv very small
ifr = 11; JLbeatPlot({JbB_70dB(ifr+1) JbI_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % linear; I wv very small
qq=[];

function RG10198;
% RG10198 >>>> cell 1 ("spikes + EPSPs")  
%      70 dB bin/contra
% -1   1-2-BFS  0.1*:0.1*:2 kHz  70 dB     4 Hz beat               20 x 1 x 6000/8500 ms B
% -2   1-3-BFS  0.1*:0.1*:2 kHz  -70|70 dB 4 Hz beat               20 x 1 x 6000/8500 ms B
%  1-2-BFS % one missing abf screws up the bookkeeping; only read up to 1100 Hz
JLreadBeats('RG10198', 1);
% ==70 dB
JLbeatPlot(JbB_70dB(1));  % 100 Hz  R = [0.28 0.93 0.31] 27 spikes; rho=0.71; ~mon peaked C wv
JLbeatPlot(JbB_70dB(2));  % 200 Hz  R = [0.52 0.74 0.39] 12 spikes; rho=0.57;
JLbeatPlot(JbB_70dB(3));  % 300 Hz  R = [0.72 0.84 0.59] 38 spikes; rho=0.43; I wv multi peaked
JLbeatPlot(JbB_70dB(4));  % 400 Hz  R = [0.78 0.93 0.65] 101 spikes; rho=0.38; strange timing I STA
JLbeatPlot(JbB_70dB(5));  % 500 Hz  R = [0.69 0.93 0.6] 121 spikes; rho=0.43; still delayed I timing in STA
JLbeatPlot(JbB_70dB(6));  % 600 Hz  R = [0.58 0.89 0.56] 60 spikes; rho=0.52; I timing better, but still > C
JLbeatPlot(JbB_70dB(7));  % 700 Hz  R = [0.46 0.84 0.41] 29 spikes; rho=0.44; same
JLbeatPlot(JbB_70dB(8));  % 800 Hz  R = [0.57 0.75 0.38] 28 spikes; rho=0.29; now I timing precedes C
JLbeatPlot(JbB_70dB(9));  % 900 Hz  R = [0.23 0.82 0.12] 18 spikes; rho=0.26; I timing is poor -> I-STA reduced
JLbeatPlot(JbB_70dB(10));  % 1000 Hz  R = [0.81 0.79 0.77] 11 spikes; rho=0.24; revival of I timing, but few spks
JLbeatPlot(JbB_70dB(11));  % 1100 Hz  R = [0.66 0.84 0.65] 13 spikes; rho=0.21; same
% linearity: note that I-alone is missing
% --70 dB
ifr = 1; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % C wv very stable
ifr = 2; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % ditto
ifr = 3; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % C "facilitated"
ifr = 4; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % C "suppressed"
ifr = 5; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % curious inversion of pk sizes
ifr = 6; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % diff in timing +76 us
ifr = 7; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % + 181 us
ifr = 8; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % + 100 us
ifr = 9; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % +49 us
ifr = 10; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % +6 us
ifr = 11; JLbeatPlot({JbB_70dB(ifr) JbC_70dB(ifr) JbC_70dB(ifr) }, 3, 'x3'); % + 10 us
% RG10198 >>>> cell 2 ("spikes + EPSPs + heartbeat")  
%      40 dB 3xbin/ipsi/2xcontra
%      70 dB 2xbin/ipsi/2xcontra
JLreadBeats('RG10198', 1);
% ==70 dB
% 50 Hz
JLbeatPlot(JbB2_70dB(1));  % 50 Hz  R = [0.28 0.48 0.055] 12 spikes; rho=0.017;
% 100 Hz; I small, but consistent; C wv large, 2 pks, APs @ 1st, larger pk
JLbeatPlot([JbB_70dB(1) JbB2_70dB(2)]); 
JLbeatPlot({JbB_70dB(1) JbB_70dB(1) JbB2_70dB(2) }, 3, 'x3');  % good consistency
% 150 Hz
JLbeatPlot(JbB2_70dB(3));  % 150 Hz  R = [0 0 0] 1 spikes; rho=0.25;
% 200 Hz
JLbeatPlot([JbB_70dB(2) JbB2_70dB(4)]);  % 200 Hz  R = [0.99 0.99 0.98] 4 spikes; rho=0.42;
JLbeatPlot({JbB_70dB(2) JbB_70dB(2) JbB2_70dB(4) }, 3, 'x3');  % good consistency
% 250 Hz
JLbeatPlot(JbB2_70dB(5));  % 250 Hz  R = [0.87 0.97 0.74] 3 spikes; rho=0.32; C&I multipeaked 
% 300 Hz
JLbeatPlot([JbB_70dB(3) JbB2_70dB(6)]);  % 300 Hz  R = [0.86 1 0.86] 3 spikes; rho=0.4;
JLbeatPlot({JbB_70dB(3) JbB_70dB(3) JbB2_70dB(6) }, 3, 'x3'); % amazingly identical
% 350 Hz
JLbeatPlot(JbB2_70dB(7));  % 350 Hz  R = [0.92 1 0.9] 3 spikes; rho=0.39;
% 400 Hz
JLbeatPlot([JbB_70dB(4) JbB2_70dB(8)]);  % 400 Hz  R = [0.83 0.85 0.72] 9 spikes; rho=0.48; 2-pkd I, also in HIST
JLbeatPlot({JbB_70dB(4) JbB_70dB(4) JbB2_70dB(8) }, 3, 'x3'); % amazingly identical
% 450 Hz
JLbeatPlot(JbB2_70dB(9));  % 450 Hz  R = [0.96 0.9 0.86] 9 spikes; rho=0.46; secondary I peak shrinking
% 500 Hz
JLbeatPlot([JbB_70dB(5) JbB2_70dB(10)]);  % 500 Hz  R = [0.94 0.9 0.83] 16 spikes; rho=0.42;
JLbeatPlot({JbB_70dB(5) JbB_70dB(5) JbB2_70dB(10) }, 3, 'x3'); % very consistent
% 550 Hz
JLbeatPlot(JbB2_70dB(11));  % 550 Hz  R = [0.91 0.9 0.78] 10 spikes; rho=0.42;  I&C ~mono peaked; I STA timing late
% 600 Hz
JLbeatPlot([JbB_70dB(6) JbB2_70dB(12)]);  % 600 Hz  R = [0.91 0.92 0.85] 13 spikes; rho=0.46; mono pk; good binaurality
JLbeatPlot({JbB_70dB(6) JbB_70dB(6) JbB2_70dB(12) }, 3, 'x3'); % very consistent
% 650 Hz
JLbeatPlot(JbB2_70dB(13));  % 650 Hz  R = [0.84 0.81 1] 2 spikes; rho=0.48;
% 700 Hz
JLbeatPlot([JbB_70dB(7) JbB2_70dB(14)]);  % 700 Hz  R = [0.69 0.82 0.67] 7 spikes; rho=0.49;
JLbeatPlot({JbB_70dB(7) JbB_70dB(7) JbB2_70dB(14) }, 3, 'x3'); % virtually identical
% 750 Hz
JLbeatPlot(JbB2_70dB(15));  % 750 Hz  R = [0.99 1 1] 2 spikes; rho=0.45;
% 800 Hz
JLbeatPlot([JbB_70dB(8) JbB2_70dB(16)]);  % 800 Hz  R = [0.86 0.62 0.61] 3 spikes; rho=0.38;
JLbeatPlot({JbB_70dB(8) JbB_70dB(8) JbB2_70dB(16) }, 3, 'x3'); % virtually identical 
% 850 Hz
JLbeatPlot(JbB2_70dB(17));  % 850 Hz  R = [0.99 0.93 0.98] 2 spikes; rho=0.32;
% 900 Hz
JLbeatPlot([JbB_70dB(9) JbB2_70dB(18)]);  % 900 Hz  R = [0 0 0] 1 spikes; rho=0.27;
JLbeatPlot({JbB_70dB(9) JbB_70dB(9) JbB2_70dB(18) }, 3, 'x3'); % virtually identical 
% 950 Hz
JLbeatPlot(JbB2_70dB(19));  % 950 Hz  R = [0 0 0] 1 spikes; rho=0.22;
% 1000 Hz
JLbeatPlot([JbB_70dB(10) JbB2_70dB(20)]);  % 1000 Hz  R = [0 0 0] 1 spikes; rho=0.2;
JLbeatPlot({JbB_70dB(10) JbB_70dB(10) JbB2_70dB(20) }, 3, 'x3'); % virtually identical 
% 1100 Hz and up
JLbeatPlot(JbB_70dB(11));  % 1100 Hz  R = [0 0 0] 0 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(12));  % 1200 Hz  R = [0 0 0] 1 spikes; rho=0.12;
JLbeatPlot(JbB_70dB(13));  % 1300 Hz  R = [0 0 0] 1 spikes; rho=0.098;
JLbeatPlot(JbB_70dB(14));  % 1400 Hz  R = [0 0 0] 0 spikes; rho=0.062;
JLbeatPlot(JbB_70dB(15));  % 1500 Hz  R = [0 0 0] 1 spikes; rho=0.04;
JLbeatPlot(JbB_70dB(16));  % 1600 Hz  R = [0 0 0] 0 spikes; rho=0.03;
JLbeatPlot(JbB_70dB(17));  % 1700 Hz  R = [1 0.87 0.89] 2 spikes; rho=0.02;
% 1800 Hz & up: spikes look spontaneous not driven
JLbeatPlot(JbB_70dB(18));  % 1800 Hz  R = [0.13 0.37 0.44] 6 spikes; rho=0.011;
JLbeatPlot(JbB_70dB(19));  % 1900 Hz  R = [0.23 0.25 0.29] 7 spikes; rho=0.0034;
JLbeatPlot(JbB_70dB(20));  % 2000 Hz  R = [0.25 0.36 0.31] 14 spikes; rho=0.0022;
% ==40 dB
% 50 Hz: noise
ifr=1; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]);% 50 Hz R = [0.23 0.13 0.11] 18 spikes; rho=-0.0021;
% 100 Hz: Analog-like waveforms, APs hardly locked.
ifr=2; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]);% 100 Hz R = [0.4 0.078 0.027] 21 spikes; rho=0.072;
% 150 Hz: Analog-like waveforms, APs hardly locked. (field?)
ifr=3; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]);% 150 Hz R=[0.05 0.11 0.018] 25 spikes; rho=0.095;
% 200 Hz: I: Analog-like waveform (field?); C looks more significant.
ifr=4; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]);% 200 Hz R=[0.16 0.28 0.1] 18 spikes; rho=0.2;
% 250 Hz: C&I wvs look less filed-like, but APs & STAs are still messy.
ifr=5; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]); % 250 Hz R = [0.28 0.46 0.3] 22 spks; rho=0.11;
% 300 Hz; Now the driving starts. Decent binaurality.
ifr=6; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]); % 300 Hz R = [0.7 0.64 0.55] 32 spks; rho=0.15;
% 350 Hz: good binaurality. Small I wv contradicts good I phase locking.
ifr=7; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]); % 350 Hz R = [0.94 0.83 0.79] 32 spks; rho=0.25;
% 400 Hz: very binaural. I wv period doubling is not consistent w I-CHIST & R. The B-STA peaks ~400 us before 0.
ifr=8; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]); % 400 Hz R = [0.94 0.93 0.88] 22 spks; rho=0.33;
% 450 Hz: very binaural, same inconsistencies as 400 Hz, though. 
ifr=9; JLbeatPlot([JbB_40dB(ifr) JbB2_40dB(ifr) JbB3_40dB(ifr)]);  % 450 Hz  R = [0.89 0.9 0.78] 29 spikes; rho=0.28;

function RG10201;
% RG10201 cell 1 ("spikes + EPSPs")  
%      60 dB bin/ipsi
% -1   1-1-BFS  250*:150*:2500 Hz 60 dB     4 Hz beat 16 x 1 x 6000/8500 ms B
% -2   1-2-BFS  250*:150*:2500 Hz 60|-60 dB 4 Hz beat 16 x 1 x 6000/8500 ms B
JLreadBeats('RG10201', 1);
% ==60 dB
JLbeatPlot(JbB_60dB(1));  % 250 Hz  R = [0.67 0.95 0.59] 5 spikes; rho=0.43; % C&I multi pks
JLbeatPlot(JbB_60dB(2));  % 400 Hz  R = [0.95 0.65 0.6] 12 spikes; rho=0.62; % ~single pks
JLbeatPlot(JbB_60dB(3));  % 550 Hz  R = [0.92 0.86 0.99] 2 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(4));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.4;
JLbeatPlot(JbB_60dB(5));  % 850 Hz  R = [0 0 0] 0 spikes; rho=0.42;
JLbeatPlot(JbB_60dB(6));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.34;
% linearity (I only)
ifr = 1; JLbeatPlot({JbB_60dB(ifr) JbB_60dB(ifr) JbI_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbB_60dB(ifr) JbI_60dB(ifr)}, 3, 'x3'); % ~linear
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbB_60dB(ifr) JbI_60dB(ifr)}, 3, 'x3'); % ~linear
%
% RG10201 >>>>>>>cell 2 ("spikes + EPSPs")  
%      40 dB 2xbin/ipsi/3xcontra
%      60 dB bin [xt freq]
%      60 dB bin/ipsi/contra
JLreadBeats('RG10201', 2);
% ==60 dB
JLbeatPlot(JbB_60dB(1));  % 50 Hz  R = [0 0 0] 1 spikes; rho=0.033;  noise
JLbeatPlot(JbB_60dB(2));  % 150 Hz  R = [0.99 0.63 0.62] 18 spikes; rho=0.43;; I mult pk (not hist); I dominates
JLbeatPlot(JbB_60dB(3));  % 250 Hz  R = [0.98 0.58 0.56] 15 spikes; rho=0.54; contradiction 2-pk I wv <> 1-pk hist
JLbeatPlot(JbB_60dB(4));  % 350 Hz  R = [0.98 0.86 0.86] 7 spikes; rho=0.46; still somewhat multipeaked
JLbeatPlot(JbB_60dB(5));  % 450 Hz  R = [0.98 0.86 0.86] 5 spikes; rho=0.57; I dominates STA 
% spikes disappear; rho still pretty high
JLbeatPlot(JbB_60dB(6));  % 550 Hz  R = [1 1 1] 2 spikes; rho=0.55;
JLbeatPlot(JbB_60dB(7));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.54;
JLbeatPlot(JbB_60dB(8));  % 750 Hz  R = [0 0 0] 1 spikes; rho=0.46;
JLbeatPlot(JbB_60dB(9));  % 850 Hz  R = [0 0 0] 0 spikes; rho=0.43;
JLbeatPlot(JbB_60dB(10));  % 950 Hz  R = [0 0 0] 0 spikes; rho=0.31;
JLbeatPlot(JbB_60dB(11));  % 1050 Hz  R = [0 0 0] 0 spikes; rho=0.22;
% ==60 dB xtende freq range
%   XXXXXXXXXXXXXXXXXXXXXXXXXXXX
%   XXXXXXXXXXXXXXXXXXXXXXXXXXXX
%   XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ==40 dB
JLbeatPlot(JbB_40dB(1));  % 50 Hz  R = [0 0 0] 0 spikes; rho=-0.0024;  NOISE
JLbeatPlot(JbB_40dB(2));  % 150 Hz  R = [0.72 0.86 0.97] 2 spikes; rho=0.0041; still noisy
JLbeatPlot(JbB_40dB(3));  % 250 Hz  R = [0 0 0] 0 spikes; rho=0.027;
JLbeatPlot(JbB_40dB(4));  % 350 Hz  R = [0.91 0.78 0.8] 8 spikes; rho=0.37; decent bin, 1-pkd I&C, steep right flank
JLbeatPlot(JbB_40dB(5));  % 450 Hz  R = [0.91 0.81 0.7] 8 spikes; rho=0.57; % ditto
JLbeatPlot(JbB_40dB(6));  % 550 Hz  R = [0.9 0.96 0.94] 5 spikes; rho=0.54; a few excellently timed spikes
% no spikes, still good rho values
JLbeatPlot(JbB_40dB(7));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.47;
JLbeatPlot(JbB_40dB(8));  % 750 Hz  R = [0 0 0] 1 spikes; rho=0.44;
JLbeatPlot(JbB_40dB(9));  % 850 Hz  R = [0.99 0.91 0.84] 2 spikes; rho=0.39;
JLbeatPlot(JbB_40dB(10));  % 950 Hz  R = [0 0 0] 0 spikes; rho=0.29;
JLbeatPlot(JbB_40dB(11));  % 1050 Hz  R = [0 0 0] 0 spikes; rho=0.21;
%========linearity====
% --60 dB (note JL's notes on stability)
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % very linear I dominates
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 6, 'x3');
% 150 Hz  R = [0.99 0.63 0.62] 18 spikes; rho=0.43;
% 150 Hz  R = [1 0.46 0.43] 9 spikes; rho=0.5;
% 150 Hz  R = [0 0 0] 0 spikes; rho=0.047;
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3');  %same
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 6, 'x3');
% 250 Hz  R = [0.98 0.58 0.56] 15 spikes; rho=0.54;
% 250 Hz  R = [0.99 0.19 0.15] 5 spikes; rho=0.56;
% 250 Hz  R = [0.3 0.98 0.13] 2 spikes; rho=0.13;
ifr = 4; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % stable shapes; some suppr
ifr = 4; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 6, 'x3');
% 350 Hz  R = [0.98 0.86 0.86] 7 spikes; rho=0.46;
% 350 Hz  R = [0.95 0.51 0.75] 2 spikes; rho=0.47;
% 350 Hz  R = [0.62 0.67 0.72] 5 spikes; rho=0.12;
ifr = 5; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % stable shapes; some suppr
ifr = 5; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 6, 'x3');
% 450 Hz  R = [0.98 0.86 0.86] 5 spikes; rho=0.57;
% 450 Hz  R = [0 0 0] 0 spikes; rho=0.52;
% 450 Hz  R = [0.72 0.97 0.78] 3 spikes; rho=0.27;
ifr = 6; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % some shape change
ifr = 6; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 6, 'x3');
% 550 Hz  R = [1 1 1] 2 spikes; rho=0.55;
% 550 Hz  R = [0 0 0] 0 spikes; rho=0.39;
% 550 Hz  R = [0 0 0] 1 spikes; rho=0.25;
ifr = 7; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % shape changes
ifr = 8; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % same
ifr = 9; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % big changes, also in rel ampl
ifr = 10; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % same
ifr = 11; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr)}, 3, 'x3'); % C strongly reduced; I ~stable
% --40 dB note JL's notes on stability
ifr = 2; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr)}, 3, 'x3'); % noisy
ifr = 3; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr)}, 3, 'x3');
% instability affects relative peak sizes ... 
ifr = 3; JLbeatPlot({JbB_40dB(ifr) JbB2_40dB(ifr)  JbB2_40dB(ifr) }, 3, 'x3');

% RG10201 >>>> cell 3 ("spikes + EPSPs")  
%      40 dB bin/contra
%      60 dB bin
JLreadBeats('RG10201', 3);
% ===60 dB  Not many spikes, but well-timed ones.
JLbeatPlot(JbB_60dB(1));  % 50 Hz  R = [0 0 0] 0 spikes; rho=0.027; Noise
JLbeatPlot(JbB_60dB(2));  % 150 Hz  R = [0 0 0] 0 spikes; rho=0.18;  rugged waveforms; no spikes
JLbeatPlot(JbB_60dB(3));  % 250 Hz  R = [0.99 0.92 0.88] 6 spikes; rho=0.32;few APs but well timed; I mult pk, HIST 1-pk
JLbeatPlot(JbB_60dB(4));  % 350 Hz  R = [0.98 0.97 0.98] 3 spikes; rho=0.45; same story
JLbeatPlot(JbB_60dB(5));  % 450 Hz  R = [0.98 0.98 1] 2 spikes; rho=0.48; same
JLbeatPlot(JbB_60dB(6));  % 550 Hz  R = [1 1 1] 2 spikes; rho=0.46; same
JLbeatPlot(JbB_60dB(7));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.4; no APs at all
JLbeatPlot(JbB_60dB(8));  % 750 Hz  R = [1 0.92 0.91] 2 spikes; rho=0.36;
JLbeatPlot(JbB_60dB(9));  % 850 Hz  R = [0 0 0] 0 spikes; rho=0.37;
JLbeatPlot(JbB_60dB(10));  % 950 Hz  R = [0 0 0] 0 spikes; rho=0.35;
JLbeatPlot(JbB_60dB(11));  % 1050 Hz  R = [0 0 0] 0 spikes; rho=0.3;
JLbeatPlot(JbB_60dB(12));  % 1150 Hz  R = [0 0 0] 0 spikes; rho=0.27;
JLbeatPlot(JbB_60dB(13));  % 1250 Hz  R = [0 0 0] 0 spikes; rho=0.22;
% nothing much happens at still higher freqs
% ===40 dB  almost no APs. C mostly dominates. Practically mono pkd.
JLbeatPlot(JbB_40dB(1));  % 50 Hz  R = [0 0 0] 0 spikes; rho=0.012;
JLbeatPlot(JbB_40dB(2));  % 150 Hz  R = [0 0 0] 0 spikes; rho=0.027; noisy wvs, no APs, C dominates
JLbeatPlot(JbB_40dB(3));  % 250 Hz  R = [0 0 0] 0 spikes; rho=0.05; I~=C still noisy
JLbeatPlot(JbB_40dB(4));  % 350 Hz  R = [0 0 0] 0 spikes; rho=0.019; I~=C. Somewhat multi pkd
JLbeatPlot(JbB_40dB(5));  % 450 Hz  R = [0 0 0] 0 spikes; rho=0.2; C dominates again
JLbeatPlot(JbB_40dB(6));  % 550 Hz  R = [0 0 0] 1 spikes; rho=0.43; 
JLbeatPlot(JbB_40dB(7));  % 650 Hz  R = [0 0 0] 0 spikes; rho=0.5; good rho ut no APs
JLbeatPlot(JbB_40dB(8));  % 750 Hz  R = [0 0 0] 1 spikes; rho=0.42;
JLbeatPlot(JbB_40dB(10));  % 950 Hz  R = [0 0 0] 0 spikes; rho=0.34;
JLbeatPlot(JbB_40dB(13));  % 1250 Hz  R = [0 0 0] 0 spikes; rho=0.18; I dominates; rho decreases
% nothing much happens at still higher freqs
% ===linearity (contra only)
ifr = 3; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % noisy; ~lin
ifr = 4; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % same
ifr = 5; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % stable shape ; C ampl enhanced
ifr = 6; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % same
ifr = 7; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % pretty linear
ifr = 8; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % linear
ifr = 9; JLbeatPlot({JbB_40dB(ifr) JbB_40dB(ifr)  JbC_40dB(ifr) }, 3, 'x3'); % ditto



function RG10204
% RG10204 >>> cell 10 ("spikes + EPSPs")  
%      50 dB bin/ipsi/contra
%      60 dB 3xbin/2xipsi/2xcontra
%      70 dB bin/ipsi/contra
%      71 dB 6xbin/7xipsi/6xcontra (70 dB, 5 freqs, 300-700 Hz)
%      FS, NTD, BN, FM, CFS
JLreadBeats('RG10204', 10);
% ---50 dB: most conditions are not binaural: either I or C dominates. Crossover ~300 Hz
JLbeatPlot(JbB_50dB(2));  % 200 Hz  R = [0.37 0.55 0.041] 27 spikes; rho=0.051; I dominates
JLbeatPlot(JbB_50dB(3));  % 300 Hz  R = [0.54 0.64 0.32] 48 spikes; rho=0.15; ~balance
JLbeatPlot(JbB_50dB(4));  % 400 Hz  R = [0.41 0.81 0.3] 71 spikes; rho=0.27;  C dominates
JLbeatPlot(JbB_50dB(5));  % 500 Hz  R = [0.32 0.81 0.2] 98 spikes; rho=0.33;  C dominates
JLbeatPlot(JbB_50dB(6));  % 600 Hz  R = [0.15 0.87 0.11] 106 spikes; rho=0.38;
JLbeatPlot(JbB_50dB(7));  % 700 Hz  R = [0.26 0.87 0.27] 147 spikes; rho=0.37;
JLbeatPlot(JbB_50dB(8));  % 800 Hz  R = [0.11 0.87 0.15] 130 spikes; rho=0.34; I absent
JLbeatPlot(JbB_50dB(9));  % 900 Hz  R = [0.13 0.86 0.081] 118 spikes; rho=0.3; C dominates
JLbeatPlot(JbB_50dB(10));  % 1000 Hz  R = [0.24 0.8 0.13] 105 spikes; rho=0.25; &c
JLbeatPlot(JbB_50dB(11));  % 1100 Hz  R = [0.28 0.75 0.34] 89 spikes; rho=0.21;
JLbeatPlot(JbB_50dB(12));  % 1200 Hz  R = [0.27 0.76 0.31] 78 spikes; rho=0.18;
JLbeatPlot(JbB_50dB(13));  % 1300 Hz  R = [0.15 0.66 0.14] 55 spikes; rho=0.14;
JLbeatPlot(JbB_50dB(14));  % 1400 Hz  R = [0.36 0.64 0.22] 44 spikes; rho=0.12;
JLbeatPlot(JbB_50dB(15));  % 1500 Hz  R = [0.32 0.48 0.2] 56 spikes; rho=0.091;
% ---60 dB
ifr=2; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]); % 200 Hz  R = [0.72 0.54 0.26] 140 spikes; rho=0.25;
% 300 Hz: Ipsi peak splitting in both CHIST & waveform
ifr=3; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]); % 300 Hz  R = [0.43 0.85 0.36] 214 spikes; rho=0.41;
% 400 Hz: I peaks fuse
ifr=4; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 400 Hz  R = [0.55 0.91 0.5] 267 spikes; rho=0.48;
% 500 Hz: C dominates, but not as extremely as @ 50 dB
ifr=5; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]); % 500 Hz  R = [0.49 0.89 0.42] 272 spikes; rho=0.48;
% 600+ Hz; same story; also note strange (early) I timing of STA
ifr=6; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 600 Hz  R = [0.51 0.88 0.47] 220 spikes; rho=0.42;
ifr=7; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]); % 700 Hz  R = [0.52 0.86 0.52] 250 spikes; rho=0.39;
ifr=8; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 800 Hz  R = [0.54 0.82 0.43] 216 spikes; rho=0.36;
ifr=9; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 900 Hz  R = [0.51 0.79 0.42] 181 spikes; rho=0.34;
ifr=10; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1000 Hz  R = [0.55 0.76 0.42] 182 spikes; rho=0.31;
ifr=11; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1100 Hz  R = [0.58 0.74 0.43] 112 spikes; rho=0.27;
ifr=12; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1200 Hz  R = [0.4 0.76 0.35] 121 spikes; rho=0.23;
ifr=13; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1300 Hz  R = [0.43 0.65 0.28] 123 spikes; rho=0.19;
% now ipsi gets relatively bigger, and I STA timing returns to normal
ifr=14; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1400 Hz  R = [0.57 0.67 0.39] 127 spikes; rho=0.16;
ifr=15; JLbeatPlot([JbB_60dB(ifr) JbB2_60dB(ifr) ]);  % 1500 Hz  R = [0.3 0.7 0.19] 79 spikes; rho=0.13;
% consistency tests
ifr=2; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % okay, but ampl varies
ifr=5; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % ipsi timing changes
ifr=8; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % okay
ifr=10; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % okay, but I amp changes
ifr=12; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % okay
ifr=15; JLbeatPlot({JbB_60dB(ifr) JbB2_60dB(ifr) JbB2_60dB(ifr) }, 3); % amp chages, but stable shapes.
% ---70 dB (first round, 100:100:1500 Hz) decent, but not brilliant binaurality. Single peaks above 400 Hz. 
%  C dominates throughout.
JLbeatPlot(JbB_70dB(2)); % 200 Hz  R = [0.75 0.36 0.15] 435 spikes; rho=0.44; many spikes! long-tailed I CHIST
JLbeatPlot(JbB_70dB(3));  % 300 Hz  R = [0.58 0.83 0.41] 220 spikes; rho=0.42; same picture
JLbeatPlot(JbB_70dB(4));  % 400 Hz  R = [0.56 0.87 0.45] 153 spikes; rho=0.41;
JLbeatPlot(JbB_70dB(5));  % 500 Hz  R = [0.77 0.86 0.64] 173 spikes; rho=0.42;
JLbeatPlot(JbB_70dB(6));  % 600 Hz  R = [0.66 0.81 0.5] 187 spikes; rho=0.38;
JLbeatPlot(JbB_70dB(7));  % 700 Hz  R = [0.7 0.81 0.59] 189 spikes; rho=0.35;
JLbeatPlot(JbB_70dB(8));  % 800 Hz  R = [0.75 0.78 0.59] 162 spikes; rho=0.33;
JLbeatPlot(JbB_70dB(9));  % 900 Hz  R = [0.7 0.78 0.6] 192 spikes; rho=0.28;
JLbeatPlot(JbB_70dB(10));  % 1000 Hz  R = [0.76 0.73 0.53] 117 spikes; rho=0.27;
JLbeatPlot(JbB_70dB(11));  % 1100 Hz  R = [0.53 0.76 0.43] 90 spikes; rho=0.23;
JLbeatPlot(JbB_70dB(12));  % 1200 Hz  R = [0.69 0.63 0.45] 91 spikes; rho=0.19;
JLbeatPlot(JbB_70dB(13));  % 1300 Hz  R = [0.59 0.65 0.34] 95 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(14));  % 1400 Hz  R = [0.42 0.73 0.34] 80 spikes; rho=0.13;
JLbeatPlot(JbB_70dB(15));  % 1500 Hz  R = [0.63 0.58 0.43] 69 spikes; rho=0.11;
% ---70 dB (second round, 300:100:700 Hz) 
% 300 Hz  R = [0.53 0.84 0.39] 1305 spikes; rho=0.41; small secondary peaks. Note similarity bwt CHIST & waveforms.
ifr=1; JLbeatPlot([JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ]);
% 400 Hz  R = [0.63 0.89 0.54] 778 spikes; rho=0.4; secondary peaks reduced
ifr=2; JLbeatPlot([JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ]); 
% 500 Hz  R = [0.72 0.86 0.61] 856 spikes; rho=0.45; secondary peaks gone 
ifr=3; JLbeatPlot([JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ]);
% 600 Hz  R = [0.7 0.78 0.54] 766 spikes; rho=0.39;
ifr=4; JLbeatPlot([JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ]);
% 700 Hz  R = [0.72 0.79 0.56] 897 spikes; rho=0.37;
ifr=5; JLbeatPlot([JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ]);
%========
% SPL effects. Interesting dance of peaks @ low freqs. C is relatively stable at all freqs.
ifr=2; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=3; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3); % timing shift
ifr=4; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3); % peak competition?
ifr=5; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=6; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3); % ipsi kicking in
ifr=7; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=8; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=9; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=10; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=11; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=12; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=13; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=14; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
ifr=15; JLbeatPlot({JbB_50dB(ifr) [JbB_60dB(ifr) JbB2_60dB(ifr) ] JbB_70dB(ifr)}, 3);
%=====
% changes over time
% strange shift of I timing, while C timing is stable. Is one input "down-regulated?"
% 200 Hz
ifr=2; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi pk shift +0.44 ms
ifr=2; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); % STA ~stable
ifr=2; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6); % Ipsi peak shift 0.074 cycle = 0.37 ms
ifr=2; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 7); % nothing wrong w AP detection..
% 300 Hz
ifr=3; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi peak timing changes by +1.05 ms
ifr=3; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); % STA is ~stable
ifr=3; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6);%timing re stim changes by 0.32 cycle=1.05 ms
ifr=3; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 7); % nothing wrong w AP detection..
% 400 Hz
ifr=4; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); %Ipsi pk shift +0.6 ms
ifr=4; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); % STA is ~stable
ifr=4; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6);%timing re stim changes by 0.22 cycle ~ 0.54 ms
% 500 Hz
ifr=5; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi pk shift 0.52 ms
ifr=5; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); % STA is stable
ifr=5; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6); % Ipsi phase shift 0.21 cycle = 0.4 ms
% 600 Hz
ifr=6; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi pk shift 0.12 ms
ifr=6; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); %STA ~ stable (but Ipsi contrib is weak)
ifr=6; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6); % Ipsi phase shift 0.1 cycle = 0.17 ms
% 700 Hz 
ifr=7; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi pk shift 0.18 ms
ifr=7; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); %STA ~ stable? Hard to tell.
ifr=7; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6); % Ipsi phase shift 0.15 cycle = 0.21 ms
% 800 Hz
ifr=8; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 3); % Ipsi pk shift 0.044 ms
ifr=8; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 5); %STA ~ stable? Hard to tell.
ifr=8; JLbeatPlot({ JbB_60dB(ifr)  JbB_60dB(ifr) JbB2_60dB(ifr)}, 6); % Ipsi phase shift 0.054 cycle = 0.068 ms

% ====linearity
% 50 dB 
ifr = 2; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % too noisy; I driven
ifr = 3; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % pretty linear
ifr = 4; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear
ifr = 5; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; C dominates
ifr = 6; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; C dominates even more
ifr = 7; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; strong C dominance
ifr = 8; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ipsi gone
ifr = 10; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; tiny, consistent,  I contrib
ifr = 12; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ditto
ifr = 15; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ditto
% 60 dB
ifr = 2; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); %linear. 
ifr = 3; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % strange Ipsi hickup
ifr = 4; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); %mutual suppr; stable shapes 
ifr = 5; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); %mutual suppr; stable shapes 
ifr = 6; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % I->C suppr strongest
ifr = 7; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % same
ifr = 8; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % strong I->C suppr
ifr = 9; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % same
ifr = 10; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3);% only strong I->C suppr
ifr = 11; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3);% same
ifr = 12; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3);% same
ifr = 13; JLbeatPlot({[JbB_60dB(ifr) JbB2_60dB(ifr) ] JbI_60dB(ifr) JbC_60dB(ifr) }, 3);% same
% 70 dB (first round)
ifr = 2; JLbeatPlot({JbB_70dB(ifr) JbI_70dB(ifr) JbC_70dB(ifr) }, 3); 


% XXXXXXXX****** compare monaural & binaural
ifr = 1; % 300 Hz  504 + 587 = 1305 spikes. I timing is degraded (0.89->0.53) by adding contra stim.
% .. is that degradation the main cuase for the small I amplitude in the
% STA? The I-CHIST gets long tails, as ifC alone as well able to evoke
% spikes, irregardless of binaural phase. Most of htese measurements show
% good linearlity in their waveforms. Does this help to undersatnd the link
% between the waveforms and the spikes? What is the prediction of vector strength and
% CHIST shapes if I & C combine in an OR way? The real bin CHIST show
% "negative" parts. Is this inhibition? Can the nonlinear way in which 
% the CHIST combine be explained in terms of the way in which the waveforms
% combine?
II7 = [JbI_71dB(ifr) JbI2_71dB(ifr) JbI3_71dB(ifr) JbI4_71dB(ifr) JbI5_71dB(ifr) JbI6_71dB(ifr) ];
BB7 = [JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ];
CC7 = [JbC_71dB(ifr) JbC2_71dB(ifr) JbC3_71dB(ifr) JbC4_71dB(ifr) JbC5_71dB(ifr) JbC6_71dB(ifr) ];
JLbeatPlot({BB7 CC7 II7}, 3); JLbeatPlot({BB7 CC7 II7}, 5); JLbeatPlot({BB7 CC7 II7}, 6);
ifr = 2; % 400 Hz  292+441=778 spikes. I timing 0.83->0.63; C timing 0.94->0.89. C dominates STA.
II7 = [JbI_71dB(ifr) JbI2_71dB(ifr) JbI3_71dB(ifr) JbI4_71dB(ifr) JbI5_71dB(ifr) JbI6_71dB(ifr) ];
BB7 = [JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ];
CC7 = [JbC_71dB(ifr) JbC2_71dB(ifr) JbC3_71dB(ifr) JbC4_71dB(ifr) JbC5_71dB(ifr) JbC6_71dB(ifr) ];
JLbeatPlot({BB7 CC7 II7}, 3); JLbeatPlot({BB7 CC7 II7}, 5); JLbeatPlot({BB7 CC7 II7}, 6);
ifr = 3; % 500 Hz  376+334=856 spikes. I timing 0.86->0.72; C timing 0.91->0.86. C&I almost equivalent.
II7 = [JbI_71dB(ifr) JbI2_71dB(ifr) JbI3_71dB(ifr) JbI4_71dB(ifr) JbI5_71dB(ifr) JbI6_71dB(ifr) ];
BB7 = [JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ];
CC7 = [JbC_71dB(ifr) JbC2_71dB(ifr) JbC3_71dB(ifr) JbC4_71dB(ifr) JbC5_71dB(ifr) JbC6_71dB(ifr) ];
JLbeatPlot({BB7 CC7 II7}, 3); JLbeatPlot({BB7 CC7 II7}, 5); JLbeatPlot({BB7 CC7 II7}, 6);
ifr = 4; % 600 Hz  368+340=766 spikes. I timing  0.76->0.70; C timing 0.83-> 0.78. C dominant.
II7 = [JbI_71dB(ifr) JbI2_71dB(ifr) JbI3_71dB(ifr) JbI4_71dB(ifr) JbI5_71dB(ifr) JbI6_71dB(ifr) ];
BB7 = [JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ];
CC7 = [JbC_71dB(ifr) JbC2_71dB(ifr) JbC3_71dB(ifr) JbC4_71dB(ifr) JbC5_71dB(ifr) JbC6_71dB(ifr) ];
JLbeatPlot({BB7 CC7 II7}, 3); JLbeatPlot({BB7 CC7 II7}, 5); JLbeatPlot({BB7 CC7 II7}, 6);
ifr = 5; % 700 Hz  526+354=897 spikes; I timing 0.84->0.72; C timing 0.86->0.79. C&I almost equivalent.
II7 = [JbI_71dB(ifr) JbI2_71dB(ifr) JbI3_71dB(ifr) JbI4_71dB(ifr) JbI5_71dB(ifr) JbI6_71dB(ifr) ];
BB7 = [JbB_71dB(ifr) JbB2_71dB(ifr) JbB3_71dB(ifr) JbB4_71dB(ifr) JbB5_71dB(ifr) JbB6_71dB(ifr) ];
CC7 = [JbC_71dB(ifr) JbC2_71dB(ifr) JbC3_71dB(ifr) JbC4_71dB(ifr) JbC5_71dB(ifr) JbC6_71dB(ifr) ];
JLbeatPlot({BB7 CC7 II7}, 3); JLbeatPlot({BB7 CC7 II7}, 5); JLbeatPlot({BB7 CC7 II7}, 6);
%
% % Previous Observations from this cell:
% %    - spike-trig av look like the loop-averaged traces. Makes sense
% %      because the spikes are phase locked to both ears
% %    - 70 dB, 400 Hz: from the periodic analys (plot 3) there is little
% %      interaction: ipsi & contra are not different in binaural condition
% %      and monaural conditions. BUT when looking at spike-trig av (plot 5),
% %      the phase locking to the ipsi ear is suppressed by the presence of
% %      the contra ear. This is also visible in the cycle histograms (plot
% %      6), where the ipsi-R drops from 0.83 (monaural) to 0.63 (binaural),
% %      while the contra-R is hardly affected (0.94 -> 0.89). It would 
% %      appear as though the occurence of spikes affects the way the ipsi 
% %      ear (but not the contra) triggers the other spikes. But a much
% %      simpler explanation is that under binaural stimulation, the contra 
% %      ear triggers spikes of its own, without the need for coincident
% %      input from ipsi. This dilutes ipsi phase locking and thereby reduces
% %      the ipsi-STA.
% ifreq=2; % 400 Hz
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},6); % CHST
% %
% %      At 300 Hz, the situation is comparable (ipsi suffers from contra):
% ifreq=1; % 300 Hz
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},6); % CHST
% %
% %      At 200 Hz, the effect is reversed: now phase locking to contra 
% %      suffers more from ipsi stim than vice versa:
% ifreq=2; % 200 Hz (no pooled data)
% pp=0; JLbeatPlot({JbB_70dB(ifreq) JbI_70dB(ifreq) JbC_70dB(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_70dB(ifreq) JbI_70dB(ifreq) JbC_70dB(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_70dB(ifreq) JbI_70dB(ifreq) JbC_70dB(ifreq)},6); % CHST
% %      The "refractory period" in the binaural CHSTs (both ipsi & contra)
% %      shows that the dips in the ipsi&contra period-analysis are true. It
% %      is as if an EPSP has a negative after effect (or is immediately 
% %      followed by an IPSP). Note that the "suppression" is now restricted
% %      to spike-related metrics; the averaged waveforms don't show it.
% %      
% %      In the above analyses, note that the shape of the CHST is highly
% %      correlated with the monaural periodic means. Together with the fact
% %      that the latter precedes the former by ~200 us, this strongly 
% %      suggests that the periodic means reflect the potential that
% %      modulates spike initiation. The correlation between the two should
% %      be anayzed quantitatively.
% %      same analysis on 60-dB data:
% ifreq=2; % 200 Hz (no pooled data)
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},6); % CHST
% %      same analysis on 50-dB data:
% ifreq=2; % 200 Hz (no pooled data)
% pp=0; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq)},6); % CHST
% %   (contra ear does not do much)
% %
% %      At higher frequencies, the across-ear effect of stim is less:
% ifreq=3; % 500 Hz
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},3); % periodic analysis
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},5); % STA 
% pp=0; JLbeatPlot({JbB_70dB_p(ifreq) JbI_70dB_p(ifreq) JbC_70dB_p(ifreq)},6); % CHST
% %      Note (again) the good correlation between the shapes of the monaural
% %      CHSTs and the period-averaged potential.
% 
% % >>>>>>>>>>>>cell 12
% clear Jb*
% APslope = [-5 -0.5 1]; Nav = 3;
% %  60 dB SPL
% for ifreq=1:15, % ifreq*100 Hz
%     JbB_60dB(ifreq) = JLbeat('RG10204', '12-1-',ifreq, Nav, APslope); % binaural
%     JbI_60dB(ifreq) = JLbeat('RG10204', '12-2-',ifreq, Nav, APslope); % ipsi
%     JbC_60dB(ifreq) = JLbeat('RG10204', '12-3-',ifreq, Nav, APslope); % contra
% end
% ifreq = 1; pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},3); % 100 Hz
% ifreq = 2; pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},3); % 200 Hz
% 
% % Nice:
% ifreq = 3; % 300 Hz
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},3); 
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},5); 
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},6); 
% pp=0; JLbeatPlot(JbB_60dB(ifreq)); 
% ifreq = 5; % 500 Hz
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},3); 
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},5); 
% pp=0; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq)},6); 
% pp=0; JLbeatPlot(JbB_60dB(ifreq)); 
% % IPSI-only evokes no spikes..
% 
% % real-time comparison of waveforms
% JLwaterfall(JbB_70dB_p);
qq=[]
% RG10204 >>>>> cell 11 ("EPSPs + very small spikes")  
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
% general: spike detection is problematic. E.g., JbB_50dB(5) average waveforms suggest that APs were missed.
%  ==60 dB
% secondary peaks (if any) disappear @ 300 Hz and up
JLbeatPlot(JbB_60dB(1)); % not responsive
JLbeatPlot(JbB_60dB(2)); % 200 Hz  R = [0 0 0] 0 spikes; rho=0.23; % C ~ multi peaks spacing ~2 ms
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0.98 0.98 0.99] 4 spikes; rho=0.41; 4 perfectly timed spikes..(not just onset)
JLbeatPlot(JbB_60dB(4)); % 400 Hz  R = [0 0 0] 1 spikes; rho=0.49; 
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.94 0.96 0.86] 5 spikes; rho=0.56; 5 perfectly timed spikes..
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.97 0.91 0.98] 2 spikes; rho=0.49;
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0 0 0] 1 spikes; rho=0.45;
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0 0 0] 1 spikes; rho=0.33;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0 0 0] 1 spikes; rho=0.27;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.18;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0 0 0] 1 spikes; rho=0.13;
%  ==50 dB 
JLbeatPlot(JbB_50dB(1));  % 100 Hz  R = [0 0 0] 1 spikes; rho=-0.0033; not responsive
JLbeatPlot(JbB_50dB(2));  % 200 Hz  R = [0.98 0.64 0.55] 3 spikes; rho=0.076; still very noise
JLbeatPlot(JbB_50dB(3));  % 300 Hz  R = [0.91 0.87 0.86] 3 spikes; rho=0.26;
JLbeatPlot(JbB_50dB(4));  % 400 Hz  R = [0.9 0.5 0.47] 11 spikes; rho=0.41; I dominates
JLbeatPlot(JbB_50dB(5));  % 500 Hz  R = [0.93 0.78 0.66] 9 spikes; rho=0.38;
JLbeatPlot(JbB_50dB(6));  % 600 Hz  R = [0.75 0.92 0.89] 5 spikes; rho=0.27; % 5 well-timed APs
JLbeatPlot(JbB_50dB(7));  % 700 Hz  R = [0.8 0.49 0.17] 6 spikes; rho=0.31;
JLbeatPlot(JbB_50dB(8));  % 800 Hz  R = [0.78 0.93 0.94] 3 spikes; rho=0.2;
JLbeatPlot(JbB_50dB(9));  % 900 Hz  R = [0.81 0.65 0.59] 11 spikes; rho=0.17;
JLbeatPlot(JbB_50dB(10));  % 1000 Hz  R = [0.63 0.97 0.8] 2 spikes; rho=0.13;
JLbeatPlot(JbB_50dB(11));  % 1100 Hz  R = [0.96 0.8 0.6] 2 spikes; rho=0.07;
%  ==40 dB 
JLbeatPlot(JbB_40dB(3));  % 300 Hz  R = [0 0 0] 1 spikes; rho=0.052;  I dominates waveform
JLbeatPlot(JbB_40dB(4));  % 400 Hz  R = [0.37 0.97 0.6] 2 spikes; rho=0.14;  I dominates
JLbeatPlot(JbB_40dB(5));  % 500 Hz  R = [0.81 0.86 1] 2 spikes; rho=0.094;
JLbeatPlot(JbB_40dB(6));  % 600 Hz  R = [0.56 0.81 0.035] 2 spikes; rho=0.067;  C slightlyy dominating
JLbeatPlot(JbB_40dB(7)); % 700 Hz  R = [0.71 0.55 0.66] 5 spikes; rho=0.1; noisy
% linearity. Maybe slight reduction in Ampl, but shapes are highly preserved.
% --60 dB
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % pretty linear
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 4; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3);
ifr = 5; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % Amp reduction both I&C
ifr = 6; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % Amp reduction both I&C
ifr = 7; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 8; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 12; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % same story
% --50 dB (note order of recordings: C-I-B)
ifr = 2; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % C weak; ipsi enhanced
ifr = 3; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % I dominant. Very linear.
ifr = 4; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % I dominant. Very linear
ifr = 5; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % I~C. Shapes stable; slight Amp reduction
ifr = 6; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % C starts to dominate. Stable shapes & sizes.
ifr = 7; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ditto
ifr = 12; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % all the way up ..
ifr = 15; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % .. to high freqs
% --40 dB
ifr = 2; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % very noisy
ifr = 3; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % noisy; Strong I dom. Decent linearity.
ifr = 4; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % I rules. Pretty linear.
ifr = 5; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % I/C balance. Linear.
ifr = 6; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % C slightly bigger. Linear.
ifr = 7; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % I~C. Linearity.
%
% RG10204 >>>> cell 12 ("spikes + EPSPs")  
%      30 dB bin/contra
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
JLreadBeats('RG10204', 12);
%  ==60 dB
% single-peaked from 300 Hz. High driven rates @ 500-700 Hz. False AP triggers?
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0 0 0] 0 spikes; rho=0.051; Only an I contrib
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0 0 0] 1 spikes; rho=0.16; I/C balance. I is multi-peaked
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.79 0.83 0.64] 13 spikes; rho=0.43;
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.71 0.92 0.65] 41 spikes; rho=0.51; % decent binaural detection
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.64 0.93 0.57] 200 spikes; rho=0.53; many spikes! C dominates
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.51 0.93 0.44] 219 spikes; rho=0.47; many spikes! C dominates
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.6 0.91 0.51] 181 spikes; rho=0.43; % many spikes; C dominates
JLbeatPlot(JbB_60dB(8)); % 800 Hz  R = [0.7 0.89 0.61] 127 spikes; rho=0.38; C dominates. Still good binaurality.
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.67 0.9 0.6] 77 spikes; rho=0.31; C dominates. Good binaurality.
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.71 0.91 0.67] 55 spikes; rho=0.27; C dominates. Good binaurality.
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.48 0.89 0.44] 52 spikes; rho=0.24;
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.55 0.85 0.45] 48 spikes; rho=0.21;
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0.59 0.83 0.48] 29 spikes; rho=0.18;
JLbeatPlot(JbB_60dB(14));  % 1400 Hz  R = [0.52 0.84 0.49] 21 spikes; rho=0.14; 
JLbeatPlot(JbB_60dB(15));  % 1500 Hz  R = [0.58 0.79 0.44] 16 spikes; rho=0.09; % C dominates all the time
%  ==50 dB 
%  note nonmonotonic response (more spikes @ 50 dB than @ 60 dB)
JLbeatPlot(JbB_50dB(1));  % 100 Hz  R = [0.73 0.65 0.033] 2 spikes; rho=-0.00019; % noise
JLbeatPlot(JbB_50dB(2));  % 200 Hz  R = [0.18 0.74 0.29] 4 spikes; rho=0.067; noisy. C/I balance
JLbeatPlot(JbB_50dB(3));  % 300 Hz  R = [0.14 0.73 0.19] 11 spikes; rho=0.18; I wv larger, but C dominates locking
JLbeatPlot(JbB_50dB(4));  % 400 Hz  R = [0.71 0.77 0.47] 5 spikes; rho=0.34; I dominates waveform
JLbeatPlot(JbB_50dB(5));  % 500 Hz  R = [0.65 0.85 0.52] 40 spikes; rho=0.39; I/C balance.
JLbeatPlot(JbB_50dB(6));  % 600 Hz  R = [0.39 0.92 0.39] 104 spikes; rho=0.42; strong C dominance
JLbeatPlot(JbB_50dB(7));  % 700 Hz  R = [0.4 0.93 0.35] 253 spikes; rho=0.45; ditto
JLbeatPlot(JbB_50dB(8));  % 800 Hz  R = [0.51 0.91 0.47] 275 spikes; rho=0.39; C still dominating
JLbeatPlot(JbB_50dB(9));  % 900 Hz  R = [0.54 0.9 0.48] 210 spikes; rho=0.36;
JLbeatPlot(JbB_50dB(10));  % 1000 Hz  R = [0.5 0.88 0.42] 139 spikes; rho=0.32;
JLbeatPlot(JbB_50dB(11));  % 1100 Hz  R = [0.44 0.84 0.35] 78 spikes; rho=0.28;
JLbeatPlot(JbB_50dB(12));  % 1200 Hz  R = [0.42 0.85 0.38] 94 spikes; rho=0.23;  still string C dominance
JLbeatPlot(JbB_50dB(13));  % 1300 Hz  R = [0.47 0.77 0.35] 62 spikes; rho=0.18;
JLbeatPlot(JbB_50dB(14));  % 1400 Hz  R = [0.3 0.78 0.21] 45 spikes; rho=0.14;
JLbeatPlot(JbB_50dB(15));  % 1500 Hz  R = [0.37 0.71 0.31] 41 spikes; rho=0.081;
%  ==40 dB 
JLbeatPlot(JbB_40dB(3));  % 300 Hz  R = [0.23 0.47 0.5] 4 spikes; rho=0.0082; noisy. C/I balance
JLbeatPlot(JbB_40dB(4));  % 400 Hz  R = [0.86 0.92 0.99] 2 spikes; rho=0.056; I dominates. 2 well-timed spikes..
JLbeatPlot(JbB_40dB(5));  % 500 Hz  R = [0.61 0.77 0.58] 7 spikes; rho=0.12; slight C dominance
JLbeatPlot(JbB_40dB(6));  % 600 Hz  R = [0.33 0.85 0.2] 13 spikes; rho=0.11; strong C dominance. Not really binaural.
JLbeatPlot(JbB_40dB(7));  % 700 Hz  R = [0.29 0.8 0.26] 34 spikes; rho=0.2;
JLbeatPlot(JbB_40dB(8));  % 800 Hz  R = [0.35 0.9 0.29] 153 spikes; rho=0.28; % high rate, but mainly C driven
JLbeatPlot(JbB_40dB(9));  % 900 Hz  R = [0.29 0.84 0.19] 114 spikes; rho=0.23; ditto
JLbeatPlot(JbB_40dB(10));  % 1000 Hz  R = [0.16 0.81 0.043] 90 spikes; rho=0.19; ditto
JLbeatPlot(JbB_40dB(11));  % 1100 Hz  R = [0.15 0.74 0.1] 71 spikes; rho=0.13;
JLbeatPlot(JbB_40dB(12));  % 1200 Hz  R = [0.27 0.72 0.21] 47 spikes; rho=0.086;
JLbeatPlot(JbB_40dB(13));  % 1300 Hz  R = [0.19 0.75 0.19] 41 spikes; rho=0.064;
JLbeatPlot(JbB_40dB(14));  % 1400 Hz  R = [0.13 0.5 0.13] 39 spikes; rho=0.041;
JLbeatPlot(JbB_40dB(15));  % 1500 Hz  R = [0.25 0.56 0.29] 18 spikes; rho=0.023;
%  ==30 dB 
% well driven by C, also decent phase locking. Not truly bin.
JLbeatPlot(JbB_30dB(3));  % 300 Hz  R = [0.4 0.18 0.21] 9 spikes; rho=-0.0086;  noisy
JLbeatPlot(JbB_30dB(4));  % 400 Hz  R = [0.29 0.048 0.47] 7 spikes; rho=-0.0022;
JLbeatPlot(JbB_30dB(5));  % 500 Hz  R = [0.56 0.46 0.37] 9 spikes; rho=-0.00093; somewhat binaural
JLbeatPlot(JbB_30dB(6));  % 600 Hz  R = [0.48 0.2 0.3] 7 spikes; rho=-0.003;
JLbeatPlot(JbB_30dB(7));  % 700 Hz  R = [0.25 0.38 0.22] 7 spikes; rho=0.012; C dominance
JLbeatPlot(JbB_30dB(8));  % 800 Hz  R = [0.43 0.69 0.27] 17 spikes; rho=0.048;
JLbeatPlot(JbB_30dB(9));  % 900 Hz  R = [0.33 0.82 0.25] 36 spikes; rho=0.038;
JLbeatPlot(JbB_30dB(10));  % 1000 Hz  R = [0.077 0.8 0.047] 45 spikes; rho=0.021;
JLbeatPlot(JbB_30dB(11));  % 1100 Hz  R = [0.15 0.61 0.065] 20 spikes; rho=0.014;
JLbeatPlot(JbB_30dB(12));  % 1200 Hz  R = [0.015 0.5 0.18] 11 spikes; rho=0.0057;
JLbeatPlot(JbB_30dB(13));  % 1300 Hz  R = [0.066 0.53 0.19] 16 spikes; rho=-0.0058;
JLbeatPlot(JbB_30dB(14));  % 1400 Hz  R = [0.29 0.49 0.044] 20 spikes; rho=-0.0059;
JLbeatPlot(JbB_30dB(15));  % 1500 Hz  R = [0.37 0.57 0.15] 18 spikes; rho=-0.0042;
% linearity. 
% --60 dB
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % noisy, but pretty linear
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % stable shapes; may slight Ampl reduction
ifr = 4; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 5; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 6; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 7; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 8; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ipsi somewhat suppressed
ifr = 9; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % slight interaural effects
ifr = 10; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % pretty linear
ifr = 11; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 12; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear
ifr = 13; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3);
% --50 dB
ifr = 2; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % noisy, but pretty linear
ifr = 3; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction; C enhancement
ifr = 4; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction; C enhancement
ifr = 5; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % same?
ifr = 6; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~linear
ifr = 7; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction.
ifr = 8; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction.
ifr = 9; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction.
ifr = 10; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes. I reduction.
ifr = 11; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~linear
ifr = 12; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~linear
ifr = 13; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3);  % ~linear
ifr = 14; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % slight C enhancement
ifr = 15; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~linear
% --40 dB
ifr = 2; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % too noisy
ifr = 3; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % ~lin
ifr = 4; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % ~lin. slight C enhancement
ifr = 5; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % ~lin
ifr = 7; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % ~lin
ifr = 10; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % lin
ifr = 12; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3);
ifr = 15; JLbeatPlot({JbB_40dB(ifr) JbI_40dB(ifr) JbC_40dB(ifr) }, 3); % lin

% RG10204 >>> cell 13 ("spikes + EPSPs", "strong beat")  
%      40 dB bin
%      50 dB bin/ipsi/contra
%      60 dB bin/ipsi/contra
% NTD
% ==60 dB  all mono peaked
JLreadBeats('RG10204', 13);
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0 0 0] 0 spikes; rho=0.082; I only
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.16 0.89 0.32] 4 spikes; rho=0.26; C kicks in
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.34 0.96 0.29] 13 spikes; rho=0.41;
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.37 0.95 0.38] 36 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.51 0.94 0.46] 213 spikes; rho=0.52;  C dominates wv and locking
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.42 0.92 0.39] 269 spikes; rho=0.51;
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.55 0.92 0.51] 242 spikes; rho=0.49;
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0.53 0.9 0.46] 210 spikes; rho=0.43;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.4 0.88 0.34] 149 spikes; rho=0.41;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.43 0.89 0.38] 111 spikes; rho=0.38;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.45 0.86 0.41] 92 spikes; rho=0.36;
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.3 0.89 0.26] 61 spikes; rho=0.32;
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0.15 0.89 0.23] 44 spikes; rho=0.31;
% ==50 dB  note nonmonotonicity ~1000 Hz (more spikes @50 dB than @60 dB)
JLbeatPlot(JbB_50dB(1));  % 100 Hz  R = [0.35 0.44 0.51] 4 spikes; rho=0.003; % too noisy
JLbeatPlot(JbB_50dB(2));  % 200 Hz  R = [0.15 0.69 0.036] 9 spikes; rho=0.081; I dom waveform, but C, locking
JLbeatPlot(JbB_50dB(3));  % 300 Hz  R = [0.35 0.53 0.48] 10 spikes; rho=0.16;
JLbeatPlot(JbB_50dB(4));  % 400 Hz  R = [0.31 0.55 0.41] 9 spikes; rho=0.12; % I-C balance
JLbeatPlot(JbB_50dB(5));  % 500 Hz  R = [0.36 0.82 0.45] 35 spikes; rho=0.28; I~C in wave, but C dominates locking
JLbeatPlot(JbB_50dB(6));  % 600 Hz  R = [0.37 0.89 0.27] 46 spikes; rho=0.27; C dominates wave & locking
JLbeatPlot(JbB_50dB(7));  % 700 Hz  R = [0.3 0.89 0.27] 232 spikes; rho=0.42; C dominates
JLbeatPlot(JbB_50dB(8));  % 800 Hz  R = [0.25 0.9 0.24] 532 spikes; rho=0.37;
JLbeatPlot(JbB_50dB(9));  % 900 Hz  R = [0.22 0.89 0.2] 437 spikes; rho=0.36; C very dominant
JLbeatPlot(JbB_50dB(10));  % 1000 Hz  R = [0.21 0.86 0.18] 268 spikes; rho=0.38;
JLbeatPlot(JbB_50dB(12));  % 1200 Hz  R = [0.21 0.85 0.16] 200 spikes; rho=0.32;
JLbeatPlot(JbB_50dB(15));  % 1500 Hz  R = [0.26 0.77 0.17] 104 spikes; rho=0.21;
% ==40 dB  I is too weak - hardly binaural.
JLbeatPlot(JbB_40dB(3));  % 300 Hz  R = [0.38 0.041 0.42] 10 spikes; rho=0.015; too noisy
JLbeatPlot(JbB_40dB(4));  % 400 Hz  R = [0.23 0.36 0.17] 13 spikes; rho=0.0083;
JLbeatPlot(JbB_40dB(5));  % 500 Hz  R = [0.081 0.42 0.25] 15 spikes; rho=0.05;
JLbeatPlot(JbB_40dB(6));  % 600 Hz  R = [0.38 0.64 0.032] 24 spikes; rho=0.038;
JLbeatPlot(JbB_40dB(7));  % 700 Hz  R = [0.18 0.75 0.076] 32 spikes; rho=0.093;
JLbeatPlot(JbB_40dB(8));  % 800 Hz  R = [0.31 0.84 0.25] 85 spikes; rho=0.23;
JLbeatPlot(JbB_40dB(9));  % 900 Hz  R = [0.17 0.87 0.14] 94 spikes; rho=0.22;
JLbeatPlot(JbB_40dB(10));  % 1000 Hz  R = [0.14 0.82 0.13] 112 spikes; rho=0.2;
JLbeatPlot(JbB_40dB(11));  % 1100 Hz  R = [0.1 0.79 0.086] 72 spikes; rho=0.15;
JLbeatPlot(JbB_40dB(12));  % 1200 Hz  R = [0.086 0.78 0.054] 59 spikes; rho=0.13;
% facilitation vs inhibition: spike counts C-I-B
% inhib, 60 dB
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5);  % 42+1=36 spikes; Ee or E0
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5); % 231+2=213 spikes; Ee or E0
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5);  % 262+0=269; Ee or E0
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5); % 210+2=242; Ee cell
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5); % 150+0=210; E0 cell ?!
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5);  % 116+0=149 spikes; E0
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5);  % 85+3=111 spikes; E0
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5); % 83+1=92 spikes; E0
ifreq = 12; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 5); % 79+0=61 spikes; E0
% inhib, 50 dB
ifreq = 7; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 221+7=232 spikes; Ee or E0
ifreq = 8; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 573+7=532 spikes; EI
ifreq = 9; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5);% 452+7=437 spikes ~EI
ifreq = 10; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 344+4 = 268 spikes EI
ifreq = 11; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 234+11=282 spikes; Ee or E0
ifreq = 12; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5);  % 199+5=200 spikes; Ee or E0
% linearity. 
% --60 dB
ifr = 2; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % noise but linear
ifr = 3; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % very linear; maybe slight C suppr
ifr = 4; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % linear; C suppr
ifr = 5; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % stable shapes; mutual suppr
ifr = 6; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ~lin; mutual supp
ifr = 7; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % stable shape; mutual suppr
ifr = 8; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 9; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 10; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
ifr = 11; JLbeatPlot({JbB_60dB(ifr) JbI_60dB(ifr) JbC_60dB(ifr) }, 3); % ditto
% --50 dB
ifr = 3; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % noisy but linear
ifr = 4; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~lin; some mutual suppr
ifr = 5; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % stable shapes; I is suppressed
ifr = 6; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~lin
ifr = 7; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; some I supp
ifr = 8; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % linear; some I supp
ifr = 9; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ~lin; I suppressed; C enhanced
ifr = 10; JLbeatPlot({JbB_50dB(ifr) JbI_50dB(ifr) JbC_50dB(ifr) }, 3); % ditto

function RG10209
% RG10209 >>> cell 2 ("spikes + EPSPs")  
%      60 dB bin/ipsi/contra
% -3   2-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -4   2-3-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -5   2-4-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
JLreadBeats('RG10209', 2);
% huge spikes. Cell is binaural , but contra is dominant.
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.29 0.76 0.2] 41 spikes; rho=0.27; C driven
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.37 0.93 0.29] 95 spikes; rho=0.46; C very dominant
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.41 0.94 0.38] 145 spikes; rho=0.5; 
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.46 0.94 0.42] 167 spikes; rho=0.51; slowly getting more binaural ..
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.48 0.93 0.42] 216 spikes; rho=0.51;
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.51 0.88 0.46] 164 spikes; rho=0.46;
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0.54 0.86 0.46] 201 spikes; rho=0.41;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.55 0.87 0.47] 193 spikes; rho=0.41;
JLbeatPlot(JbB_60dB(10)); % 1000 Hz  R = [0.46 0.87 0.36] 168 spikes; rho=0.42;
JLbeatPlot(JbB_60dB(11)); % 1100 Hz  R = [0.47 0.86 0.41] 148 spikes; rho=0.4;
JLbeatPlot(JbB_60dB(12)); % 1200 Hz  R = [0.4 0.82 0.34] 119 spikes; rho=0.4;
JLbeatPlot(JbB_60dB(13)); % 1300 Hz  R = [0.35 0.78 0.27] 99 spikes; rho=0.39;
JLbeatPlot(JbB_60dB(14));  % 1400 Hz  R = [0.3 0.78 0.27] 83 spikes; rho=0.39;
JLbeatPlot(JbB_60dB(15));  % 1500 Hz  R = [0.29 0.73 0.28] 57 spikes; rho=0.36;  % C dominates waveform
% linearity test
ifreq = 2; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % linear w some mutual suppr
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % linear w some mutual suppr
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % "clipping" of I & C waveforms
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % "clipping" of I & C waveforms
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % same clipping
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suprr; waveforms ~stable
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suprr; waveforms ~stable
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % ditto
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % ditto
ifreq = 12; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % ditto
ifreq = 15; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % now it's more linear

% RG10209 cell 3 ("spikes + EPSPs" "isolation changed during recording")  
%      60 dB 3xbin/ipsi/contra
% XXX-6   3-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B artifacts 
% -7   3-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B isolation changes during recording, better use 3-3-
% -8   3-3-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B signal still changing
% -9   3-4-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B signal more stable, but smaller
% -10  3-5-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
%
JLreadBeats('RG10209', 3);
% stability problems: both AP size and # APs changes drastically over recordings!
% 200 Hz: double peaked C locking -> vector strength underestimates locking
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.91 0.53 0.51] 10 spikes; rho=0.4;
JLbeatPlot(JbB2_60dB(2));  % 200 Hz  R = [0.86 0.32 0.15] 27 spikes; rho=0.36;
ifreq = 2; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % shrinking but consistent
% 300 Hz
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.94 0.57 0.5] 37 spikes; rho=0.44;
JLbeatPlot(JbB2_60dB(3));  % 300 Hz  R = [0.93 0.54 0.43] 49 spikes; rho=0.39;
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % shrinking but consistent
% 400 Hz
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.94 0.6 0.51] 182 spikes; rho=0.53;
JLbeatPlot(JbB2_60dB(4));  % 400 Hz  R = [0.95 0.52 0.45] 240 spikes; rho=0.42;
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % shrinking but consistent
% 500 Hz
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.95 0.7 0.63] 308 spikes; rho=0.49;
JLbeatPlot(JbB2_60dB(5)); % 500 Hz  R = [0.95 0.71 0.65] 208 spikes; rho=0.39;
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % waveform changed (ipsi a lot)
% 600 Hz
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.94 0.81 0.72] 96 spikes; rho=0.36;
JLbeatPlot(JbB2_60dB(6));  % 600 Hz  R = [0.93 0.77 0.68] 241 spikes; rho=0.32;
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % ipsi changes
% 700 Hz
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.9 0.86 0.78] 24 spikes; rho=0.29;
JLbeatPlot(JbB2_60dB(7));  % 700 Hz  R = [0.76 0.85 0.6] 84 spikes; rho=0.28;
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % shrinking (maybe asymmetric)
% 800 Hz
JLbeatPlot(JbB_60dB(8)); % 800 Hz  R = [0.93 0.78 0.72] 51 spikes; rho=0.28;
JLbeatPlot(JbB2_60dB(8));  % 800 Hz  R = [0.88 0.78 0.71] 167 spikes; rho=0.24;
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 900 Hz  
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.92 0.8 0.76] 77 spikes; rho=0.33; % impressive binaurality @ 900 Hz
JLbeatPlot(JbB2_60dB(9));  % 900 Hz  R = [0.9 0.79 0.71] 119 spikes; rho=0.23;
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 1000 Hz  
JLbeatPlot(JbB_60dB(10));  
JLbeatPlot(JbB2_60dB(10)); 
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 1000 Hz  
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.87 0.76 0.72] 121 spikes; rho=0.38; 
JLbeatPlot(JbB2_60dB(10));  % 1000 Hz  R = [0.86 0.67 0.61] 98 spikes; rho=0.24;
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % asymmetric distortion
% 1100 Hz  
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.84 0.75 0.7] 50 spikes; rho=0.33;
JLbeatPlot(JbB2_60dB(11));  % 1100 Hz  R = [0.84 0.73 0.66] 131 spikes; rho=0.28;
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 1200 Hz  
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.86 0.74 0.68] 46 spikes; rho=0.31;
JLbeatPlot(JbB2_60dB(12)); % 1200 Hz  R = [0.77 0.7 0.59] 99 spikes; rho=0.38;
ifreq = 12; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 1300 Hz  
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0.87 0.63 0.57] 14 spikes; rho=0.27; false triggers
JLbeatPlot(JbB2_60dB(13));  % 1300 Hz  R = [0.68 0.48 0.49] 64 spikes; rho=0.29; false triggers
ifreq = 13; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % little change
% 1400 Hz  
JLbeatPlot(JbB_60dB(14));  % 1400 Hz  R = [0.86 0.72 0.69] 16 spikes; rho=0.24; FALSE TRIGGERS
JLbeatPlot(JbB2_60dB(14)); % 1400 Hz  R = [0.79 0.79 0.5] 15 spikes; rho=0.3;
ifreq = 14; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
% 1500 Hz  
JLbeatPlot(JbB_60dB(15));  % 1500 Hz  R = [0.76 0.68 0.73] 19 spikes; rho=0.17;
JLbeatPlot(JbB2_60dB(15));  % 1500 Hz  R = [0.75 0.67 0.51] 14 spikes; rho=0.21;
ifreq = 15; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % mainly shrinking 
%----
% linearity tests:  % difficult to assess w all this shrinking
ifreq = 2; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % linear
ifreq = 2; JLbeatPlot({JbB2_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);
ifreq = 5; JLbeatPlot({JbB2_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % changed ..
ifreq = 8; JLbeatPlot({JbB2_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % .. ~linear
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % big I change ..
ifreq = 11; JLbeatPlot({JbB2_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % .. more like shrinking

% RG10209 >>> cell 4 ("spikes + EPSPs", "beats")  BEAUTY: Many spikes, stable; good binaurality.
%      60 dB 2xbin/ipsi/contra 
% -11  4-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -12  4-2-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -13  4-3-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
% -14  4-4-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% 2 x NTD
JLreadBeats('RG10209', 4);
% multiple peaks at low freqs
% 200 Hz
ifreq=2; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 200 Hz  R = [0.78 0.62 0.49] 98 spikes; rho=0.45; 
ifreq=2; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3) ;% amazingly stable
% 300 Hz. Contra locks better. Note CHIST ~ waveform similarity
ifreq =3; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]); % 300 Hz  R = [0.56 0.88 0.45] 305 spikes; rho=0.53;
ifreq=3; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % again very stable
% 400 Hz. Peak splitting in both ears; CHIST ~ waveform similarity I/C/B. Bin. hist has 4 peaks.  
ifreq=4; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 400 Hz  R = [0.58 0.58 0.27] 169 spikes; rho=0.56;
ifreq=4; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % again very stable
% 500 Hz. Peaks start to merge.
ifreq=5; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 500 Hz  R = [0.79 0.78 0.62] 230 spikes; rho=0.55;
ifreq=5; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % amazingly stable
% 600 Hz. mono-peaked CHIST. % strange I timing of STA.
ifreq=6; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 600 Hz  R = [0.9 0.81 0.74] 463 spikes; rho=0.48;
ifreq=6; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % very stable
% 700 Hz; STA: I timing on 2nd, somewhat weaker peak; C timing a bit late.
ifreq=7; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 700 Hz  R = [0.87 0.84 0.74] 511 spikes; rho=0.39;
ifreq=7; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3);% very stable
% 800 Hz; both I&C waveforms are double-peaked. STA timed on slightly smaller
%         2nd I peak; biggest C peak
ifreq=8; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 800 Hz  R = [0.87 0.76 0.7] 337 spikes; rho=0.26;
ifreq=8; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % stable; slight C change
% 900 Hz; I still double peaked..
ifreq=9; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]); % 900 Hz  R = [0.81 0.77 0.67] 239 spikes; rho=0.24;
ifreq=9; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % okay
% 1000 Hz I peaks are fusing; repeat has consistent changes: I waveform (relative peak sizes),
%         CHIST shape (skewedness) and STA (re peak sizes)
ifreq=10; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]); % 1000 Hz  R = [0.81 0.78 0.66] 186 spikes; rho=0.37;
ifreq=10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % decent. I peak height change
ifreq=10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 6); % CHIST effect of I peak change
ifreq=10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 5); % STA effect of I peak change
% 1100 Hz single peaks from here
ifreq=11; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]); % 1100 Hz  R = [0.79 0.75 0.62] 137 spikes; rho=0.48;
ifreq=11; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) }, 3); % stable
% 1200 Hz
ifreq=12; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 1200 Hz  R = [0.75 0.73 0.54] 98 spikes; rho=0.52;
ifreq=12; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq)}, 3); % stable
% 1300 Hz
ifreq=13; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 1300 Hz  R = [0.72 0.57 0.6] 129 spikes; rho=0.48;
ifreq=13; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq)}, 3); % stable
% 1400 Hz
ifreq=14; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 1400 Hz  R = [0.56 0.44 0.51] 79 spikes; rho=0.44;
ifreq=14; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq)}, 3);  % stable
% 1500 Hz
ifreq=15; JLbeatPlot([JbB_60dB(ifreq) JbB2_60dB(ifreq)]);  % 1500 Hz  R = [0.48 0.5 0.43] 55 spikes; rho=0.31;
ifreq=15; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq)}, 3); % stable
% across-frequency comparisons. Any systematic patterns??
Xfreq = {[JbB_60dB(2) JbB2_60dB(2)] [JbB_60dB(3) JbB2_60dB(3)] [JbB_60dB(4) JbB2_60dB(4)] [JbB_60dB(5) JbB2_60dB(5)] [JbB_60dB(6) JbB2_60dB(6)] [JbB_60dB(7) JbB2_60dB(7)]  };
JLbeatPlot(Xfreq,3); subsaf
JLbeatPlot(Xfreq,5); subsaf % STA
JLbeatPlot(Xfreq,6); subsaf % CHIST
Xfreq2 = {[JbB_60dB(8) JbB2_60dB(8)] [JbB_60dB(9) JbB2_60dB(9)] [JbB_60dB(10) JbB2_60dB(10)] [JbB_60dB(11) JbB2_60dB(11)] [JbB_60dB(12) JbB2_60dB(12)] [JbB_60dB(13) JbB2_60dB(13)]  };
JLbeatPlot(Xfreq2,3); subsaf
JLbeatPlot(Xfreq2,5); subsaf % STA
JLbeatPlot(Xfreq2,6); subsaf % CHIST
% 200 Hz  R = [0.78 0.62 0.49] 98 spikes; rho=0.45;
% 300 Hz  R = [0.56 0.88 0.45] 305 spikes; rho=0.53;
% 400 Hz  R = [0.58 0.58 0.27] 169 spikes; rho=0.56;
% 500 Hz  R = [0.79 0.78 0.62] 230 spikes; rho=0.55;
% 600 Hz  R = [0.9 0.81 0.74] 463 spikes; rho=0.48;
% 700 Hz  R = [0.87 0.84 0.74] 511 spikes; rho=0.39;
% 800 Hz  R = [0.87 0.76 0.7] 337 spikes; rho=0.26;
% 900 Hz  R = [0.81 0.77 0.67] 239 spikes; rho=0.24;
% 1000 Hz  R = [0.81 0.78 0.66] 186 spikes; rho=0.37;
% 1100 Hz  R = [0.79 0.75 0.62] 137 spikes; rho=0.48;
% 1200 Hz  R = [0.75 0.73 0.54] 98 spikes; rho=0.52;
% 1300 Hz  R = [0.72 0.57 0.6] 129 spikes; rho=0.48;
%===
% linearity. Are high freqs a more stringent test? More Noisy?
ifreq=2; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);  % very linear
ifreq=3; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);  % ditto
ifreq=4; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % amazingly linear
ifreq=5; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % ditto
ifreq=6; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % slight mutaul suppr
ifreq=7; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % slight mutaul suppr
ifreq=8; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % pretty linear
ifreq=9; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr
ifreq=10; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % small changes
ifreq=11; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr
ifreq=12; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr
ifreq=13; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % I->C suppr
ifreq=14; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % C suppr; I change
ifreq=14; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr
ifreq=15; JLbeatPlot({[JbB_60dB(ifreq) JbB2_60dB(ifreq)] JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr

% RG10209 >>> cell 5 ("spikes + EPSPs, lost cell after 500 Hz")  
%      60 dB bin
% -17  5-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              7* x 1 x 6000/10000 ms B
JLreadBeats('RG10209', 5);
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.48 0.82 0.53] 17 spikes; rho=0.44; usual multpeak stuff
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.57 0.76 0.5] 19 spikes; rho=0.58; peaks fusing
JLbeatPlot(JbB_60dB(4)); % 400 Hz  R = [0.83 0.87 0.69] 90 spikes; rho=0.6;
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.86 0.9 0.75] 193 spikes; rho=0.59;
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.75 0.91 0.65] 197 spikes; rho=0.55;
% 700 Hz has an artefact due to losing the cell
JLbeatPlot(JbB_60dB(7)); % 700 Hz  R = [0.88 0.88 0.79] 196 spikes; rho=0.41;

% RG10209 cell 6 ("lost cell")  
%      60 dB bin
% -18  6-2-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              3* x 1 x 6000/10000 ms B
JLreadBeats('RG10209', 6);
JLbeatPlot(JbB_60dB(1));%100 Hz  R = [0.62 0.61 0.3] 9 spikes; rho=0.015; NOISE
JLbeatPlot(JbB_60dB(2));%200 Hz  R = [0.23 0.97 0.2] 147 spikes; rho=0.31; very sharp C tuning. Spikes not I tuned.
JLbeatPlot(JbB_60dB(3));%300 Hz  R = [0.57 0.87 0.6] 7 spikes; rho=0.48; same phenomenon

% RG10209 cell 7 ("spikes + EPSPs", "instable signal")  
%      60 dB bin/ipsi/contra
% -19  7-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat              15 x 1 x 6000/10000 ms B
% -20  7-2-BFS 100*:100*:1500 Hz 60|-60 dB 4 Hz beat              9* x 1 x 6000/10000 ms B
% -21  7-3-BFS 100*:100*:1500 Hz -60|60 dB 4 Hz beat              15 x 1 x 6000/10000 ms B
error('Bookkeeping problems');
% Multiple ABF recordings matching 'RG10209/7-1-BFS'.
% No ABF recordings matching 'RG10209/7-2-BFS'.
% No ABF recordings matching 'RG10209/7-3-BFS'.
qq=[];

function RG10214
% RG10214 >>> cell 1 ("spikes + EPSPs" "lost cell during recording")  
%      60 dB bin
% -1   1-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B% RG216B >>>  cell 1 ("mini spikes, large EPSPS")  
JLreadBeats('RG10214', 1);
% contra dominates at 200 Hz and up
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.75 0.99 0.75] 7 spikes; rho=0.65;
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.99 0.71 0.67] 9 spikes; rho=0.68; ipsi triggered on 2nd, smaller pk
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.9 0.92 0.84] 18 spikes; rho=0.66; mono-peaked contra dominates
% spike triggering becomes somewhat problematic at 500 Hz and up
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.8 0.96 0.71] 4 spikes; rho=0.61;  ipsi/contra balance
JLbeatPlot(JbB_60dB(6)); % 600 Hz  R = [0.91 1 0.89] 2 spikes; rho=0.57; ipsi/contra balance
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.55;  ipsi dominates
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0 0 0] 0 spikes; rho=0.53;  ipsi dominates
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0 0 0] 1 spikes; rho=0.48;  ipsi dominates
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.39;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0 0 0] 0 spikes; rho=0.28; ipsi&contra balance
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0 0 0] 0 spikes; rho=0.18; contra dom.
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0 0 0] 0 spikes; rho=0.094; contra dominates
%
% RG10214 >>> cell 2 ("strange waveform", "beats")  
%      60 dB 2xbin/ipsi/contra
% -2   2-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -3   2-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -4   2-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -5   2-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
JLreadBeats('RG10214', 2);
%  mainly contra driven. Doesn't look too useful.
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.82 0.79 0.55] 12 spikes; rho=0.084;
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.27 0.82 0.38] 33 spikes; rho=0.11;
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.59 0.67 0.94] 30 spikes; rho=0.11;
% consistent repeat
JLbeatPlot(JbB2_60dB(4));  % 400 Hz  R = [0.36 0.43 0.92] 49 spikes; rho=0.12;
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.45 0.39 0.91] 30 spikes; rho=0.065;
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.65 0.65 0.87] 24 spikes; rho=0.04;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.048;
% linearity test 
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % contra only
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % well linear
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % still okayish

% RG10214 >>> cell 3 ("EPSPS no clear Aps")  
%      60 dB bin/ipsi/contra
% -6   3-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -7   3-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -8   3-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
JLreadBeats('RG10214', 3); % CHECK AP thr!
%
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0.97 0.68 0.81] 3 spikes; rho=0.4; sinusoidal
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.6 0.35 0.54] 5 spikes; rho=0.44; Multiple peaks, but not pronounced
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.74 0.81 0.5] 16 spikes; rho=0.57;Multiple peaks, but not pronounced
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.82 0.9 0.69] 62 spikes; rho=0.63; triggering to 2nd weak ipsi peak
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.5 0.91 0.49] 51 spikes; rho=0.56; contra dominates
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.6 0.85 0.46] 53 spikes; rho=0.55; equal-sized ipsi peaks
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.83 0.82 0.74] 48 spikes; rho=0.49; 2nd ipsi peak now bigger
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0.82 0.75 0.63] 43 spikes; rho=0.43;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.71 0.68 0.43] 30 spikes; rho=0.43;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.6 0.75 0.5] 28 spikes; rho=0.42;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.46 0.78 0.36] 42 spikes; rho=0.44;
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.33 0.82 0.38] 19 spikes; rho=0.41;
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0.7 0.68 0.45] 15 spikes; rho=0.37;
JLbeatPlot(JbB_60dB(14));  % 1400 Hz  R = [0.65 0.71 0.57] 17 spikes; rho=0.29;
JLbeatPlot(JbB_60dB(15));  % 1500 Hz  R = [0.83 0.69 0.58] 18 spikes; rho=0.2;
% change of multiple ipsi peaks w freq
JLbeatPlot({JbB_60dB(4) JbB_60dB(5) JbB_60dB(6) JbB_60dB(7) },3); xlim([0 10]); subsaf
% linearity test 
ifreq = 2; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % ipsi somewhat changed
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % very linear
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % pretty linear
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % pretty linear
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % pretty linear
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % contra somewhat facilitated?
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % contra somewhat facilitated?
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % contra facilitated
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % contra facilitated
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);

% RG10214 >>> cell 4 ("spikes + EPSPs", "beats and inhibition", "heart beats?")  
%      60 dB 2xbin/ipsi/contra
% -10  4-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
% -11  4-2-BFS  100*:100*:1500 Hz 60|-60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -12  4-3-BFS  100*:100*:1500 Hz -60|60 dB 4 Hz beat            15 x 1 x 6000/10000 ms B
% -13  4-4-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
JLreadBeats('RG10214', 4);
% spikes look nice
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0.47 0.48 0.15] 202 spikes; rho=0.39;
% note similarity of waveforms & CHIST in next rec
JLbeatPlot(JbB_60dB(2)); % 200 Hz  R = [0.15 0.63 0.081] 195 spikes; rho=0.54; all very multipeaked; 
JLbeatPlot(JbB_60dB(3));  % 300 Hz  R = [0.7 0.61 0.38] 151 spikes; rho=0.54; very multipeaked again
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.92 0.76 0.68] 232 spikes; rho=0.51; timing on 2nd, big ipsi peak
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.89 0.84 0.73] 161 spikes; rho=0.46; ipsi 1s peak reduced
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.84 0.85 0.7] 141 spikes; rho=0.45; "textbook"
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.82 0.89 0.71] 119 spikes; rho=0.45; "textbook"
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0.76 0.8 0.62] 60 spikes; rho=0.42;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.73 0.77 0.52] 35 spikes; rho=0.3;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.61 0.81 0.49] 17 spikes; rho=0.22;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.64 0.74 0.53] 21 spikes; rho=0.15;
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.38 0.66 0.38] 11 spikes; rho=0.14;
JLbeatPlot(JbB_60dB(13));  % 1300 Hz  R = [0.36 0.62 0.23] 15 spikes; rho=0.11;
JLbeatPlot(JbB_60dB(14));  % 1400 Hz  R = [0.5 0.72 0.4] 12 spikes; rho=0.092;
JLbeatPlot(JbB_60dB(15));  % 1500 Hz  R = [0.57 0.79 0.46] 21 spikes; rho=0.066;
% reproducibility?
ifreq=2; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) },3); % 200 Hz stable
ifreq=4; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) },3); % ipsi changed
ifreq=6; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) },3); % everything shrunk
ifreq=8; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) },3) ; % ditto
ifreq=10; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB2_60dB(ifreq) },3) % ditto
% linearity
ifreq = 2; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; ipsi change
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; bg ipsi change
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); ylim([-1 0]);% mutual suppr & shape change
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); ylim([-1 0]);% ditto
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); ylim([-1 0]);% big mutual suppr
ifreq = 12; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); ylim([-1 0]);% big mutual suppr

% RG10214 >>> cell 5 ("losing cell")  strange waveforms
%      60 dB bin
% -14  5-1-BFS  100*:100*:1500 Hz 60 dB     4 Hz beat            15 x 1 x 6000/10000 ms B
JLreadBeats('RG10214', 5);
% spikes look strange; not very useful
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.65 0.43 0.4] 7 spikes; rho=0.025;
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.14 0.12 0.98] 8 spikes; rho=0.073; 
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.63 0.28 0.7] 13 spikes; rho=0.04;
JLbeatPlot(JbB_60dB(8)); % 800 Hz  R = [0.73 0.51 0.6] 13 spikes; rho=0.065;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.2 0.58 0.21] 5 spikes; rho=0.049;

% RG10214 >>> cell 6 ("spikes + EPSPs")  
%      10 dB bin/ipsi/contra
%      20 dB bin/ipsi/contra
%      30 dB bin/ipsi/contra
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB 3xbin/ipsi/contra
%      70 dB bin/ipsi/contra
%      80 dB bin/ipsi/contra
%   NTD
JLreadBeats('RG10214', 6);
%
JLbeatPlot(JbB_60dB(1));  % 100 Hz  R = [0.092 0.91 0.14] 5 spikes; rho=0.18; contra driven
JLbeatPlot(JbB_60dB(2));  % 200 Hz  R = [0.38 0.97 0.37] 27 spikes; rho=0.39; contra driven
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0.76 0.84 0.62] 22 spikes; rho=0.47;
JLbeatPlot(JbB_60dB(4));  % 400 Hz  R = [0.81 0.98 0.77] 22 spikes; rho=0.5;
JLbeatPlot(JbB_60dB(5));  % 500 Hz  R = [0.95 0.77 0.7] 49 spikes; rho=0.5;
JLbeatPlot(JbB_60dB(6));  % 600 Hz  R = [0.91 0.85 0.73] 17 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(7));  % 700 Hz  R = [0.88 0.9 0.78] 31 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(8));  % 800 Hz  R = [0.85 0.92 0.76] 23 spikes; rho=0.44;
JLbeatPlot(JbB_60dB(9));  % 900 Hz  R = [0.71 0.89 0.62] 18 spikes; rho=0.42;
JLbeatPlot(JbB_60dB(10));  % 1000 Hz  R = [0.75 0.95 0.62] 9 spikes; rho=0.4;
JLbeatPlot(JbB_60dB(11));  % 1100 Hz  R = [0.74 0.87 0.65] 11 spikes; rho=0.35;
JLbeatPlot(JbB_60dB(12));  % 1200 Hz  R = [0.61 0.77 0.74] 3 spikes; rho=0.31;
% 50 dB SPL
JLbeatPlot(JbB_50dB(2));  % 200 Hz  R = [0.47 0.92 0.31] 11 spikes; rho=0.22;
JLbeatPlot(JbB_50dB(3));  % 300 Hz  R = [0.55 0.82 0.44] 18 spikes; rho=0.21; multi peaked
JLbeatPlot(JbB_50dB(4));  % 400 Hz  R = [0.77 0.94 0.73] 144 spikes; rho=0.54; textbook
JLbeatPlot(JbB_50dB(5));  % 500 Hz  R = [0.66 0.96 0.63] 220 spikes; rho=0.54; many more spikes than @60 dB
JLbeatPlot(JbB_50dB(6));  % 600 Hz  R = [0.83 0.92 0.77] 165 spikes; rho=0.53;
JLbeatPlot(JbB_50dB(7));  % 700 Hz  R = [0.85 0.91 0.79] 110 spikes; rho=0.53; textbook
JLbeatPlot(JbB_50dB(8));  % 800 Hz  R = [0.78 0.88 0.75] 58 spikes; rho=0.5; again nice
JLbeatPlot(JbB_50dB(9));  % 900 Hz  R = [0.8 0.86 0.69] 51 spikes; rho=0.46; still nice
JLbeatPlot(JbB_50dB(10));  % 1000 Hz  R = [0.68 0.87 0.6] 41 spikes; rho=0.42;
JLbeatPlot(JbB_50dB(11));  % 1100 Hz  R = [0.61 0.79 0.62] 24 spikes; rho=0.37;
JLbeatPlot(JbB_50dB(12));  % 1200 Hz  R = [0.75 0.82 0.55] 14 spikes; rho=0.32;
JLbeatPlot(JbB_50dB(13));  % 1300 Hz  R = [0.55 0.76 0.48] 18 spikes; rho=0.25;
JLbeatPlot(JbB_50dB(14));  % 1400 Hz  R = [0.5 0.86 0.65] 15 spikes; rho=0.19;
JLbeatPlot(JbB_50dB(15));  % 1500 Hz  R = [0.53 0.85 0.69] 8 spikes; rho=0.14;
% 40 dB SPL
JLbeatPlot(JbB_40dB(2));  % 200 Hz  R = [0.68 0.83 0.51] 5 spikes; rho=0.053;
JLbeatPlot(JbB_40dB(3));  % 300 Hz  R = [0.46 0.51 0.52] 8 spikes; rho=0.1;
JLbeatPlot(JbB_40dB(4));  % 400 Hz  R = [0.83 0.62 0.17] 3 spikes; rho=0.12;
JLbeatPlot(JbB_40dB(5));  % 500 Hz  R = [0.48 0.93 0.46] 43 spikes; rho=0.35;
JLbeatPlot(JbB_40dB(6));  % 600 Hz  R = [0.67 0.95 0.65] 127 spikes; rho=0.48;
JLbeatPlot(JbB_40dB(7));  % 700 Hz  R = [0.78 0.93 0.73] 121 spikes; rho=0.51;
JLbeatPlot(JbB_40dB(8));  % 800 Hz  R = [0.73 0.88 0.65] 80 spikes; rho=0.51;
JLbeatPlot(JbB_40dB(9));  % 900 Hz  R = [0.72 0.89 0.71] 45 spikes; rho=0.48;
JLbeatPlot(JbB_40dB(10));  % 1000 Hz  R = [0.7 0.84 0.56] 41 spikes; rho=0.43;
JLbeatPlot(JbB_40dB(11));  % 1100 Hz  R = [0.77 0.69 0.5] 28 spikes; rho=0.37;
JLbeatPlot(JbB_40dB(12));  % 1200 Hz  R = [0.72 0.93 0.72] 7 spikes; rho=0.32;
JLbeatPlot(JbB_40dB(13));  % 1300 Hz  R = [0.44 0.73 0.71] 9 spikes; rho=0.25;
JLbeatPlot(JbB_40dB(14));  % 1400 Hz  R = [0.69 0.79 0.75] 12 spikes; rho=0.19;
JLbeatPlot(JbB_40dB(15));  % 1500 Hz  R = [0.68 0.84 0.63] 7 spikes; rho=0.14;
% 30 dB SPL
JLbeatPlot(JbB_30dB(2));  % 200 Hz  R = [0.76 0.86 0.98] 2 spikes; rho=0.01;
JLbeatPlot(JbB_30dB(3));  % 300 Hz  R = [0.49 0.32 0.31] 11 spikes; rho=0.016;
JLbeatPlot(JbB_30dB(4));  % 400 Hz  R = [0.41 0.18 0.41] 7 spikes; rho=0.025;
JLbeatPlot(JbB_30dB(5));  % 500 Hz  R = [0.63 0.19 0.27] 7 spikes; rho=0.049;
JLbeatPlot(JbB_30dB(6));  % 600 Hz  R = [0.26 0.89 0.24] 38 spikes; rho=0.22;  contra driven
JLbeatPlot(JbB_30dB(7));  % 700 Hz  R = [0.56 0.86 0.49] 64 spikes; rho=0.35;  contra dominates
JLbeatPlot(JbB_30dB(8));  % 800 Hz  R = [0.81 0.88 0.69] 96 spikes; rho=0.45;
JLbeatPlot(JbB_30dB(9));  % 900 Hz  R = [0.57 0.86 0.44] 58 spikes; rho=0.41;
JLbeatPlot(JbB_30dB(10));  % 1000 Hz  R = [0.68 0.83 0.63] 32 spikes; rho=0.4;
JLbeatPlot(JbB_30dB(11));  % 1100 Hz  R = [0.68 0.68 0.58] 22 spikes; rho=0.37;
JLbeatPlot(JbB_30dB(12));  % 1200 Hz  R = [0.72 0.84 0.64] 19 spikes; rho=0.3;
JLbeatPlot(JbB_30dB(13));  % 1300 Hz  R = [0.54 0.73 0.53] 13 spikes; rho=0.22;
JLbeatPlot(JbB_30dB(14));  % 1400 Hz  R = [0.54 0.89 0.43] 9 spikes; rho=0.18;
JLbeatPlot(JbB_30dB(15));  % 1500 Hz  R = [0.25 0.73 0.5] 12 spikes; rho=0.15;
% 20 dB SPL
JLbeatPlot(JbB_20dB(2));  % 200 Hz  R = [0.66 0.19 0.34] 5 spikes; rho=0.0032;
JLbeatPlot(JbB_20dB(5));  % 500 Hz  R = [0.49 0.56 0.33] 5 spikes; rho=0.0034;
JLbeatPlot(JbB_20dB(6));  % 600 Hz  R = [0.62 0.73 0.26] 5 spikes; rho=0.042;
JLbeatPlot(JbB_20dB(7));  % 700 Hz  R = [0.4 0.71 0.42] 17 spikes; rho=0.094;
JLbeatPlot(JbB_20dB(8));  % 800 Hz  R = [0.4 0.88 0.42] 27 spikes; rho=0.21; % contra driven
JLbeatPlot(JbB_20dB(9));  % 900 Hz  R = [0.56 0.8 0.49] 37 spikes; rho=0.19; contra dominates
JLbeatPlot(JbB_20dB(10));  % 1000 Hz  R = [0.57 0.74 0.48] 26 spikes; rho=0.23; as binaural as it gets @ 20 dB
JLbeatPlot(JbB_20dB(11));  % 1100 Hz  R = [0.45 0.78 0.49] 17 spikes; rho=0.21;'
JLbeatPlot(JbB_20dB(12));  % 1200 Hz  R = [0.22 0.69 0.19] 19 spikes; rho=0.15;
JLbeatPlot(JbB_20dB(13));  % 1300 Hz  R = [0.69 0.89 0.79] 8 spikes; rho=0.11;
JLbeatPlot(JbB_20dB(14)); % 1400 Hz  R = [0.26 0.78 0.26] 12 spikes; rho=0.11;
JLbeatPlot(JbB_20dB(15));  % 1500 Hz  R = [0.25 0.87 0.032] 5 spikes; rho=0.097;
% 10 dB SPL. not driven; no consistent waveform, either.
JLbeatPlot(JbB_10dB(2));  % 200 Hz  R = [0 0 0] 1 spikes; rho=-0.0042;
JLbeatPlot(JbB_10dB(6));  % 600 Hz  R = [0 0 0] 0 spikes; rho=-0.0033;
JLbeatPlot(JbB_10dB(7));  % 700 Hz  R = [0 0 0] 0 spikes; rho=0.00059;
JLbeatPlot(JbB_10dB(8));  % 800 Hz  R = [0 0 0] 0 spikes; rho=0.0047;
JLbeatPlot(JbB_10dB(9));  % 900 Hz  R = [0 0 0] 1 spikes; rho=0.0071;
JLbeatPlot(JbB_10dB(12));  % 1200 Hz  R = [0 0 0] 0 spikes; rho=-0.0021;
% 70 dB SPL.
JLbeatPlot(JbB_70dB(1));  % 100 Hz  R = [0.26 0.55 0.31] 8 spikes; rho=0.25; not much
JLbeatPlot(JbB_70dB(2));  % 200 Hz  R = [0.51 0.9 0.4] 122 spikes; rho=0.63; binaural, but C dominates
JLbeatPlot(JbB_70dB(3));  % 300 Hz  R = [0.25 0.52 0.11] 55 spikes; rho=0.49; not very binaural, but note CHST
JLbeatPlot(JbB_70dB(4)); % 400 Hz  R = [0.59 0.7 0.39] 27 spikes; rho=0.43;  C dominates
JLbeatPlot(JbB_70dB(5));  % 500 Hz  R = [0.69 0.84 0.47] 12 spikes; rho=0.47;
JLbeatPlot(JbB_70dB(6));  % 600 Hz  R = [0.79 0.9 0.69] 39 spikes; rho=0.45; more binaural
JLbeatPlot(JbB_70dB(7));  % 700 Hz  R = [0.73 0.87 0.66] 41 spikes; rho=0.47;
JLbeatPlot(JbB_70dB(8));  % 800 Hz  R = [0.63 0.91 0.53] 48 spikes; rho=0.45;
JLbeatPlot(JbB_70dB(9));  % 900 Hz  R = [0.72 0.89 0.65] 36 spikes; rho=0.41;
JLbeatPlot(JbB_70dB(10)); % 1000 Hz  R = [0.59 0.77 0.54] 26 spikes; rho=0.35;
JLbeatPlot(JbB_70dB(11));  % 1100 Hz  R = [0.71 0.92 0.68] 20 spikes; rho=0.3;
JLbeatPlot(JbB_70dB(12));  % 1200 Hz  R = [0.64 0.71 0.52] 17 spikes; rho=0.26;
JLbeatPlot(JbB_70dB(13));  % 1300 Hz  R = [0.57 0.81 0.49] 6 spikes; rho=0.21;
JLbeatPlot(JbB_70dB(14));  % 1400 Hz  R = [0.41 0.73 0.35] 9 spikes; rho=0.17;
JLbeatPlot(JbB_70dB(15));  % 1500 Hz  R = [0.4 0.68 0.15] 4 spikes; rho=0.12;
% 80 dB SPL
JLbeatPlot(JbB_80dB(1));  % 100 Hz  R = [0.24 0.96 0.2] 16 spikes; rho=0.51; C driven
JLbeatPlot(JbB_80dB(2));  % 200 Hz  R = [0.39 0.57 0.13] 137 spikes; rho=0.6; multipeaked binaurality. CHIST (mon+bin)
JLbeatPlot(JbB_80dB(3));  % 300 Hz  R = [0.6 0.75 0.45] 74 spikes; rho=0.51; CHIST does not show multi peaks of waveform
JLbeatPlot(JbB_80dB(4));  % 400 Hz  R = [0.91 0.88 0.8] 57 spikes; rho=0.48;
JLbeatPlot(JbB_80dB(5));  % 500 Hz  R = [0.85 0.79 0.62] 20 spikes; rho=0.49;
JLbeatPlot(JbB_80dB(6));  % 600 Hz  R = [0.75 0.87 0.52] 34 spikes; rho=0.48;
JLbeatPlot(JbB_80dB(7));  % 700 Hz  R = [0.69 0.89 0.56] 40 spikes; rho=0.47;
JLbeatPlot(JbB_80dB(8));  % 800 Hz  R = [0.73 0.84 0.62] 38 spikes; rho=0.46;
JLbeatPlot(JbB_80dB(9));  % 900 Hz  R = [0.88 0.82 0.76] 30 spikes; rho=0.4;
JLbeatPlot(JbB_80dB(10)); % 1000 Hz  R = [0.54 0.75 0.65] 10 spikes; rho=0.35;
JLbeatPlot(JbB_80dB(11)); % 1100 Hz  R = [0.67 0.88 0.58] 12 spikes; rho=0.29;
JLbeatPlot(JbB_80dB(12)); % 1200 Hz  R = [0.7 0.88 0.51] 11 spikes; rho=0.25;
JLbeatPlot(JbB_80dB(13)); % 1300 Hz  R = [0.41 0.86 0.31] 8 spikes; rho=0.22;
JLbeatPlot(JbB_80dB(14)); % 1400 Hz  R = [0.77 0.92 0.45] 2 spikes; rho=0.17;
JLbeatPlot(JbB_80dB(15));  % 1500 Hz  R = [0.51 0.31 0.54] 5 spikes; rho=0.12;

% repeats @ 60 dB
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 3); % relative pk-height changes
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 5); % much fewer spikes in the end 
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 6);
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 3);  % marked C changes
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 3); % overall ampl grows a lot
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbB2_60dB(ifreq) JbB3_60dB(ifreq) }, 3); % ditto
% pooling @ 60 dB
ifreq=4; JLbeatPlot([JbB_60dB(4) JbB2_60dB(4) JbB3_60dB(4) ]); % 400 Hz  R = [0.7 0.89 0.54] 141 spikes; rho=0.52;
% linearity @ 60 dB
ifreq = 3; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; shapes unaffected
ifreq = 4; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; shapes unaffected
ifreq = 5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; shapes unaffected
ifreq = 6; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mutual suppr; shapes unaffected
ifreq = 7; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % mut suppr; sliht shape change
ifreq = 8; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % stronger nonlin interaction
ifreq = 9; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3); % I "rectified by C", C suppressed
ifreq = 10; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);
ifreq = 11; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) }, 3);
% linearity @ other SPLs
% 80 dB
ifreq = 3; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3);%pretty linear; maybe some mutual suppr
ifreq = 5; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3);%ipsi suppressed; waveforms ~stable
ifreq = 7; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % mutual suppr & ipsi changed
ifreq = 7; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 5); % other view on nonlin @ 700 Hz
ifreq = 8; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % strong suppr & wave changes
ifreq = 8; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 5); % other view on nonlin @ 800 Hz
ifreq = 9; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % suppr & ipsi changes
ifreq = 10; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % strong I->C suppr
ifreq = 11; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % strong I->C suppr
ifreq = 12; JLbeatPlot({JbB_80dB(ifreq) JbI_80dB(ifreq) JbC_80dB(ifreq) }, 3); % strong I->C suppr+ wave changes
% 40 dB 
ifreq = 3; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % very linear
ifreq = 5; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % very linear
ifreq = 6; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
ifreq = 7; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
ifreq = 8; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
ifreq = 10; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
ifreq = 12; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
ifreq = 14; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) }, 3); % linear
% 20 dB
ifreq = 7; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3); % linear; slight facilitation of C
ifreq = 8; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3); % linear; slight facilitation of C
ifreq = 9; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3); % linear; slight facilitation of C
ifreq = 10; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3);% linear; facilitation of C
ifreq = 12; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3); % linear; facilitation of C
ifreq = 14; JLbeatPlot({JbB_20dB(ifreq) JbI_20dB(ifreq) JbC_20dB(ifreq) }, 3); % linear; facilitation of C
% SPL effects
% 500 Hz: quite big waveform changes
ifreq = 5; JLbeatPlot({JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_70dB(ifreq) JbB_80dB(ifreq) }, 3); subsaf
% 700 Hz: not dramatic; note C phase lag @ 80 dB; invariance also shown by CHISTs
ifreq = 7; JLbeatPlot({JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_70dB(ifreq) JbB_80dB(ifreq) }, 3); subsaf
ifreq = 7; JLbeatPlot({JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_70dB(ifreq) JbB_80dB(ifreq) }, 6); subsaf
% nonmonotonic changes in spike rate % Rbeat
% 700 Hz  30 dB; R = [0.56 0.86 0.49] 64 spikes; rho=0.35;
% 700 Hz  40 dB; R = [0.78 0.93 0.73] 121 spikes; rho=0.51;
% 700 Hz  50 dB; R = [0.85 0.91 0.79] 110 spikes; rho=0.53;
% 700 Hz  60 dB; R = [0.88 0.9 0.78] 31 spikes; rho=0.48;
% 700 Hz  70 dB; R = [0.73 0.87 0.66] 41 spikes; rho=0.47;
% 700 Hz  80 dB; R = [0.69 0.89 0.56] 40 spikes; rho=0.47;
qq=[]

function RG10216B
% ======================RG10216B ==================
% RG216B >>>  cell 1 ("mini spikes, large EPSPS")  
%      70 dB 2xbin/contra/ipsi % few spikes, but good phase locking
JLreadBeats('RG10216B', 1);
% few spikes, but good phase locking:
JLbeatPlot(JbB_70dB(5)); % 500 Hz  R = [0.93 0.91 0.83] 8 spikes; rho=0.52;
JLbeatPlot([JbB_70dB(5) JbB2_70dB(5)]); % 500 Hz  R = [0.95 0.83 0.79] 18 spikes; rho=0.51;  % Long STA latency..
JLbeatPlot(JbB_70dB(6)); % 600 Hz   R = [0.92 0.84 0.78] 12 spikes
JLbeatPlot(JbB_70dB(7)); % 700 Hz  R = [0.9 0.9 0.88] 12 spikes; rho=0.28; smaller ipsi peak triggers spikes??
JLbeatPlot(JbB_70dB(8)); % 800 Hz  R = [0.81 0.86 0.8] 4 spikes; rho=0.24;
JLbeatPlot(JbB_70dB(9)); % 900 Hz  R = [0.95 0.82 0.94] 3 spikes; rho=0.39; % contra is increasingly dominant
% linear binaural combination:
ifreq=5; JLbeatPlot({JbB_70dB(ifreq) JbI_70dB(ifreq) JbC_70dB(ifreq) },3); % 500 Hz

% RG216B >>> cell 2 ("mini spikes, large EPSPS")
%      40 dB bin/contra/ipsi 
%      50 dB bin/contra/ipsi 
%      60 dB 5xbin/contra/ipsi 
%      70 dB bin/contra/ipsi 
JLreadBeats('RG10216B', 2);
% decent phase locking, good binaural phase sensitivity
ifreq =2; JLbeatPlot(JbB_70dB(ifreq)); % 200 Hz  R = [0.56 0.82 0.48] 10 spikes; rho=0.46; mult contra pks
ifreq =3; JLbeatPlot(JbB_70dB(ifreq)); % 300 Hz  R = [0.82 0.83 0.7] 13 spikes; rho=0.53; mult contra pks
ifreq =4; JLbeatPlot(JbB_70dB(ifreq)); % 400 Hz  R = [0.86 0.8 0.56] 11 spikes; rho=0.58; 
ifreq =5; JLbeatPlot(JbB_70dB(ifreq)); % 500 Hz  R = [0.72 0.66 0.63] 19 spikes; rho=0.52; long STA latency
ifreq =6; JLbeatPlot(JbB_70dB(ifreq)); % 600 Hz  R = [0.71 0.71 0.8] 15 spikes; rho=0.51; ditto
ifreq =7; JLbeatPlot(JbB_70dB(ifreq)); % 700 Hz  R = [0.83 0.82 0.73] 11 spikes; rho=0.39; IPSI triggers via 2nd small peak
ifreq =8; JLbeatPlot(JbB_70dB(ifreq)); % 800 Hz  R = [0.29 0.34 0.52] 7 spikes; rho=0.28;
% at lower SPLs, contra is dominant
ifreq = 5; JLbeatPlot(JbB_40dB(ifreq)); % 500 Hz  R = [0 0 0] 0 spikes; rho=0.12;
ifreq = 5; JLbeatPlot(JbB_50dB(ifreq)); % 500 Hz  R = [0.84 0.9 0.57] 3 spikes; rho=0.42;
ifreq = 6; JLbeatPlot(JbB_50dB(ifreq)); % 600 Hz  R = [0.97 0.037 0.29] 2 spikes; rho=0.3;
% pool to get good averaging
ifreq=4; JLbeatPlot([JbB_60dB(ifreq) JbB_60dB2(ifreq) JbB_60dB3(ifreq) JbB_60dB4(ifreq) JbB_60dB5(ifreq)]);
% ... 400 Hz  R = [0.74 0.87 0.61] 52 spikes; rho=0.55; Note very good binaural linearity
ifreq=5; JLbeatPlot([JbB_60dB(ifreq) JbB_60dB2(ifreq) JbB_60dB3(ifreq) JbB_60dB4(ifreq) JbB_60dB5(ifreq)]); 
% ... 500 Hz  R = [0.76 0.87 0.71] 67 spikes; rho=0.53; very linear as well.
% Linearity confirmed by true linearity tests (using monaural stim):
ifreq=5; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 500 Hz
ifreq=5; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 500 Hz
ifreq=5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq) JbC_60dB(ifreq) },3); % 500 Hz
% a hint of nonlinear interaural interaction at 70 dB SPL
ifreq=5; JLbeatPlot({JbB_70dB(ifreq) JbI_70dB(ifreq) JbC_70dB(ifreq) },3); % 500 Hz
% asymmetric STA timing
ifreq = 3; JLbeatPlot(JbB_70dB(ifreq)); % 300 Hz  R = [0.82 0.83 0.7] 13 spikes; rho=0.53; contra is multipeaked
ifreq = 4; JLbeatPlot(JbB_70dB(ifreq)); % 400 Hz  R = [0.86 0.8 0.56] 11 spikes; rho=0.58;
ifreq = 5; JLbeatPlot(JbB_70dB(ifreq)); % 500 Hz  R = [0.72 0.66 0.63] 19 spikes; rho=0.52;
ifreq = 6; JLbeatPlot(JbB_70dB(ifreq)); % 600 Hz  R = [0.71 0.71 0.8] 15 spikes; rho=0.51;
ifreq = 7; JLbeatPlot(JbB_70dB(ifreq)); % 700 Hz  R = [0.83 0.82 0.73] 11 spikes; rho=0.39;
ifreq = 8; JLbeatPlot(JbB_70dB(ifreq)); % 800 Hz  R = [0.29 0.34 0.52] 7 spikes; rho=0.28;

function RG10219;
% =======================RG10219==========================
% >>> cell 1 ("spikes & EPSPs")
%      30 dB contra/ipsi
%      40 dB bin/ipsi/contra
%      50 dB bin/ipsi/contra
%      60 dB bin
%      70 dB bin
% classical coincidence detector, e.g. 50 dB, 500 Hz: 
%     (contra, ipsi, bin)  -> (15,0,50) spikes
%   500 Hz, 50 dB: APs occur during steep descent of waveform. Is that due 
%   to backpropagation latency, or might the downward slope (which is 
%   the most accurately timed feature of the waveform) actually be
%   triggering the APs? Compare to other conditions & cells.

JLreadBeats('RG10219', 1);

JLwaterfall(JbB_40dB);
JLwaterfall(JbB_50dB);
JLwaterfall(JbB_60dB);
JLwaterfall(JbB_70dB);
JLwaterfall(JbI_50dB);

% excellent binaural timing, indep of spike count
ifreq=5; JLbeatPlot(JbB_50dB(ifreq)); % 400 Hz  R = [1 0.96 0.94] 3 spikes; rho=0.3;
ifreq=6; JLbeatPlot(JbB_50dB(ifreq)); % 450 Hz  R = [0.9 0.93 0.85] 44 spikes; rho=0.5; IPSI lead in STA
ifreq=7; JLbeatPlot(JbB_50dB(ifreq)); % 500 Hz  R = [0.91 0.94 0.85] 44 spikes; rho=0.52;  IPSI lead in STA
ifreq=8; JLbeatPlot(JbB_50dB(ifreq)); % 550 Hz  R = [0.89 0.93 0.8] 25 spikes; rho=0.45;  IPSI lead in STA
ifreq=9; JLbeatPlot(JbB_50dB(ifreq)); % 600 Hz  R = [0.94 0.94 0.77] 6 spikes; rho=0.4;
ifreq=10; JLbeatPlot(JbB_50dB(ifreq)); % 650 Hz  R = [0.95 0.88 0.83] 7 spikes; rho=0.36;
ifreq=11; JLbeatPlot(JbB_50dB(ifreq)); % 700 Hz  R = [0.95 0.98 0.94] 7 spikes; rho=0.35;
ifreq=12; JLbeatPlot(JbB_50dB(ifreq)); % 750 Hz  R = [0 0 0] 0 spikes; rho=0.31;
ifreq=13; JLbeatPlot(JbB_50dB(ifreq)); % 800 Hz  R = [0 0 0] 1 spikes; rho=0.34;
% at 50 dB SPL, linearity (bin resp = ipsi + contra) reigns ...
pp=0; ifreq=1; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 200 Hz
pp=0; ifreq=2; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 250 Hz
pp=0; ifreq=3; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 300 Hz
pp=0; ifreq=4; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 350 Hz
pp=0; ifreq=5; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 400 Hz
pp=0; ifreq=6; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 450 Hz
pp=0; ifreq=7; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 500 Hz
pp=0; ifreq=8; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 550 Hz
pp=0; ifreq=9; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 600 Hz
pp=0; ifreq=10; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 650 Hz
pp=0; ifreq=11; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 700 Hz
pp=0; ifreq=12; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 750 Hz
% ... except maybe at the highest stim freq, where ipsi response is tiny
pp=0; ifreq=13; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); % 800 Hz

% at 40 dB SPL, S/N is poorer, but response again looks linear
pp=0; ifreq=1; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 200 Hz
pp=0; ifreq=2; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 250 Hz
pp=0; ifreq=3; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 300 Hz
pp=0; ifreq=4; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 350 Hz
pp=0; ifreq=5; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 400 Hz
pp=0; ifreq=6; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 450 Hz
pp=0; ifreq=7; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 500 Hz
pp=0; ifreq=8; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 550 Hz
pp=0; ifreq=9; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 600 Hz
pp=0; ifreq=10; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 650 Hz
pp=0; ifreq=11; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 700 Hz
pp=0; ifreq=12; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 750 Hz
pp=0; ifreq=13; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % 800 Hz

% across level comparisons
ifreq=1; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 200 Hz
ifreq=2; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 250 Hz
ifreq=3; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 300 Hz
ifreq=4; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 350 Hz
ifreq=5; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 400 Hz
ifreq=6; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 450 Hz
ifreq=7; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 500 Hz
ifreq=8; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 550 Hz
ifreq=9; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 600 Hz
ifreq=10; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 650 Hz
ifreq=11; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 700 Hz
ifreq=12; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 750 Hz
ifreq=13; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_70dB(ifreq) },3); % 800 Hz

% RG10219 >>>  cell 2  % hardly spiking, but binaural waveforms
JLreadBeats('RG10219', 2);
% this looks somewhat binaurally-nonlinear ... but note the time lapse between bin and mon recs
ifreq=2; JLbeatPlot({JbB_50dB(ifreq) JbC_50dB(ifreq) JbI_50dB(ifreq)},3)
% and nonlinearit is not supported by spectrum:
ifreq=2; JLbeatPlot(JbB_50dB(ifreq),2)
ifreq=2; JLbeatPlot(JbB_60dB(ifreq))
% no spikes
JLbeatPlot(JbB_60dB(2)); % 200 Hz  R = [0 0 0] 0 spikes; rho=0.37; secondary peaks
JLbeatPlot(JbB_60dB(3)); % % 300 Hz  R = [0 0 0] 0 spikes; rho=0.32; large secondary peaks
JLbeatPlot(JbB_50dB(2)); % % 200 Hz  R = [0 0 0] 0 spikes; rho=0.35; weak secondary peaks
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0 0 0] 0 spikes; rho=0.25; secondary peaks
ifreq = 2; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) }, 3) % 200 Hz SPL dep of shape
ifreq = 3; JLbeatPlot({JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) }, 3) % 300 Hz SPL dep of shape

% >>>  RG10219 cell 3 ("spikes + EPSPS", weak binaurality; mostly ipsi driven; problematic spike triggering)  
%      30 dB bin
%      40 dB bin
%      50 dB bin
%      60 dB bin/ipsi/contra
%      70 dB bin
%      80 dB bin/ipsi/contra
JLreadBeats('RG10219', 3);
% 60 dB :hardly driven
ifreq=4; JLbeatPlot(JbB_60dB(ifreq)); % 400 Hz  R = [0.96 0.14 0.14] 15 spikes; rho=0.47;  %IPSI locked, not bin ..
ifreq=4; JLbeatPlot(JbI_60dB(ifreq)); % 400 Hz  R = [0.96 0.32 0.4] 6 spikes; rho=0.48; .. confirmed by IPSI stim
ifreq=5; JLbeatPlot(JbB_60dB(ifreq));  % 500 Hz  R = [0.94 0.38 0.33] 14 spikes; rho=0.51; again IPSI locked
% only @ higher freqs does the cell become somewhat binaural: (but note very low spike counts)
ifreq=8; JLbeatPlot(JbB_60dB(ifreq)); % 800 Hz  R = [0.93 0.93 1] 3 spikes; rho=0.43;
ifreq=9; JLbeatPlot(JbB_60dB(ifreq)); % 900 Hz  R = [0.96 0.47 0.36] 5 spikes; rho=0.37;
% hardly driven @ 30 dB SPL, except around 1000 Hz:
ifreq=10; JLbeatPlot(JbB_30dB(ifreq)); % 1000 Hz  R = [0.34 0.37 0.44] 11 spikes; rho=0.087;
ifreq=9; JLbeatPlot(JbB_40dB(ifreq)); % 900 Hz  R = [0 0 0] 1 spikes; rho=0.32;
ifreq=10; JLbeatPlot(JbB_40dB(ifreq)); % 1000 Hz  R = [0 0 0] 1 spikes; rho=0.25;


% >>>  RG10219 cell 4 ("mini spikes and big EPSPs" )  
%      40 dB bin
%      50 dB bin
%      60 dB bin
%      70 dB bin
% good binaurality
JLreadBeats('RG10219', 4);
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0.87 0.4 0.29] 16 spikes; rho=0.48; nice sharp single peaks
JLbeatPlot(JbB_60dB(4)); % 400 Hz  R = [0.85 0.9 0.72] 26 spikes; rho=0.66;
JLbeatPlot(JbB_60dB(5)); % 500 Hz  R = [0.92 0.75 0.78] 25 spikes; rho=0.62;
JLbeatPlot(JbB_60dB(6)); % % 600 Hz  R = [0.92 0.9 0.8] 13 spikes; rho=0.56; but spike waveforms and latency 
%          of STA suggest that the timing of the APs is problematic: bcause EPSP>AP   
JLbeatPlot(JbB_60dB(7)); % 700 Hz  R = [0.86 0.99 0.82] 3 spikes; rho=0.5; 3 well timed spikes. Same timing problem as above, though
JLbeatPlot(JbB_60dB(8)); % % 800 Hz  R = [0.84 0.95 0.78] 6 spikes; rho=0.42; spike triggering okay? Still good binaurality
% no spikes at 70 dB (nonmonotonicity reminding of Thomas' recordings)
JLbeatPlot(JbB_70dB(4)); % 400 Hz  R = [0 0 0] 0 spikes; rho=0.5;
% 50 dB SPL
JLbeatPlot(JbB_50dB(4));% 400 Hz  R = [0.95 0.59 0.52] 7 spikes; rho=0.58; ipsi rules
JLbeatPlot(JbB_50dB(5)) % 500 Hz  R = [0.91 0.93 0.94] 6 spikes; rho=0.61; a few well timed spikes; STA timing looks good
% 40 dB - no spikes
JLbeatPlot(JbB_40dB(5)); % 500 Hz  R = [0 0 0] 0 spikes; rho=0.25;

% >>>> RG10219 cell 5 ("spikes and EPSPs, but slowly losing cell" )  
%      30 dB bin
%      40 dB bin
%      50 dB bin
%      60 dB bin
JLreadBeats('RG10219', 5);
% 30 & 40 dB: good binaural sensitivity. Few spikes, but well timed ones.
JLbeatPlot(JbB_30dB(4)); % 400 Hz  R = [0.93 0.92 0.74] 5 spikes; rho=0.28;
JLbeatPlot(JbB_30dB(5)); % 500 Hz  R = [0.91 0.98 0.88] 4 spikes; rho=0.3;
%
JLbeatPlot(JbB_40dB(2)); % 200 Hz  R = [0.2 0.13 0.11] 16 spikes; rho=0.032;
JLbeatPlot(JbB_40dB(3)); % 300 Hz  R = [0.89 0.9 0.86] 7 spikes; rho=0.3;
JLbeatPlot(JbB_40dB(4)); % 400 Hz  R = [0.89 0.97 0.87] 14 spikes; rho=0.42;
JLbeatPlot(JbB_40dB(5)); % 500 Hz  R = [0.88 0.87 0.76] 11 spikes; rho=0.42;
JLbeatPlot(JbB_40dB(6)); % 600 Hz  R = [0.92 0.89 0.71] 3 spikes; rho=0.45;
JLbeatPlot(JbB_40dB(7)); % 700 Hz  R = [0 0 0] 1 spikes; rho=0.42;
% at 50 & 60 dB, spike rates do not increase. Still excellent binaural sensitivity
JLbeatPlot(JbB_50dB(2)); % 200 Hz  R = [0.67 0.75 0.44] 25 spikes; rho=0.45; (more spike re 60 dB, poorer ITD tuning)
JLbeatPlot(JbB_50dB(3)); % 300 Hz  R = [0.86 0.93 0.73] 10 spikes; rho=0.48;
JLbeatPlot(JbB_50dB(4)); % 400 Hz  R = [0.97 0.96 0.98] 5 spikes; rho=0.42;
JLbeatPlot(JbB_40dB(5)); % 500 Hz  R = [0 0 0] 1 spikes; rho=0.4;
%
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0.97 0.97 0.92] 16 spikes; rho=0.51;
JLbeatPlot(JbB_60dB(4)); % 400 Hz  R = [0.93 0.96 0.96] 4 spikes; rho=0.38;
JLbeatPlot(JbB_60dB(5)); % 500 Hz  R = [0.98 0.91 0.92] 8 spikes; rho=0.46;
JLbeatPlot(JbB_60dB(6)); % 600 Hz  R = [0 0 0] 1 spikes; rho=0.48;
JLbeatPlot(JbB_60dB(7)); % 700 Hz  R = [0 0 0] 0 spikes; rho=0.44;
JLbeatPlot(JbB_60dB(8)); % 800 Hz  R = [0 0 0] 1 spikes; rho=0.35;
JLbeatPlot(JbB_60dB(9)); % 900 Hz  R = [0 0 0] 0 spikes; rho=0.34;
JLbeatPlot(JbB_60dB(10)); % 1000 Hz  R = [0 0 0] 1 spikes; rho=0.24;

% >>>>> cell 6 ("spikes and EPSPs", "lost cell during BFS-6-2")  
%      60 dB bin/ipsi
% signs of nonlinearity; lots of spikes, but spike thr looks okay
JLreadBeats('RG10219', 6);
JLbeatPlot(JbB_60dB(2)); % 200 Hz  R = [0.63 0.67 0.39] 78 spikes; rho=0.53; CHST & MonAvgd waveforms look similar
JLbeatPlot(JbB_60dB(3)); % 300 Hz  R = [0.82 0.91 0.71] 210 spikes; rho=0.62; clear compression, see plots 2&3
JLbeatPlot(JbB_60dB(4)); % 400 Hz  R = [0.86 0.91 0.76] 115 spikes; rho=0.56;
JLbeatPlot(JbB_60dB(5)); % 500 Hz  R = [0.85 0.9 0.79] 87 spikes; rho=0.47;  % STA asymmetry "contra leading"
JLbeatPlot(JbB_60dB(6)); % 600 Hz  R = [0.79 0.86 0.66] 50 spikes; rho=0.44; (CONTRA dominant)
JLbeatPlot(JbB_60dB(7)); % 700 Hz  R = [0.84 0.82 0.8] 31 spikes; rho=0.47; % STA asymmetry "contra leading"
JLbeatPlot(JbB_60dB(8)); % 800 Hz  R = [0.83 0.82 0.64] 22 spikes; rho=0.4; % STA asymmetry "contra leading"
JLbeatPlot(JbB_60dB(9)); % 900 Hz  R = [0.91 0.63 0.68] 7 spikes; rho=0.32;  STA asymmetry "contra leading"
JLbeatPlot(JbB_60dB(10)); % 1000 Hz  R = [0.81 0.57 0.48] 8 spikes; rho=0.24;
JLbeatPlot(JbB_60dB(11)); % 1100 Hz  R = [0.74 0.8 0.54] 8 spikes; rho=0.2;
JLbeatPlot(JbB_60dB(12)); % 1200 Hz  R = [0.46 0.46 0.32] 8 spikes; rho=0.15;
% confirmation of the nonlinear binaural interaction, BUT note that cell was getting lost during ipsi rec.
% too bad contra is missing ..
ifreq=3; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq)  JbI_60dB(ifreq) }, 3)% 300 Hz 
ifreq=5; JLbeatPlot({JbB_60dB(ifreq) JbI_60dB(ifreq)  JbI_60dB(ifreq) }, 3)% 500 Hz 
% multiple peaks at low freqs
JLbeatPlot(JbI_60dB(2)) % 200 Hz
JLbeatPlot(JbC_60dB(2)) % 200 Hz % in contra & ipsi. The raw waveforms are clean and quite deterministic

function RG10223;
% =======================RG10223==========================
% RG10223 >>> cell 1 ("spikes & EPSPs" "gain was only 10, set to 50 after this recording")
% -3   1-3-BFS 100*:100*:2500 Hz 60 dB     4 Hz beat 25 x 1 x 6000/9500 ms B
%  Filter settings? (looks like HP but JL reports DC. Yeah looks like DC, see DC shifts btw repeats)
%      60 dB bin
JLreadBeats('RG10223', 1);
% not driven below 800 Hz, but nicely binaural & linear:
JLbeatPlot(JbB_60dB(5)); % 500 Hz; no spikes; rho=0.78; note high rho value
JLbeatPlot(JbB_60dB(6)); % 600 Hz; R = [0 0 0] 0 spikes; rho=0.72;
JLbeatPlot(JbB_60dB(7)); % 700 Hz;  R = [0 0 0] 0 spikes; rho=0.66;
% driven w (few, tiny spikes) @ 800 Hz and up. Note strange ipsi phase locking (dip) EI cell??
JLbeatPlot(JbB_60dB(8)); % 800 Hz;  R = [0.72 0.88 0.83] 10 spikes; rho=0.58;
JLbeatPlot(JbB_60dB(9)); % 900 Hz  R = [0.94 0.99 0.97] 5 spikes; rho=0.54
JLbeatPlot(JbB_60dB(10)); % 1000 Hz  R = [0.5 0.77 0.66] 8 spikes; rho=0.46
JLbeatPlot(JbB_60dB(11)); % 1100 Hz  R = [0.9 0.89 0.82] 5 spikes; rho=0.34; here ipsi fires at peak
JLbeatPlot(JbB_60dB(12)); % 1200 Hz  R = [0.93 0.84 0.93] 3 spikes; rho=0.19; ditto
JLbeatPlot(JbB_60dB(13)); % 1300 Hz  R = [0 0 0] 0 spikes; rho=0.17
JLbeatPlot(JbB_60dB(14)); % 1400 Hz  R = [0.86 0.9 0.69] 4 spikes; rho=0.23
JLbeatPlot(JbB_60dB(15)); % 1500 Hz  R = [0.15 0.88 0.26] 4 spikes; rho=0.27
JLbeatPlot(JbB_60dB(16)); % 1600 Hz  R = [0 0 0] 1 spikes; rho=0.27
JLbeatPlot(JbB_60dB(17)); % 1700 Hz   R = [0 0 0] 0 spikes; rho=0.21

% RG10223 >>> cell 2 ("spikes & EPSPs"; from 2-2-BFS: "isolation getting worse")
%      50 dB 2xbin/ipsi/contra
% -4   2-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -5   2-2-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -6   2-3-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -7   2-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 2);
% multiple peaks, absent or weak coincidence detection @ 200 Hz
JLbeatPlot(JbB_50dB(2)); % 200 Hz  R = [0.4 0.47 0.11] 102 spikes; rho=0.33
JLbeatPlot(JbB_50dB(3)); % 300 Hz   R = [0.7 0.4 0.32] 49 spikes; rho=0.38; strange contra timing
% good binaurality at 400 Hz and up. 
JLbeatPlot(JbB_50dB(4)); % 400 Hz   R = [0.91 0.85 0.75] 111 spikes; rho=0.59; strange contra timing
JLbeatPlot(JbB_50dB(5)); % 500 Hz   R = [0.86 0.86 0.74] 124 spikes; rho=0.56; contra timing on smaller 2nd "peak"?
JLbeatPlot(JbB_50dB(6)); % 600 Hz   R = [0.92 0.88 0.82] 170 spikes; rho=0.45; again contra timing slightly smaller 2nd peak
JLbeatPlot(JbB_50dB(7)); % 700 Hz   R = [0.85 0.86 0.74] 124 spikes; rho=0.43; contra timing on *larger* 2nd peak
JLbeatPlot(JbB_50dB(8)); % 800 Hz   R = [0.85 0.82 0.72] 67 spikes; rho=0.38; now contra timing is "normal"
JLbeatPlot(JbB_50dB(9)); % 900 Hz   R = [0.8 0.74 0.72] 38 spikes; rho=0.34; contra timing stil normal
JLbeatPlot(JbB_50dB(10)); % 1000 Hz   R = [0.8 0.65 0.49] 23 spikes; rho=0.31; reduced binaurality @ higher freqs
JLbeatPlot(JbB_50dB(11)); % 1100 Hz   R = [0.76 0.6 0.44] 18 spikes; rho=0.25; reduced binaurality @ higher freqs
JLbeatPlot(JbB_50dB(12)); % 1200 Hz   R = [0.73 0.59 0.35] 20 spikes; rho=0.25; ditto
JLbeatPlot(JbB_50dB(13)); % 1300 Hz   R = [0.63 0.71 0.62] 17 spikes; rho=0.23; not bad for 1300 Hz ..
JLbeatPlot(JbB_50dB(14)); % 1400 Hz   R = [0.99 0.98 1] 3 spikes; rho=0.2; again this sparse, but accurate coding
JLbeatPlot(JbB_50dB(15)); % 1500 Hz   R = [0 0 0] 1 spikes; rho=0.13; but now it's over
% evolution of double STA peak w freq. It's always the same 2nd contra
%    peak that drives the APs, whether it is smaller or bigger than 1st peak
JLbeatPlot({JbB_50dB(3) JbB_50dB(4) JbB_50dB(5)},5); % 300/400/500 Hz
JLbeatPlot({JbB_50dB(4) JbB_50dB(5) JbB_50dB(6)},5); % 400/400/500 Hz
JLbeatPlot({JbB_50dB(5) JbB_50dB(6) JbB_50dB(7)},5); % 500/600/700 Hz
JLbeatPlot({JbB_50dB(3) JbB_50dB(4) JbB_50dB(5) JbB_50dB(6) JbB_50dB(7) JbB_50dB(8) },5); ylim([-0.2 0.23]); subsaf
% how does his relate to timing re stim? Look at elative timing ipsi/contra
JLbeatPlot({JbB_50dB(3) JbB_50dB(4) JbB_50dB(5) JbB_50dB(6) JbB_50dB(7) JbB_50dB(8) },6);  subsaf
% the following looks like nonlinear ipsi/contra combination, but recall
% worsening recording conditions ("isolation getting worse"). Strangely
% enough, the contra recording (the last one) sometimes looks best
ifreq = 3; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
ifreq = 4; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
ifreq = 5; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
% the repeat of the bin rec (2-4-BFS) is more conistent with the monaurals:
ifreq = 4; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
ifreq = 5; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
ifreq = 6; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3);
% .. but it has no APs.
% One almost wonders if 2-1-BFS and the repeat 2-4-BFS are really the same cell:
ifreq = 2; JLbeatPlot({JbB_50dB(ifreq) JbB_50dB(ifreq) JbB2_50dB(ifreq) },3);%rho =[0.33 0.66]; contra is diff
ifreq = 3; JLbeatPlot({JbB_50dB(ifreq) JbB_50dB(ifreq) JbB2_50dB(ifreq) },3);%rho =[0.38 0.52]; strange differences
ifreq = 4; JLbeatPlot({JbB_50dB(ifreq) JbB_50dB(ifreq) JbB2_50dB(ifreq) },3);%rho=[0.59 0.76]; ipsi is different; contra same
ifreq = 6; JLbeatPlot({JbB_50dB(ifreq) JbB_50dB(ifreq) JbB2_50dB(ifreq) },3);%rho=[0.45 0.62]
ifreq = 8; JLbeatPlot({JbB_50dB(ifreq) JbB_50dB(ifreq) JbB2_50dB(ifreq) },3);%rho=[0.38 0.2]; everything has shrunk
% There are similarities, but it almost looks as if the weights of the
% different inputs have changed.

% RG10223 >>> cell 3 ("spikes + EPSPs")  
%      50 dB 2xbin/ipsi/contra
% -8   3-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -9   3-2-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -10  3-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -11  3-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 3);
% multiple peaks, absent or weak coincidence detection @ 200 Hz
JLbeatPlot(JbB_50dB(2)); % 200 Hz   R = [0.63 0.58 0.33] 52 spikes; rho=0.45
JLbeatPlot(JbB_50dB(3)); % 300 Hz   R = [0.58 0.38 0.18] 143 spikes; rho=0.27; note similarity btw CHST & av waveforms
% better binaurality at 400 Hz and up, many spikes
JLbeatPlot(JbB_50dB(4)); % 400 Hz   R = [0.63 0.87 0.49] 315 spikes; rho=0.49; strange timing to ipsi; 2nd peak?
JLbeatPlot(JbB_50dB(5)); % 500 Hz   R = [0.84 0.76 0.57] 550 spikes; rho=0.42; nonlinear bin comb, see also spectrum
JLbeatPlot(JbB_50dB(6)); % 600 Hz   R = [0.85 0.78 0.61] 532 spikes; rho=0.34
JLbeatPlot(JbB_50dB(7)); % 700 Hz   R = [0.81 0.77 0.58] 506 spikes; rho=0.28
% weird changes over time: averaged waveforms identical, but APs are completely gone
ifreq=2; JLbeatPlot({JbB_50dB(ifreq) JbB2_50dB(ifreq) JbB2_50dB(ifreq) },3); % 200 Hz
% contra "suppresses" ipsi, bit hinterms of synchr & av amplitude
ifreq=3; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 300 Hz; discrepancy btw ipsi ampl &..
ifreq=3; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 6); % 300 Hz .. its spike triggering ability
ifreq=4; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 400 Hz
ifreq=4; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 6); % 400 Hz
% at 200 Hz, ipsi & contra are more equivalent
ifreq=2; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 5); % 200 Hz
ifreq=2; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 6); % 200 Hz

% RG10223 >>> cell 4 ("spikes + EPSPs")  
%      50 dB 2xbin/ipsi/contra
%      40 dB bin/ipsi
% -12  4-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -13  4-2-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -14  4-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -15  4-4-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -16  4-5-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -17  4-6-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 4);
% at 200 Hz, contra dominates waveform and has two peaks. Ipsi has better locking, though 
JLbeatPlot(JbB_50dB(2)); % 200 Hz R = [0.77 0.5 0.36] 15 spikes; rho=0.44;
JLbeatPlot(JbB_50dB(3)); % 300 Hz R = [0.83 0.23 0.031] 40 spikes; rho=0.33;  not binaural 
JLbeatPlot(JbB2_50dB(3)); % 300 Hz R = [0.79 0.25 0.13] 55 spikes; rho=0.31; hardly binaural
% cell gets binaural @ 400 Hz and up
JLbeatPlot(JbB_50dB(4)); % 400 Hz  R = [0.88 0.88 0.74] 137 spikes; rho=0.57
% STA tming vries over time
JLbeatPlot(JbB_50dB(5)); % 500 Hz  R = [0.9 0.81 0.72] 114 spikes; rho=0.57
JLbeatPlot(JbB2_50dB(5)); % 500 Hz   R = [0.91 0.75 0.67] 204 spikes; rho=0.48
% ipsi contrib grows over time
JLbeatPlot(JbB_50dB(6)); % 600 Hz  R = [0.89 0.86 0.77] 107 spikes; rho=0.49
JLbeatPlot(JbB2_50dB(6)); % 600 Hz R = [0.88 0.83 0.73] 147 spikes; rho=0.42

JLbeatPlot(JbB_50dB(7)); % 700 Hz  R = [0.84 0.82 0.73] 81 spikes; rho=0.44
JLbeatPlot(JbB_50dB(8)); % 800 Hz  R = [0.84 0.76 0.61] 49 spikes; rho=0.42; asymmetric STA timing 
JLbeatPlot(JbB_50dB(9)); % 900 Hz  R = [0.85 0.74 0.68] 33 spikes; rho=0.35; asymmetric STA timing 
JLbeatPlot(JbB_50dB(10)); % 1000 Hz   R = [0.85 0.77 0.59] 24 spikes; rho=0.28
JLbeatPlot(JbB_50dB(11)); % 1100 Hz   R = [0.76 0.65 0.41] 18 spikes; rho=0.19; binaurality is fading
JLbeatPlot(JbB_50dB(12)); % 1200 Hz   R = [0.78 0.7 0.61] 20 spikes; rho=0.12;
JLbeatPlot(JbB_50dB(13)); % 1300 Hz   R = [0.87 0.65 0.56] 14 spikes; rho=0.087; 
JLbeatPlot(JbB_50dB(14)); % 1400 Hz   R = [0.6 0.67 0.64] 23 spikes; rho=0.073; small rho, but still binaural
JLbeatPlot(JbB_50dB(15)); % 1500 Hz   R = [0.66 0.49 0.27] 17 spikes; rho=0.07; binaurality almost gone
% ========changes over time==========
% ipsi & contra combine pretty linearly, but some mutual suppression,
% mostly contra suppresses ipsi. Absolute timing of peaks & Relative peak sizes vary over time.
% (Abs timing gets delayed). Also there seems to be sizeable DC shifts over time.
ifreq=4; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 400 Hz
ifreq=4; JLbeatPlot([JbB_50dB(ifreq) JbB2_50dB(ifreq)]); % 400 Hz; DC shift over time
ifreq=5; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 500 Hz
ifreq=5; JLbeatPlot([JbB_50dB(ifreq) JbB2_50dB(ifreq)]); % 500 Hz % DC shift over time
ifreq=6; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 600 Hz
ifreq=7; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 700 Hz
ifreq=8; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 800 Hz
ifreq=9; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) JbB2_50dB(ifreq)}, 3); subsaf; % 900 Hz
ifreq=9; JLbeatPlot([JbB_50dB(ifreq) JbB2_50dB(ifreq)]); % 900 Hz % DC shift over time
ifreq=5; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 3); % 500 Hz; early ipsi peak enhanced in bin
ifreq=6; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 3); % 600 Hz; same
ifreq=7; JLbeatPlot({JbB_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) }, 3); % 700 Hz;
% across-SPL comparison
ifreq=3; JLbeatPlot({JbB2_50dB(ifreq) JbB_40dB(ifreq)}, 3); % 300 Hz huge diff; rho=[0.31 0.077];
ifreq=5; JLbeatPlot({JbB2_50dB(ifreq) JbB_40dB(ifreq)}, 3); % 500 Hz marginal diff; rho=[0.48 0.44]
ifreq=7; JLbeatPlot({JbB2_50dB(ifreq) JbB_40dB(ifreq)}, 3); % 700 Hz marginal diff; rho=[0.41 0.4]


% RG10223 >>> cell 5 ("spikes + EPSPs, but losing cells slowly after 1000 Hz", isolation getting worse)  
%      50 dB bin
%      40 dB bin/ipsi/contra
% -18  5-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -19  5-2-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -20  5-3-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -21  5-4-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 5);
% No spikes except (?) @ 40 dB, f~500 Hz. 
JLbeatPlot(JbB_50dB(2)); % 200 Hz;  R = [0 0 0] 0 spikes; rho=0.64; linear, multi-peaked 
JLbeatPlot(JbB_50dB(3)); % 300 Hz; R = [0 0 0] 0 spikes; rho=0.47; 
JLbeatPlot(JbB_40dB(3)); % 300 Hz; R = [0 0 0] 0 spikes; rho=0.055; 10 dB less is a lot ...
JLbeatPlot(JbB_50dB(4)); % 400 Hz  R = [0 0 0] 0 spikes; rho=0.73; nice waveform bit no spikes
JLbeatPlot(JbB_40dB(4)); % 400 Hz  R = [0.57 0.92 0.44] 14 spks; rho=0.55; dubious spikes, but well timed to cntr
JLbeatPlot(JbB_50dB(5)); % 500 Hz  R = [0 0 0] 0 spikes; rho=0.68; ditto
JLbeatPlot(JbB_40dB(5)); % 500 Hz  R = [0.48 0.88 0.45] 28 spks; rho=0.65; dubious spikes
JLbeatPlot({JbB_40dB(5) JbB_50dB(5) JbB_50dB(5)}, 3); % 500 Hz 40 dB vs 50 dB; timing & size ipsi changes
JLbeatPlot(JbB_50dB(6)); % 600 Hz  R = [0 0 0] 0 spikes; rho=0.47; ipsi dominates
JLbeatPlot(JbB_40dB(6)); % 600 Hz   R = [0.15 0.85 0.53] 3 spikes; rho=0.59;
JLbeatPlot(JbB_50dB(7)); % 700 Hz  R = [0 0 0] 0 spikes; rho=0.34; ipsi dominates, contra double-pkd
JLbeatPlot(JbB_50dB(8)); % 800 Hz  R = [0 0 0] 0 spikes; rho=0.25; ipsi dominates, contra double-pkd
JLbeatPlot(JbB_50dB(9)); % 900 Hz  R = [0 0 0] 0 spikes; rho=0.21; ipsi dominates, contra double-pkd
JLbeatPlot(JbB_50dB(10)); % 1000 Hz  R = [0 0 0] 0 spikes; rho=0.17; ipsi dom.; 1 contra peak shrinking
JLbeatPlot(JbB_50dB(11)); % 1100 Hz  R = [0 0 0] 0 spikes; rho=0.13; ipsi reduced as well; cntra peaks fuse
JLbeatPlot({JbB_50dB(9) JbB_50dB(10) JbB_50dB(11) }, 3); % 3 previous wavforms together
JLbeatPlot(JbB_50dB(12)); % 1200 Hz ARTEFACT 
clampfit(JbB_50dB(12)); % 1200 Hz;  ARTEFACT: Vrec plunges at 4000 ms
clampfit(JbB_50dB(13)); % huge drift
clampfit(JbB_50dB(14)); % somewhat stabilizing
JLbeatPlot(JbB_50dB(14)); % 1400 Hz   R = [0 0 0] 0 spikes; rho=0.12; monopeaked
JLbeatPlot(JbB_50dB(15)); % 1500 Hz   R = [0 0 0] 0 spikes; rho=0.11; monopeaked
% with its rho-value of 0.73, this is a candidate for event analysis:
clampfit(JbB_50dB(4)); % but what are the negative excursions? IPSPs? missed spikes?

% RG10223 >>>cell 6 ("recordings after a Gigaseal, R=40Mohm, not standard waveforms")  
%      20 dB bin
%      30 dB bin
%      40 dB bin/ipsi/contra
%      50 dB bin
%      60 dB bin
%      80 dB bin
% -22  6-1-BFS 100*:100*:1500 Hz 60 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -23  6-2-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -24  6-3-BFS 100*:100*:1500 Hz 40 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -25  6-4-BFS 100*:100*:1500 Hz 40|-40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -26  6-5-BFS 100*:100*:1500 Hz -40|40 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -27  6-6-BFS 100*:100*:1500 Hz 30 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -28  6-7-BFS 100*:100*:1500 Hz 20 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% XXX -29  6-8-BFS 100*:100*:1500 Hz 80 dB     4 Hz beat 1* x 1 x 6000/9500 ms B aborted due to Voltage clamp
% -30  6-9-BFS 100*:100*:1500 Hz 80 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 6);

% spike detection is problematic with these tiny events. Not a very binaural cell, either 
JLbeatPlot(JbB_50dB(2)); % 200 Hz; R = [0.14 0.83 0.07] 18 spikes; rho=0.49; probl not real spikes
JLbeatPlot(JbB_50dB(3)); % 300 Hz;  R = [0.16 0.61 0.25] 26 spikes; rho=0.24; % doubly peaked inputs
JLbeatPlot(JbB_50dB(4)); % 400 Hz;  R = [0.11 0.93 0.089] 79 spikes; rho=0.58; % only ipsi is doubly peaked now
JLbeatPlot(JbB_50dB(5)); % 500 Hz; R = [0.11 0.92 0.049] 65 spikes; rho=0.57;
JLbeatPlot(JbB_50dB(6)); % 600 Hz;  R = [0.069 0.84 0.077] 39 spikes; rho=0.44;
JLbeatPlot(JbB_50dB(7)); % 700 Hz;  R = [0.18 0.82 0.2] 46 spikes; rho=0.27; doubly-pkd ipsi, "spikes" all synched to contra
JLbeatPlot(JbB_50dB(8)); % 800 Hz;  R = [0.14 0.8 0.16] 56 spikes; rho=0.14;
JLbeatPlot(JbB_50dB(11)); % 1100 Hz;  R = [0.24 0.79 0.22] 43 spikes; rho=0.091;
% 80 dB; still problematic spike detection
JLbeatPlot(JbB_80dB(1)); % 100 Hz   R = [0.29 0.72 0.35] 26 spikes; rho=0.8; high rho, multiple peaks
JLbeatPlot(JbB_80dB(2)); % 200 Hz  R = [0.47 0.67 0.18] 154 spikes; rho=0.74; high rho, multiple peaks
JLbeatPlot(JbB_80dB(3)); % 300 Hz   R = [0.17 0.49 0.18] 49 spikes; rho=0.6; multiple peaks
JLbeatPlot(JbB_80dB(5)); % 500 Hz  R = [0.16 0.85 0.15] 28 spikes; rho=0.57;
JLbeatPlot(JbB_80dB(7)); % 700 Hz   R = [0.17 0.73 0.18] 18 spikes; rho=0.29; fading multpeaked ipsi
JLbeatPlot(JbB_80dB(9)); % 900 Hz  R = [0.057 0.25 0.59] 8 spikes; rho=0.12; messy
% 20 dB
JLbeatPlot(JbB_20dB(2)); % 200 Hz   R = [0.48 0.74 0.077] 4 spikes; rho=-0.0023; below thr
JLbeatPlot(JbB_20dB(4)); % 400 Hz   R = [0.15 0.18 0.27] 6 spikes; rho=-0.0019; a glance of a response
JLbeatPlot(JbB_20dB(6)); % 600 Hz   R = [0.53 0.46 0.32] 5 spikes; rho=0.049; faintly contra driven
JLbeatPlot(JbB_20dB(8)); % 800 Hz   R = [0.53 0.46 0.32] 5 spikes; rho=0.049; faintly contra driven
% 30 dB
JLbeatPlot(JbB_30dB(2)); % 200 Hz  R = [0 0 0] 0 spikes; rho=0.015;
JLbeatPlot(JbB_30dB(4)); % 400 Hz  R = [0 0 0] 0 spikes; rho=0.079;
JLbeatPlot(JbB_30dB(6)); % 600 Hz  R = [0.99 0.65 0.76] 2 spikes; rho=0.23; % not really driven, but consistent waveform 
JLbeatPlot(JbB_30dB(8)); % 800 Hz R = [0 0 0] 1 spikes; rho=0.13; fading again. Contra dominant. 
% 40 dB
JLbeatPlot(JbB_40dB(2)); % 200 Hz   R = [0.13 0.29 0.67] 6 spikes; rho=0.11;
JLbeatPlot(JbB_40dB(4)); % 400 Hz   R = [0.56 0.53 0.61] 7 spikes; rho=0.15;
JLbeatPlot(JbB_40dB(6)); % 600 Hz   R = [0.13 0.67 0.13] 10 spikes; rho=0.35; contra dominant
JLbeatPlot(JbB_40dB(8)); % 800 Hz  R = [0.45 0.68 0.3] 12 spikes; rho=0.1; fading again
% change of waveforms w SPL; potentially interesting but hard to interpret ...
ifreq=2; JLbeatPlot({JbB_20dB(ifreq) JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_80dB(ifreq) },3); axis([0 7 -0.5 0.5]); subsaf
ifreq=4; JLbeatPlot({JbB_20dB(ifreq) JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_80dB(ifreq) },3); axis([0 7 -0.3 0.3]); subsaf
ifreq=5; JLbeatPlot({JbB_20dB(ifreq) JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_80dB(ifreq) },3); axis([0 7 -0.3 0.3]); subsaf
ifreq=6; JLbeatPlot({JbB_20dB(ifreq) JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq) JbB_60dB(ifreq) JbB_80dB(ifreq) },3); axis([0 7 -0.3 0.3]); subsaf
ifreq=7; JLbeatPlot({JbB_20dB(ifreq) JbB_30dB(ifreq) JbB_40dB(ifreq) JbB_50dB(ifreq)  JbB_80dB(ifreq) },3); axis([0 7 -0.3 0.3]); subsaf
% linearity of binaural combination
ifreq=2; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % good consistency
ifreq=3; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % good consistency
ifreq=4; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % reasonably consistent
ifreq=5; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % good consistency
ifreq=6; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % good consistency, incl weird multi peaks
ifreq=7; JLbeatPlot({JbB_40dB(ifreq) JbI_40dB(ifreq) JbC_40dB(ifreq) },3); % good consistency, incl weird multi peaks


% RG10223 cell 7 (""failed patch attempt: small upward spikes", "lots of inhibition, not clear beats")  
%      50 dB 2xbin/ipsi/contra
% -31  7-1-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -32  7-2-BFS 100*:100*:1500 Hz 50 dB     4 Hz beat 15 x 1 x 6000/9500 ms B
% -33  7-3-BFS 100*:100*:1500 Hz 50|-50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
% -34  7-4-BFS 100*:100*:1500 Hz -50|50 dB 4 Hz beat 15 x 1 x 6000/9500 ms B
JLreadBeats('RG10223', 7, 1)
JLbeatPlot(JbB_50dB(3)); % 300 Hz  R = [0 0 0] 0 spikes; rho=0.21;
JLbeatPlot(JbB_50dB(4)); % 400 Hz  R = [0 0 0] 0 spikes; rho=0.41;
JLbeatPlot(JbB_50dB(5)); % 500 Hz  R = [0 0 0] 0 spikes; rho=0.55; % sharp peaks
JLbeatPlot(JbB_50dB(6)); % 600 Hz  R = [0 0 0] 0 spikes; rho=0.53;
JLbeatPlot(JbB_50dB(7)); % 700 Hz  R = [0 0 0] 0 spikes; rho=0.48;
JLbeatPlot(JbB_50dB(8)); % 800 Hz  R = [0.96 0.17 0.11] 2 spikes; rho=0.43;
JLbeatPlot(JbB_50dB(9)); % 900 Hz  R = [0 0 0] 0 spikes; rho=0.48;
% decent consistency across repeats
ifreq = 3; JLbeatPlot({JbB_50dB(ifreq) JbB2_50dB(ifreq) },3); % 300 Hz
ifreq = 4; JLbeatPlot({JbB_50dB(ifreq) JbB2_50dB(ifreq) },3); % 400 Hz
ifreq = 5; JLbeatPlot({JbB_50dB(ifreq) JbB2_50dB(ifreq) },3); % 500 Hz
ifreq = 6; JLbeatPlot({JbB_50dB(ifreq) JbB2_50dB(ifreq) },3); % 600 Hz
% linearity of binaural combination
ifreq=2; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.15 0]); % slight reduction of contra
ifreq=3; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.15 0]); % linear
ifreq=4; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.2 0]); % linear
ifreq=5; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.25 0]); % linear
ifreq=6; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.25 0]); % linear
ifreq=7; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 20 -0.17 0]); % linear
ifreq=8; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 10 -0.17 0]); % linear
ifreq=9; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 10 -0.17 0]); % linear
ifreq=10; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 10 -0.17 0]); % linear
ifreq=11; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 10 -0.17 0]); % linear
ifreq=12; JLbeatPlot({JbB2_50dB(ifreq) JbI_50dB(ifreq) JbC_50dB(ifreq) },3); axis([0 10 -0.17 0]);  %linear













