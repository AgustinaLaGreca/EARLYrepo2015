% notes on JSS stimuli
% ============linear zwuis============
%  - with Ra=5 Mohm, Cm=60 pF, Rm=100 Mohm, the action, ie, Cm dependence,
%    is at 10-5000 Hz. Ampl and phase provide complementary info
%  - lin freq range: 10,20,50,100,200,500,1000,2000,5000 Hz
%  - take RMS of zwuis = 20 mV
%  - sweep with different DC offsets: -120:20:120 mV (13 cond)
%  - total length of stimulus waveform: 13 * 1 = 13 s  (1.3 Msamples)

% ===========nonlinear zwuis==========
%  - restrict range to 130 Hz; 500 Hz; 850 Hz respectively. (?)
%  - use 5 comp uberzwuis w unique DPs of order 1,2,3.
%  - use 120 mV max Voltage (?), 0 mV offset.
25.7804   
121.7542
707.9458


% ====6 components=====
% best_freq = [0    19    45    76    93   118]
% diff(best_freq) equals [19    26    31    17    25]
% ii+best _freq works for e.g. ii=41
%
% ====5 components=====
% best_freq = [0    28    52    77    97]
% diff(best_freq) equals [ 28    24    25    20]
%
% best_freq = [0    19    33    43    59]
% diff(best_freq) equals [19    14    10    16]
% ii+best _freq works for e.g. ii=18
% Use this to generate:
% freq_Hz = (18+best_freq)/1.6 ...
%    = [11.2500   23.1250   31.8750   38.1250   48.1250]; % 1.6-s cycle
% freq_Hz = (18+best_freq)*2 ...
%    = [36    74   102   122   154]; % 0.5-s cycle
% freq_Hz = (18+best_freq)*12 ...
%    = [216   444   612   732   924]; % 84-ms cycle






