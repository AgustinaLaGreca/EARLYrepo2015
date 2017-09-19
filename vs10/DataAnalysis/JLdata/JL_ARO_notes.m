% JL ARO notes

%---------Possible figures:---------
% 1. The inputs to MSO neurons can be measured with juxtacellular recordings.
% A. Search stimuli (field potentials).
% B. Histology.
% C. Spontaneous responses during juxtacellular recordings. 
% D. Analysis of kinetics of spontaneous events. 
% E. Tone responses (monaural). 
% F. Analysis showing that the size of the eEPSP predicts the triggering of 
%    action potentials during monaural stimulation.
% G. Comparison spontaneous EPSPs with suprathreshold APs?
% H. Population data?
% 
% 2. Input at low frequencies can consist of a stereotyped sequence of short events. 
% A. Monaural input example.
% B. Averaged (phase-locked) response. 
% C. Principal component analysis. 
% D. Waveform analysis, providing more evidence for the (lack of?) 
%    interdependence of these different peaks. 
% E. Population data?
% F. Comparison with SBC firing behaviour?
% 
% 3. Linear interaction between inputs from both ears during binaural beat presentations.
% A. Binaural beat example response.
% B. Averaged response for each ear.
% C. Comparison with monaural inputs (if not shown in Fig. 2).
% D. Deviation from linearity. 
% E. Spike triggering.
% F. Population data on linearity.
% 
% 4. Mechanism of suppression of out-of-phase inputs. (if applicable)
% A. Vector strength?
% B. Waveform analysis?
% C. Model?
% D. ?
% 
% Not illustrated: tuning, inhibitie, (auto)correlations, etc. 
% What is still unknown about the cell physiological mechanism underlying ITD detection by MSO neurons?
% 1.	What is the time course of the EPSPs in vivo? To what extent is this time course determining the ITD tuning curves?
% 2.	What is the timing of the monaural inputs? Does this determine phase-locking?
% 3.	How many inputs have to be simultaneously active to trigger a spike? How constant is the spike threshold? 
% 4.	How linear is the summation of postsynaptic potentials? Monaural vs. binaural. What is the role of low-voltage activated K channels in summation? The role of the other dendrite as a sink in the case of monaural stimulation? How linear is the relation between probability of triggering a spike and the average membrane potential?
% 5.	Does summation of EPSPs lead to better phase-locking?
% 6.	Why are spike counts low during ‘worst’ ITD? (What is the role of low-voltage activated K channels?)
% 7.	Are monaural response patterns due to supra-threshold monaural excitatory inputs, or the summation of sub-threshold responses with spontaneous EPSPs originating from inputs of the other ear, as suggested by (Colburn et al., 1990). 
% 8.	What is the effect of postsynaptic spikes on EPSPs and their binaural interactions?
% 9.	Is there a relation between the EPSP time course and CF?
% 10.	Is IPD tuning sharpest at CF (for similar firing rates)?
% 11.	What is the variance in the response to tones?
% 12.	What is the function of cells in the MSO that do not show phase-locking?
% 13.	Is there a rostrocaudal gradient of best ITDs in the gerbil MSO, as previously found for the cat MSO (Yin and Chan, 1990)?
% 
% What is still unknown about the role of inhibition at the gerbil MSO?
% 1.	What is the timing of inhibitory inputs relative to excitatory inputs? Is inhibition well-timed, can it precede excitation, is it long-lasting, etc.
% 2.	Are inhibitory sidebands due to direct inhibitory inputs, or are they created more peripherally?
% 3.	Does inhibition work mostly by shunting excitation, or by an effect on the membrane potential?
% 
% Which questions could be addressed in simulations of in vivo juxtacellular recordings at the MSO?
% 1.	How valid is the assumption that juxtacellular recordings directly provide information about the membrane potential?
% 2.	How would activation of the low-voltage-activated K currents change the relation between intra- and extracellular potentials?
% 3.	Does it matter where the different monaural inputs that are apparently activated by a tone end up in the dendritic tree?
% 4.	Investigate the significance of having brief phase-locked inputs at low frequencies. For example, a brief EPSP will limit activation of LVA K-channels. Would it be possible to get good phase-locking for a slowly rising and decaying current, using for example the iceberg effect or the presence of noise (Gai et al., 2009)?
% 5.	What could be the effect of the LVA K-current at the axon (see e.g. (Goldberg et al., 2008)? Do they dynamically regulate AP threshold?



% ANF responses do not show machine-like multiple onset spikes

% consistency of onset and offset multiple peaks (30 reps FS)
ABFplot({'RG09178','1-21-FS', 4}, 2:30, 1); % 88 Hz
ABFplot({'RG09178','1-21-FS', 5}, 2:30, 1); % 108 Hz
ABFplot({'RG09178','1-21-FS', 6}, 2:30, 1); % 134 Hz
ABFplot({'RG09178','1-21-FS', 7}, 2:30, 1); % 165 Hz
% PCA of these same tonal responses
JLtonePCA('RG09178', '1-21-FS', 4)
JLtonePCA('RG09178', '1-21-FS', 5)
JLtonePCA('RG09178', '1-21-FS', 6)
JLtonePCA('RG09178', '1-21-FS', 7)
JLtonePCA('RG09178', '1-21-FS', 8)


% onset and offset responses: one or multiple peaks. Always present, and slow re typical EPSP.
%  RG10197, cell 8: 'big spikes' /'immense spikes; trivial to trigger them.'
% -------CONTRA tones--------
ABFplot({'RG10197','8-6',1},2:11,[1 4]); % two peaks at onset and offset - always together
ABFplot({'RG10197','8-6',2},2:11,[1 4]); % 2 peaks @onset; 1 @ offset
ABFplot({'RG10197','8-6',3},2:11,[1 4]);
% these onset/offset peaks are not different from phase-locked responses:
ABFplot({'RG10197','8-6',14},2:11,[1 4]); % 134 Hz
ABFplot({'RG10197','8-6',24},2:11,[1 4]); % 267 Hz
ABFplot({'RG10197','8-6',34},2:11,[1 4]); % 267 Hz
ABFplot({'RG10197','8-6',34},2:11,[1 4]); % 536 Hz
% ------IPSI tones (same cell)-------
ABFplot({'RG10197','8-7',1},2:11,[1 3]); 
ABFplot({'RG10197','8-7',2},2:11,[1 3]); 
ABFplot({'RG10197','8-7',3},2:11,[1 3]);
ABFplot({'RG10197','8-7', 4},2:11,[1 3]);
ABFplot({'RG10197','8-7', 14},2:11,[1 3]);
ABFplot({'RG10197','8-7', 24},2:11,[1 3]);
ABFplot({'RG10197','8-7', 34},2:11,[1 3]);
ABFplot({'RG10197','8-7',34},2:11,[1 3]); % 536 Hz

% RG10204, cell 10 'spikes + EPSPs' / 'AP thr is amazingly stable over recordings'
% =====Estimate "typical magnitude" of EPSps in spontaneous activity====
% Early rec: '10-4-BFS'
for ii=1:10, [Yspont, dt] = JLspont(JbB_60dB(ii), 100); STDspont(ii) = std(Yspont); end
% STDspont =
%     0.1906    0.1986    0.1997    0.1981    0.1965    0.1965    0.2033    0.2021    0.1963    0.2022
% Late rec: '10-100-BFS'
for ii=1:10, [Yspont, dt] = JLspont(JbI_70dB(ii), 100); STDspont(ii) = std(Yspont); end
% STDspont =
%     0.2521    0.2575    0.2534    0.2608    0.2540    0.2510    0.2491    0.2518    0.2569    0.2708
% Look at example of spont activity; orde of magn of spont peaks seems OK
JLspontPeaks(JbI_70dB(3));
% now look at the peak groups at onset and offset; looks like ~10 inputs
ABFplot({'RG10204'  '10-8-FS' 1}, 2:11, 1); % 60 dB bin
ABFplot({'RG10204'  '10-8-FS' 2}, 2:11, 1); % 60 dB bin 
ABFplot({'RG10204'  '10-9-FS' 1}, 2:11, 1); % 60 I; no clear multiple peaks

ABFplot({'RG10204'  '10-10-FS' 1}, 2:11, 1); xlim([70 80]) % 60 dB C; here the offsets show "failures"
ABFplot({'RG10204'  '10-9-FS' 1}, 2:11, 1); xlim([70 80]) % 60 dB I; % hardly a response
ABFplot({'RG10204'  '10-9-FS' 2}, 2:11, 1); xlim([70 80]) % 60 dB I;
ABFplot({'RG10204'  '10-9-FS' 3}, 2:11, 1); xlim([70 80]) % 60 dB I; variability

ABFplot({'RG10204'  '10-10-FS' 2}, 2:11, 1); xlim([70 80]) % 60 dB C; no clear "failures", but variabililty
ABFplot({'RG10204'  '10-18-FS' 3}, 2:11, 1); xlim([70 80]) % 50 dB B;

% RG10204  <10-8-FS>    75.4::9000 Hz     60 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-9-FS>    75.4::9000 Hz     60|-60 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-10-FS>   75.4::9000 Hz     -60|60 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-11-FS>   75.4::9000 Hz     -60|40 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-12-FS>   75.4::9000 Hz     40|-40 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-13-FS>   75.4::9000 Hz     40 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-14-FS>   75.4::9000 Hz     30 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-15-FS>   75.4::9000 Hz     -30|30 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-16-FS>   75.4::9000 Hz     30|-30 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-17-FS>   75.4::9000 Hz     -50|50 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-18-FS>   75.4::9000 Hz     50 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-19-FS>   75.4::9000 Hz     50|-50 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-20-FS>   75.4::9000 Hz     70 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-21-FS>   75.4::9000 Hz     70 dB     no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-22-FS>   75.4::9000 Hz     70|-70 dB no mod                     24 x 10 x 70/100 ms    B
% RG10204  <10-23-FS>   75.4::9000 Hz     -70|70 dB no mod                     24 x 10 x 70/100 ms    B

%========STall database===========

%===========APs aligned on rising EPSP slope
JLreadBeats('RG10204', 10);
[D DS E]=readTKabf(JbB_60dB(6)); 
dt = D.dt_ms;
[CH RW AP EPSP] = deristats_tk(readTKaBF(JbB_60dB(6)), 0.05, -10, 10, 1.5);
[APsn, Tap] = getsnips({RW.getsam() dt}, [AP.tEPSPmaxRate], [-1 1]);
plot(Tap, APsn); % catplot version
%-->find best predictor of APs
s = JbB_70dB(5);
[D DS E]=readTKabf(s); 
APcrit = JL_APcrit(s);
dt = D.dt_ms;
[CH RW AP EPSP] = deristats_tk(readTKaBF(s), 0.05, APcrit(1));
Twin = [-0.8 1];
[APsn, Tap] = getsnips({RW.getsam() dt}, denan([AP.tEPSPmaxRate]), Twin);
[dAPsn, dTap] = getsnips({RW.getdsam() dt}, denan([AP.tEPSPmaxRate]), Twin);
[ddAPsn, ddTap] = deal(diff(dAPsn)/dt, dTap(2:end));
f2; h1=subplot(2,1,1); hist([AP.EPSPmaxRate],15); h2=subplot(2,1,2); hist([EPSP.maxUpRate],15); linkaxes([h1 h2], 'x')
% reasonable crit
f2; h1=subplot(2,1,1); hist([AP.VpeakEPSP],15); h2=subplot(2,1,2); hist([EPSP.Vpeak],15); linkaxes([h1 h2], 'x')
[EPSPsn, Tep] = getsnips({RW.getsam() dt}, [EPSP.tmaxUpRate], Twin);
[dEPSPsn, dTep] = getsnips({RW.getdsam() dt}, [EPSP.tmaxUpRate], Twin);
[ddEPSPsn, ddTep] = deal(diff(dEPSPsn)/dt, dTep(2:end));
f3; clf;
xplot(Tap, APsn, 'r'); % catplot version
xplot(Tep, EPSPsn, 'b'); % catplot version
xplot(Tap, APsn(:,rand(1,size(APsn,2))>0.5), 'r'); % catplot version
f4; clf;
xplot(EPSPsn, dEPSPsn, 'b'); 
xplot(APsn, dAPsn, 'r'); 
f5; clf;
xplot(dEPSPsn(2:end,:), ddEPSPsn, 'b'); 
xplot(dAPsn(2:end,:), ddAPsn, 'r'); 
f6;
plot(Tap, mean(APsn, 2), 'r', 'linewidth', 2); % catplot version
xplot(Tap, mean(APsn, 2) + std(APsn, [], 2), 'r'); % catplot version
xplot(Tap, mean(APsn, 2) - std(APsn, [], 2), 'r'); % catplot version
xplot(Tep, mean(EPSPsn, 2), 'b', 'linewidth', 2); % catplot version
xplot(Tep, mean(EPSPsn, 2) + std(EPSPsn, [], 2), 'b'); % catplot version
xplot(Tep, mean(EPSPsn, 2) - std(EPSPsn, [], 2), 'b'); % catplot version
f7;
plot(Tap, mean(dAPsn, 2), 'r', 'linewidth', 2); % catplot version
xplot(Tap, mean(dAPsn, 2) + std(dAPsn, [], 2), 'r'); % catplot version
xplot(Tap, mean(dAPsn, 2) - std(dAPsn, [], 2), 'r'); % catplot version
xplot(Tep, mean(dEPSPsn, 2), 'b', 'linewidth', 2); % catplot version
xplot(Tep, mean(dEPSPsn, 2) + std(dEPSPsn, [], 2), 'b'); % catplot version
xplot(Tep, mean(dEPSPsn, 2) - std(dEPSPsn, [], 2), 'b'); % catplot version
f8;
DT = 0.1; % ms smooth window
[sdAPsn, sdEPSPsn] = deal(smoothen(dAPsn, DT,dt), smoothen(dEPSPsn, DT,dt));
plot(Tap, mean(sdAPsn, 2), 'r', 'linewidth', 2); % catplot version
xplot(Tap, mean(sdAPsn, 2) + std(sdAPsn, [], 2), 'r'); % catplot version
xplot(Tap, mean(sdAPsn, 2) - std(sdAPsn, [], 2), 'r'); % catplot version
xplot(Tep, mean(sdEPSPsn, 2), 'b', 'linewidth', 2); % catplot version
xplot(Tep, mean(sdEPSPsn, 2) + std(sdEPSPsn, [], 2), 'b'); % catplot version
xplot(Tep, mean(sdEPSPsn, 2) - std(sdEPSPsn, [], 2), 'b'); % catplot version


% Illustrating Jeffress ..
JLreadBeats('RG10204', 10);
JlbeatPlot(JbB_60dB(2),3)
JlbeatPlot(JbB_60dB(6),3)  % look for small residuals

% stimulus suppresses the variance below baseline
JLbeatVar(JbI_70dB(2));
% look for cells with quiet baselines


% multiple peak demo
 JLreadBeats('RG10201', 2);
 JLbeatplot(JbI_60dB(2));
 ABFplot(JbI_60dB(2),1,1); % multi peakedness visible in raw data

% variance revisited
JLreadBeats('RG10204', 10);
ifreq=2; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1); 
ifreq=3; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1); 
ifreq=4; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1); 
ifreq=5; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1); 
ifreq=6; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1); 
% at higher freqs, systematic trends are:
%     - little delay between mean & var (typically <100 us)
%     - good corr between mean & var evaluated at stim freq
%     - often also a good corr evaluated at "absent stim freq"
%     - stimulation increases the variance, but often reduces the
%       across-cycle variance, meaning that a lot of the variance is
%       due to a systematic stimulus effect
%     - within a cycle, the low (hyperpolarized) part of the wave
%       corresponds to a reduction of variance below the spont value. This
%       is specially true for monaural stimulation. 
%     - the slope of var-vs-mean plots: ~0.1 mV (900-Hz data below)
%     
ifreq=7; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1);
ifreq=8; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1);
ifreq=9; JLbeatVar(JbI_70dB(ifreq),1); JLbeatVar(JbC_70dB(ifreq),1); JLbeatVar(JbB_70dB(ifreq),1);

% upshot: - stimulation increases overall variance re spont (see non-stim ear)
%         - after subtraction of the variance due to stimulus, the
%           remaining variance is *below* the spontaneous value (see stim ear)
%         - within the stim cycle, the variance increases with membrane
%           potential
%         - peak across-cycle variance lags the peaks of mean waveform
Y = JLbeatVar(JbC_60dB(2),1); % peak var leads, but dip var lags
Y = JLbeatVar(JbI_60dB(2),1); % overall var lag


% other cell
JLreadBeats('RG10201', 2);
ifreq=4; JLbeatPlot(JbB_60dB(ifreq)); 
ifreq=4; Y=JLbeatVar(JbB_60dB(ifreq),1)
JLbeatVar(JbB_60dB(2));
JLbeatPlot(JbB_60dB(2));

% variance w multiple peak
JLreadBeats('RG10201', 2);
JLbeatVar(JbI_60dB(2));
JLbeatVar(JbC_60dB(2));
JLbeatVar(JbB_60dB(2));

JLbeatVar(JbB_60dB(3));
JLbeatPlot(JbB_60dB(3), 3);


% binaural beauty
JLreadBeats('RG10204', 10);
JLbeatPlot(JbB_70dB(5));

% spike predictions, diff cells. Find more examples!
JLreadBeats('RG10196', 4)
ifreq=10; Y = JLbeatVar(JbB_80dB(ifreq),1); JLpredictSpikes(Y, nan);
JLreadBeats('RG10204', 10);
ifreq=3; Y = JLbeatVar(JbB_70dB(ifreq),1); JLpredictSpikes(Y, nan);
ifreq=5; Y = JLbeatVar(JbB_70dB(ifreq),1); JLpredictSpikes(Y, nan);

% yet another one: immense spikes, suggesting recording near axon hillock
JLreadBeats('RG10197', 8); % 70 dB

% yet another one ..
JLreadBeats('RG10198', 1);  % 70 dB

% and another ..
JLreadBeats('RG10204', 12); % 60 dB

% and
JLreadBeats('RG10204', 13);
ifreq=17; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

% and 
JLreadBeats('RG10209', 2); 
ifreq=4; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);i

% and 
JLreadBeats('RG10209', 4);
ifreq=2; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

% and 
JLreadBeats('RG10209', 5);
ifreq=4; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
% JLreadBeats('RG10214', 3); % CHECK AP thr! APs not well detectable
% ifreq=4; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

JLreadBeats('RG10214', 4);
ifreq=2; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10214', 6);
ifreq=4; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10219', 1);
ifreq=9; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10219', 6);
ifreq=2; Y = JLbeatVar(JbB_60dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10223', 2);
ifreq=2; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10223', 3);
ifreq=2; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10223', 4);
ifreq=6; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);

% and
JLreadBeats('RG10223', 6);
ifreq=4; Y = JLbeatVar(JbB_50dB(ifreq),1); JLpredictSpikes(Y, nan);


%========binaural interaction===========
JLreadBeats('RG10196', 4); % 200-800 Hz in 50-Hz steps
ifreq=2; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[-10 0] us
ifreq=3; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[20 -40] us
ifreq=4; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.99]; tau=[20 -30] us
ifreq=5; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 -10] us
ifreq=6; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=7; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)) % rho=[1.00 0.99]; tau=[10 0] us
ifreq=8; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[1.00 0.91]; tau=[20 40] us
ifreq=9; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.78]; tau=[40 380] us
ifreq=10; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.93]; tau=[70 190] us
ifreq=11; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.94]; tau=[80 130] us
ifreq=12; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.96]; tau=[40 110] us
ifreq=13; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[1.00 0.96]; tau=[-10 120] us
ifreq=14; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.97]; tau=[-40 130] us
ifreq=15; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[-50 130] us
ifreq=16; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.98 0.97]; tau=[-70 100] us

JLreadBeats('RG10214', 4); % 200  &up, 100-Hz steps
ifreq=2; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.93 0.96]; tau=[-110 -80] us
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.93 0.98]; tau=[-120 -80] us
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.84 0.99]; tau=[-140 -20] us
ifreq=5; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.67 0.99]; tau=[1720 -120] us
ifreq=6; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.91 0.96]; tau=[-379 -209] us
ifreq=7; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.91 0.98]; tau=[-400 -140] us
ifreq=8; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.74 0.99]; tau=[760 -110] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[200 -130] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq));
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq));
ifreq=12; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq));
ifreq=13; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq));

JLreadBeats('RG10204', 10); % 200 Hz & up; 100-Hz steps
ifreq=2; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 0.99]; tau=[-10 50] us
ifreq=3; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us
ifreq=4; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 20] us
ifreq=5; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[0 20] us
ifreq=6; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 10] us
ifreq=7; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-50 20] us
ifreq=8; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[0.99 1.00]; tau=[-60 10] us
ifreq=9; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[0.99 1.00]; tau=[-40 10] us
ifreq=10; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 10] us
ifreq=11; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us
ifreq=12; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 10] us
ifreq=13; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us


JLreadBeats('RG10204', 12); % 200 Hz an up; 100-Hz steps
ifreq=2; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[0 0] us
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 -30] us
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us
ifreq=5; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 -30] us
ifreq=6; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us
ifreq=7; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us
ifreq=8; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[-40 0] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 0] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 0] us
ifreq=12; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[30 0] us
ifreq=13; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[20 0] us
ifreq=14; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[30 0] us
ifreq=15; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[40 0] us

JLreadBeats('RG10204', 13); % 200 Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.99]; tau=[0 -60] us
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[10 -30] us
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.98]; tau=[20 -50] us
ifreq=5; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[30 0] us
ifreq=6; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[40 -10] us
ifreq=7; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[30 -10] us
ifreq=8; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 10] us
ifreq=12; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=13; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us
ifreq=14; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us
ifreq=15; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us

JLreadBeats('RG10209', 2); %  % 200 Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-30 -40] us
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-20 -40] us
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 -40] us
ifreq=5; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.99]; tau=[-20 -40] us
ifreq=6; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.99]; tau=[-30 -40] us
ifreq=7; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-30 -50] us
ifreq=8; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-30 -60] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-40 -60] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.98 0.99]; tau=[-30 -50] us
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-30 -50] us
ifreq=12; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[-30 -40] us
ifreq=13; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-40 -30] us
ifreq=14; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.99]; tau=[-40 -30] us
ifreq=15; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 0.99]; tau=[-40 -20] us

JLreadBeats('RG10209', 4);   % 200 Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us MLT
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us
ifreq=5; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 -10] us
ifreq=6; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-40 0] us
ifreq=7; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq));  % rho=[0.97 1.00]; tau=[-50 0] us
ifreq=8; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.85 0.99]; tau=[-50 0] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.86 1.00]; tau=[-10 0] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.98 1.00]; tau=[-20 10] us
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[-60 0] us
ifreq=12; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-70 0] us
ifreq=13; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[-80 0] us
ifreq=14; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.98 1.00]; tau=[-50 -10] us 
ifreq=15; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[0 -10] us

JLreadBeats('RG10214', 6);   % 200 Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.99]; tau=[0 0] us
ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 -10] us MLT
ifreq=4; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 0] us
ifreq=5; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[0 0] us
ifreq=6; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 0] us
ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 0] us
ifreq=8; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-30 -10] us
ifreq=9; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 -10] us
ifreq=10; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 -10] us
ifreq=11; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 -20] us
ifreq=12; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 -10] us
ifreq=13; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 -10] us
ifreq=14; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
ifreq=15; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[0 -10] us

JLreadBeats('RG10219', 1);   % 200 Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.99]; tau=[20 10] us
ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[0 0] us
ifreq=4; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[0 -10] us MLT
ifreq=5; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 -10] us
ifreq=6; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us
ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 0] us
ifreq=8; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 0] us
ifreq=9; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.97 1.00]; tau=[-30 0] us
ifreq=10; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.97 1.00]; tau=[-20 -10] us
ifreq=11; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 0] us
ifreq=12; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.98 1.00]; tau=[-30 0] us

JLreadBeats('RG10223', 2);   % 200 Hz & up, 100-Hz steps
% XXXXXXX====reject based on lack of stability (see Jeannette's remark in
% log). Instability is also apparent from systematic DC shifts in V_m.
% ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.85]; tau=[-100 270] us
% ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.96 0.68]; tau=[10 450] us
% ifreq=4; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.90 0.95]; tau=[-220 240] us
% ifreq=5; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=6; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=8; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=9; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=10; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.93]; tau=[10 508] us
% ifreq=11; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=12; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=13; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=14; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=15; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 

JLreadBeats('RG10223', 3); % 200-Hz & up, 100-Hz steps
ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-30 -10] us
ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us
ifreq=4; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 0.99]; tau=[-30 -10] us
ifreq=5; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 0.99]; tau=[-80 -30] us
ifreq=6; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.98]; tau=[-80 -30] us
ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.98 0.99]; tau=[-70 0] us
ifreq=8; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-60 10] us
ifreq=9; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-50 20] us
ifreq=10; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.99]; tau=[-60 30] us
ifreq=11; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq));
ifreq=12; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq));

JLreadBeats('RG10223', 4);  % 200-Hz & up, 100-Hz steps
% XXXXXXX====reject based on lack of stability (see Jeannette's remark in
% log). Instability is also apparent from systematic DC shifts in V_m.
% ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-150 -60] us
% ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-180 -30] us
% ifreq=4; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=5; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=6; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=8; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=9; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=10; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=11; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=12; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.98 1.00]; tau=[-130 0] us
% ifreq=13; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=14; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 
% ifreq=15; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); 

JLreadBeats('RG10223', 6);  % 200-Hz & up, 100-Hz steps
% XXX again huge DC shifts
ifreq=2; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); % rho=[1.00 0.99]; tau=[50 20] us
ifreq=3; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); % rho=[0.99 0.96]; tau=[-50 30] us
ifreq=4; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us
ifreq=5; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=6; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=7; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=8; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); % rho=[0.97 1.00]; tau=[10 -20] us
ifreq=9; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=10; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=11; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=12; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=13; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=14; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 
ifreq=15; JLbinint(JbI_40dB(ifreq), JbC_40dB(ifreq), JbB_40dB(ifreq)); 


% ==== "typical examples" for binaural interaction ===============
JLreadBeats('RG10196', 4); % 200-800 Hz in 50-Hz steps
% multiple peaks & everything
ifreq=2; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[-10 0] us
ifreq=3; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[20 -40] us
ifreq=6; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us
% but an asymmetric delay at still higher freqs
ifreq=15; JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq)); % rho=[0.99 0.98]; tau=[-50 130] us


JLreadBeats('RG10204', 10); % 200 Hz & up; 100-Hz steps
ifreq=3; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[0 10] us
ifreq=6; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 10] us
ifreq=11; JLbinint(JbI_70dB(ifreq), JbC_70dB(ifreq), JbB_70dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us

JLreadBeats('RG10204', 12); % 200 Hz an up; 100-Hz steps
ifreq=4; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us
ifreq=11; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[0 0] us

JLreadBeats('RG10204', 13); % 200 Hz & up, 100-Hz steps
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[10 -30] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[10 0] us

JLreadBeats('RG10209', 2); %  % 200 Hz & up, 100-Hz steps
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 0.99]; tau=[-20 -40] us
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.98 0.99]; tau=[-30 -50] us

JLreadBeats('RG10209', 4);   % 200 Hz & up, 100-Hz steps
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[1.00 1.00]; tau=[-10 0] us MLT
ifreq=10; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.98 1.00]; tau=[-20 10] us

JLreadBeats('RG10214', 4); % 200  &up, 100-Hz steps
% multi peaks, good corr but quite a delay
ifreq=3; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.93 0.98]; tau=[-120 -80] us
ifreq=9; JLbinint(JbI_60dB(ifreq), JbC_60dB(ifreq), JbB_60dB(ifreq)); % rho=[0.99 1.00]; tau=[200 -130] us

JLreadBeats('RG10214', 6);   % 200 Hz & up, 100-Hz steps
ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 1.00]; tau=[-20 -10] us MLT
ifreq=7; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 0] us

JLreadBeats('RG10219', 1);   % 200 Hz & up, 100-Hz steps
% also DC shift, good corr, mult peaks
ifreq=2; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[0.99 0.99]; tau=[20 10] us

JLreadBeats('RG10223', 3); % 200-Hz & up, 100-Hz steps
% DC shift, but good correlation & multiple peaks
ifreq=3; JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq)); % rho=[1.00 1.00]; tau=[-20 -10] us


%=========DATA SELECTION BINAURAL PERIODICTY ANALYSIS==============
JLreadBeats('RG10214', 6);
% single peaks
ifreq=4; Y= JLbeatVar(JbB_50dB(ifreq),1); P = JLpredictSpikes(Y,nan,1); B = JLbinint(JbI_50dB(ifreq), JbC_50dB(ifreq), JbB_50dB(ifreq));
% mlutiple peaks
ifreq=2; Y= JLbeatVar(JbB_80dB(ifreq),1); P = JLpredictSpikes(Y,nan,1); B = JLbinint(JbI_80dB(ifreq), JbC_80dB(ifreq), JbB_80dB(ifreq));






