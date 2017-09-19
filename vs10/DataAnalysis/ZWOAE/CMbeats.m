function [QQ, qq] = CMbeats(expID, recID, ibeats)
% [QQ, qq] = CMbeats(expID, recID, ibeats)

% Some examples for this function.
%
%                               F1      F2   Fbeat  L1=L2
% expID = 190; recID = 1358;    %332   282    50    95
% expID = 190; recID = 1390;    %332   282    50    95
% expID = 190; recID = 1423;    %468   418    50    95
% expID = 190; recID = 1454;    %468   530    62    95
% expID = 190; recID = 1513;    %468   418    50    80

extraplotFlag = 0;

if nargin<1, expID=[]; end
if nargin<3, ibeats=1; end

if isempty(expID),
    load('190_1423'); % this load qq and QQ
else,
    QQ = rawCMimport(expID, recID,0); %get the data, no plotting
    qq = ZWOAEimport(expID, recID); %get the stim params
end
Resp = LoopMean(QQ.signal, numel(QQ.signal)/2); %CM signal; single period
Resp =  Resp - mean(QQ.DC_mean); %correct for DC pre- and post-recording
Resp = Resp/max(Resp); %normalize
dt = QQ.dt;
df = 1/(dt*numel(Resp));
% reconstruct stimulus by comb filtering the response
Cspec = fft(Resp); %spectrum of signal (complex)
i1 = 1+round(qq.Fzwuis/df); % indices of stim components ...
i2 = 1+round(qq.Fsingle/df); % ... in the spectrum

StimSpec = Cspec*0; % zero complex spectrum 
StimSpec([i1 i2]) = Cspec([i1 i2]); % only select stim tones from resp spectrum
Sifft = ifft(StimSpec);
Stim = 2*real(Sifft); %stimulus time signal
StimEnv = 2*abs(Sifft); %stimulus envelope

% divide stimulus and response into "carrier cycles", i.e. chunks between
% consecutive stim minima
imin = localmax(1,-Stim); % indices of local minima of stimuli
imaxEnv = localmax(1,StimEnv); % indices of local maxima of stim envelope
% apply circular shift that makes stim & response start at stim minimum
% near envelope maximum
istart = find(imin>imaxEnv(1),1,'first');
ishift = imin(istart);
[Stim, Resp, StimEnv] = deal(circshift(Stim,1-ishift), circshift(Resp,1-ishift), circshift(StimEnv,1-ishift));
istartCyc = localmax(1, -Stim); % indices of start of carrier cycles
iendCycle = [istartCyc(2:end)-1; numel(Stim)]; % ditto ends
NcarCycle = numel(istartCyc); % # carrier cycles
%======================
%======================

% take time derivative of stim and response after gentle smoothing
Nwindow = 5; %window length in smoothen.m
smResp = smoothen(Resp, -Nwindow); % smooth the response while respecting its periodic character
dResp = diff(smResp([1:end 1]))/dt; % take time derivative
dStim = diff(Stim([1:end 1]))/dt; % time derivative of stimulus

% plot Stim & Resp in proper time
fh1 = figure; set(gcf,'units', 'normalized', 'position', [0.0102 0.256 0.983 0.615]);
set(fh1,'units', 'normalized', 'position', [0.00859 0.075 0.983 0.615])
TT = timeaxis(Stim, dt); Tmax = max(TT);
ah1 = subplot(2,1,1);
for icyc=1:NcarCycle,
    iend = iendCycle(icyc) + double(icyc<NcarCycle);
    ir = istartCyc(icyc):iend; % index range of current cycle
    xplot(TT(ir), Stim(ir), [ploco(icyc) ':']);
    xplot(TT(ir), Resp(ir), ploco(icyc));
    if icyc==1, xlim([0 Tmax]); end
end
TracePlotInterface(fh1);
% lissajous
%fh2 = figure; set(fh2,'units', 'normalized', 'position', [0.00313 0.0725 0.991 0.341])
ah2 = subplot(2,1,2);
for icyc=1:NcarCycle,
    istart = istartCyc(icyc);
    iend = iendCycle(icyc) + double(icyc<NcarCycle);
    ir = istart:iend; % index range of current cycle
    xoffset = dt*mean([istart iend]);
    xplot(xoffset+Stim(ir), Resp(ir), ploco(icyc));
    if icyc==1, xlim([0 Tmax]); end
end
ah2 = gca;

fh3 = figure; 
set(fh3,'units', 'normalized', 'position', [0.0102 0.639 0.983 0.233]);
TT = timeaxis(Stim, dt); Tmax = max(TT);
ah3 = gca;
% Gain stuff
for icyc=1:NcarCycle,
    istart = istartCyc(icyc);
    iend = iendCycle(icyc) + 1; iend = min(iend,numel(dStim));
    ir = istart:iend; % index range of current cycle
    dR = dResp(ir); 
    dS = dStim(ir); 
    Time = TT(ir)+dt/2; % 1/2 sample shift because of diff
    Gain = dR./abs(dS);
    % determine max gain while excluding shallow part of stimulus
    qok = (dS>0.5*max(dS)) | (dS<0.5*min(dS));
    Gain = Gain(qok);
    Time = Time(qok);
%     xplot(Time, 3*Gain, [ploco(icyc) '.-']);
%     xplot(Time, dR(qok), [ploco(icyc) '-']);
%     xplot(Time, dS(qok), [ploco(icyc) ':']);
    [maxGain(icyc), imax] = max(Gain)
    [minGain(icyc), imin] = min(Gain)
    TmaxGain(icyc) = Time(imax);
    TminGain(icyc) = Time(imin);
end
% xplot(TT, 0*TT, 'k-');


% % plot time derivatives of Stim & Resp in proper time
% for icyc=1:NcarCycle,
%     iend = iendCycle(icyc) + double(icyc<NcarCycle);
%     ir = istartCyc(icyc):iend; % index range of current cycle
%     xplot(TT(ir), dStim(ir), [ploco(icyc) ':']);
%     xplot(TT(ir), dResp(ir), ploco(icyc));
%     if icyc==1, xlim([0 Tmax]); end
% end
% TracePlotInterface(fh3);


linkaxes([ah1 ah2 ah3], 'x');
xlim(ah1, [0 20])
delete(fh3);


% find out how gain and stimulus envelope are related
fh2 = figure; % time plot of gains & envelope, including lead/lag info
set(gcf,'units', 'normalized', 'position', [0.032 0.474 0.573 0.419]);
Interp_maxGain = interp1(TmaxGain, maxGain, TT);
Interp_minGain = interp1(TminGain, minGain, TT);
qok = ~isnan(Interp_maxGain) & ~isnan(Interp_minGain);
TT = TT(qok);
Interp_maxGain = Interp_maxGain(qok);
Interp_minGain = Interp_minGain(qok);
plot(TT, Interp_maxGain, 'b')
xplot(TT, -Interp_minGain, 'r')
StimEnv = StimEnv(qok);
StimEnv = StimEnv/max(StimEnv);
xplot(TT, StimEnv, 'g');
meanBeatDur = dt*mean(diff(imaxEnv)); % approx beat dur
maxLag = round(meanBeatDur/4/dt); % # mag lag [samples] ~ 1/4 beat cycle
[MaxCp, ilagp] = maxcorr(-Interp_maxGain, StimEnv, maxLag);
[MaxCn, ilagn] = maxcorr(Interp_minGain, StimEnv, maxLag);
Lagp = dt*ilagp; Lagn = dt*ilagn;
legend({['Up-gain lags env by ' num2str(Lagp) ' ms'], ...
    ['Down-gain lags env by ' num2str(Lagn) ' ms']})

% Envelope vs Gain Lissajous
figure
plot(StimEnv, Interp_maxGain);
xplot(StimEnv, -Interp_minGain, 'r');
set(gcf,'units', 'normalized', 'position', [0.637 0.435 0.341 0.451])
xlabel('Stimulus envelope');
ylabel('Gain at steepest I/O portion');
