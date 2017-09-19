function S = JLbeat(ExpID, RecID, icond, Nxmean, AP, doPlot, seriesID);
% JLbeat - analysis of binaural beat from JL's recordings
%   S = JLbeat(ExpID, RecID, icond, Nxmean, APparam)
%      APparam = [minSlopeAcceptance -Tbefore Tafter minSlopeCutout]

if nargin<2, RecID=''; end
if nargin<3, icond=nan; end
if nargin<4, Nxmean=1; end
if nargin<5, AP = -10; end % -10 V/s ; 
if nargin<6, doPlot=(nargout<1); end
if nargin<7, seriesID = '??'; end

if numel(AP)<3,
    AP = [AP(1) -1 1.4]; %truncation 0.35 ms before peak to 1.00 ms after peak
end
if numel(AP)<4, % cut-out limit equal to acceptance limit
    AP = [AP AP(1)];
end
ThrAPdownRate = AP([1 4]); % first thr is acceptance as AP; second (smaller abs) thr for being cut out
APwindow = AP(2:3); % AP truncation window

MyFlag('LastRecordingAnalyzed', {ExpID, RecID, icond});
MaxOrd = 9;

AvMethod = 'mean'; % default

Ttrans = 100; % ms first & last part during stim ignored

[D, DS, L]=readTKABF(ExpID, RecID, icond);
[RecID, icond, icell] = deal(D.RecID, D.icond, DS.icell);
[ABFname, ABFdir] = deal(L.Recs{icond}, L.Dir);
if ~isempty(strfind(lower(ABFname), 'fake')),
    error(['Missing data as indicated by data file name ''' ABFname '''.']);
end
if ~isequal('BFS', DS.stimtype),
    error([DS.title ' is not a BFS dataset.']);
end
SPL1 = DS.SPL(1);
SPL2 = DS.SPL(end);
TTT = [num2str(icond) ': ' num2str(DS.xval(icond)) ' ' DS.x.Unit '  ' DS.filename '' DS.info ' <' ABFname '> ' ];


Fsam = 0.999993*D.AD(1).Fsam_kHz; % kHz (% correct for different clock rates of A/D and D/A)
dt = 1/Fsam; % sample period in ms
preDur = 500; % ms pre-stim by convention
burstDur = DS.burstdur(1); % ms tone duration
postDur = DS.repdur-burstDur; % ms post-tone dur
Freq1 = DS.fcar(icond,1); Period1 = 1e3/Freq1; % cycle in ms
Freq2 = DS.fcar(icond,2); Period2 = 1e3/Freq2; % cycle in ms
Nsam1 = round(1e3*Fsam/Freq1);
Nsam2 = round(1e3*Fsam/Freq2);
Tbeat = 1e3/DS.Fbeat(icond); % ms beat dur
Tcar = 1e3/Freq1;
NcarCycleBeat = round(Tbeat/Tcar); % make sure that beat cycle ...
Tbeat = NcarCycleBeat*Tcar; % ... contains integer # car cycles


stim1 = D.AD(3).samples; % ipsi stimulus waveform
stim2 = D.AD(4).samples; % contra stimulus waveform
rec = D.AD(1).samples; % recording
%rec = local_no_50Hz(rec, dt);

% timestamp & truncate APs
%[rec2, SPT] = APtruncate2(rec, ThrAPdownRate, dt, [-0.35 1.00]);
[rec, SPTraw, Tsnip, Snip, APslope] = APtruncate2(rec, ThrAPdownRate, dt, APwindow);

Time = timeaxis(stim1, dt);
tstart = preDur + Ttrans; % ms start of analysis
tend = preDur + burstDur; % ms end of analysis
Nbeats = floor((tend-tstart)/Tbeat); % # complete beats
NsamBeat = round(Tbeat*Fsam); % # samples in one beat
isam = round(tstart/dt) + (1:(Nbeats*NsamBeat)); % sample range
% apply analysis window
Y = [stim1, stim2, rec];
Y = Y(isam,:);

[dum, dum, AC] = LoopMean(Y(:,3),NsamBeat);
Y = LoopMean(Y,NsamBeat);

% extra averaging of adjacent carrier cycles
DTcar = 1e3/mean(DS.fcar(icond,:)); % approx carrier cycle
nsamCarCycle = round(DTcar/dt); % # samples corresponding to approx carrier cycle
Z = 0;
for ishift=(1:Nxmean)-floor((Nxmean+1)/2),
    Z = Z + circshift(Y,ishift*nsamCarCycle);
end
Y = Z;
Y = Z/Nxmean;
% circular shift of Y to move beat max to middle of trace 
Tshift = Ttrans-Tbeat/2;
NsamShift = round(Tshift/dt);
Y = circshift(Y,[NsamShift 0]);

% apply same analysis window and start-time convention to the APs returned 
% by deristats above.
SPT = SPTraw(SPTraw>tstart)-tstart; % delete pre-stim & trans intervals
SPT = SPT(SPT<=Nbeats*Tbeat); % the same range as applied to Y
SPT = mod(SPT+Tshift, Tbeat);


% attempt at decomposition in ipsi and contra
BeatDur = NsamBeat*dt;
df = 1e3/BeatDur; % Hz freq spacing
[Sb, Fdpb] = DPfreqs([Freq1 Freq2],MaxOrd);
R0 = Y(:,3); R0=R0-mean(R0);
Stim1 = Y(:,1); Stim2 = Y(:,2);
Sfactor = std(2*R0)/sqrt(var(Stim1)+var(Stim2)); % scale factor giving StimX ~same magnitude as R0
Stim1 = Sfactor*Stim1; Stim2 = Sfactor*Stim2; 
Spr = fft(R0);
[Z1, Z2] = deal(0*Spr); 
NmaxHarm = 0.4*1e3*Fsam/max([Freq1 Freq2]);
ih1 = 1+round(Freq1*(1:NmaxHarm)/df);
ih2 = 1+round(Freq2*(1:NmaxHarm)/df);
Z1(ih1) = Spr(ih1);
Z2(ih2) = Spr(ih2);
idpb = 1+round([Sb.freq]/df); % all DPs from freq1 and freq2 together
Z3 = Z1+Z2; Z3(idpb) = Spr(idpb); % add binaural interactions
Y1 = 2*real(ifft(Z1));
Y2 = 2*real(ifft(Z2));
Y3 = 2*real(ifft(Z3));
Yc = Y1+Y2; % prediction by adding ipsi & contra
Yoffset = 2.5*std(R0); % vertical plotting offset to avoid overlap

% do some bookkeeping, and save params & data in struct
if isequal(SPL1,SPL2), StimType = 'B'; % binaural
elseif SPL1<=0, StimType = 'C'; % contra only
elseif SPL2<=0, StimType = 'I'; % ipsi only
else, StimType = 'X'; % special
end

Fbeat = abs(Freq1-Freq2);
AnaDur = 1e3*Nbeats/Fbeat; % duration of analysis window
APthrSlope = ThrAPdownRate(1);
CutoutThrSlope = ThrAPdownRate(end);
APwindow_start = APwindow(1);
APwindow_end = APwindow(2);
iGerbil = str2double(ExpID(5:7));
UniqueCellIndex = 100*iGerbil+icell; % unique identifier of the cell
UniqueSeriesIndex = 100*UniqueCellIndex - DS.iseq; % unique identifier of the binbeat sweep
UniqueRecordingIndex = 100*UniqueSeriesIndex + icond; % unique identifier of the single recording
[JLcomment, MHcomment] = deal('');
S = CollectInStruct(ExpID, RecID, seriesID, iGerbil, icond, icell, '-', ...
    UniqueCellIndex, UniqueSeriesIndex, UniqueRecordingIndex, ABFname, ABFdir, JLcomment, MHcomment, '-', ...
    Nxmean, APthrSlope, APwindow_start, APwindow_end, CutoutThrSlope, '-', ...
    StimType, Freq1, Freq2, Fbeat, burstDur, SPL1, SPL2, '-', ...
    Fsam, dt, Nbeats, AnaDur, Ttrans, MaxOrd, TTT, Period1, Period2, '-', ...
    Stim1, Stim2, '-', ...
    R0, Y1, Y2, Y3, Yc, df, Yoffset, Sb, Fdpb, SPTraw, SPT, AC, Tshift, '-', ...
    Tsnip, Snip, APslope, '-');

if doPlot,
    JLbeatPlot(S);
end

  
%======================================
function rec = local_no_50Hz(rec, dt);
Dur = dt*numel(rec);
if (rem(Dur,20)<0.01*Dur) || (20-rem(Dur,20)<0.01*Dur),
    df = 1e3/Dur; % spacing in Hz
    Zsp = fft(rec);
    i50 = 1+round(50/df); % index of 50-Hz comp
    Zsp([i50 end+2-i50])=0; % end+2-i50 is the neg 50-Hz comp
    rec = real(ifft(Zsp));
end


