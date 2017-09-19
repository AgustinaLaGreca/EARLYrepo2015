function SMS=FM2SMS(idfSeq,calib);

% FM2SMS - convert FM menu parameters to SMS XXX
if nargin<2, calib = 'none'; end;

stimname = IDFstimname(idfSeq.stimcntrl.stimtype);

[chan sc s1 s2 scmn chanNum order] = idfstimfields(idfSeq);

%------------frequencies---------------
% collect anchor frequencies 
Flow = [s1.fmcarrlo s2.fmcarrlo];
Fhigh = [s1.fmcarrhi s2.fmcarrhi];

% compute safe samplefreqs and filter indices
[samFreq, filterIndex] = safeSampleFreq(max(max(Flow,Fhigh))); % max over channels
% collect info in struct var
upDur = [s1.sweepup s2.sweepup];
holdDur = [s1.sweephold s2.sweephold];
downDur = [s1.sweepdown s2.sweepdown];
riseDur = [s1.rise s2.rise];
fallDur = [s1.fall s2.fall];
delay = [s1.delay s2.delay]; % per-channel, POSITIVE delays
repDur= sc.interval; 
SWEEP = CollectInStruct(chan, Flow, Fhigh, ...
   repDur, upDur, holdDur, downDur, ...
   riseDur, fallDur, delay, samFreq, filterIndex);

% -------------------presentation--------------
Nrep = sc.repcount;
SPL = [s1.spl s2.spl];
Nsubseq = 1; % by definition for FM stim
% collect in struct
PRES = CollectInStruct(Nrep, SPL, Nsubseq, chan);

% ---------PRP info------------
playOrder = CreatePlayorder(Nsubseq, order);
xdummy = [0 0];
plotInfo = CreatePlotInfo('Nope (AU)', ...
   xdummy(:,1), 'lin', NaN,'BurstOnly');
PRP = CreatePRPinfo(plotInfo, playOrder);

% ----------global info-----------
GlobalInfo = CreateGlobalInfo(stimname, calib, xdummy);

StimCat = 'linfm';
createdby = mfilename;
SMS = CollectInStruct(StimCat, Nsubseq, SWEEP, PRES, PRP, GlobalInfo, createdby);
