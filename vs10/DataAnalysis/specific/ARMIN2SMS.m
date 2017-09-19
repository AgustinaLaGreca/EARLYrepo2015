function SMS=ARMIN2SMS(pp);
% ARMIN2SMS - convert ARMIN menu parameters to SMS XXX

stimname = 'armin';

prepareARMINstim(pp); 

global AKbuffer

StimParam = pp;

% -------------------presentation--------------
Nsubseq = length(pp.flipfreqs);
Nrep = pp.reps;
chan = channelChar(pp.active);
if length(pp.SPL)==1, pp.SPL = [1 1]*pp.SPL; end;
SPL = pp.SPL;
% collect in struct
PRES = CollectInStruct(Nrep, Nsubseq, chan, SPL);

% ---------PRP info------------
playOrder = CreatePlayorder(Nsubseq, pp.order);
PlotMode = pp.flipstepunit;
plotInfo = CreatePlotInfo('Flip freq (Hz)', pp.flipfreqs(playOrder), PlotMode, NaN, 'BurstOnly');
PRP = CreatePRPinfo(plotInfo, playOrder);

% ----------global info-----------
global CALIB SGSR
var1v = repmat(pp.flipfreqs,1,2);
var2v = 0*var1v;
GI = CreateGlobalInfo('', CALIB.ERCfile, var1v, var2v, 1);
GI = rmfield(GI,'cmenu');
GI.stimName = stimname;
GI.nonPDP11 = 1;
GI.variedParam = 'flipFreq';
GI.StimParams = pp;
GI.SamFreqs = SGSR.samFreqs;
GI.maxSampleRatio = SGSR.maxSampleRatio;
GlobalInfo = GI;
createdby = mfilename;

StimCat = 'armin';
SMS = CollectInStruct(StimCat, StimParam, Nsubseq, PRES, GlobalInfo, PRP, createdby);
