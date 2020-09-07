function [SPL, S] = ILD2SPL(SPL,ILD, EXP);
% ITD2delay - convert ITD to per-channel delay
%    ITD2delay(ITD, EXP) returns a 2-element array [dleft, dright] or scalar
%    specifying the delays that realize the specified ITD. The per-channel
%    delays dleft and dright are non-negative delays, one of which equals 
%    zero. EXP is the experiment, whose Preference field is needed to convert 
%    the ipsi/contra-based ITD spec to the left/right-based delay spec.
%
%    A scalar is returned only when ActiveAudioChan(EXP) is not Both, i.e., 
%    when only a single DA channel is available. A two-element array is 
%    returned whenever both DA channels are available, even when only one 
%    is active. The scalar value of the delay is evaluated by selecting 
%    first the dual-channel value and then selecting the available channel. 
%    This is done to ensure complete consistency between true binaural 
%    measurements and measurements in which monaural responses are used to 
%    construct a quasi binaural response, e.g., through crosscorrelation.
%
%    When ITD is a column vector, ITD2delay returns a Nx2 matrix or column
%    vector.
%
%    ITD2delay is a helper function for waveform calculators like toneStim
%    and noiseStim.
%
%    [Delay, S] = ITD2delay(ITD, EXP) also returns a struct with fields L
%    and R, containing the delays for channels L and R, if present.
%
%    See also IPD2phaseShift, IFD2freqShift, toneStim, noiseStim,
%    makeStimFS.

% first compute Delay assuming ILD>0 means LeftLeading ..
SPLOrg = SPL;
SPL(:,1) = SPLOrg+(ILD/2);
SPL(:,2) = SPLOrg-(ILD/2);

% SPL = [0*SPL SPL]-min(0*SPL,SPL)*[1 1];

% .. but ITD might really mean RightLeading 
IpsiIsLeft = isequal('Left', EXP.RecordingSide);
switch EXP.ITDconvention
    case 'IpsiLeading', % swap if ipsi=right 
        DoSwap = ~IpsiIsLeft;
    case 'ContraLeading', % swap if ipsi=left 
        DoSwap = IpsiIsLeft;
   case 'LeftLeading',
        DoSwap = false;
    case 'RightLeading',
        DoSwap = true;
    otherwise,
        error('Invalid ITDconvention value in Experiment.');
end 
if DoSwap,
    SPL = fliplr(SPL);
end
% reduce to Delay scalar if only a single DA channel is available
if isequal('Left', EXP.AudioChannelsUsed), 
    SPL = SPL(:,1); % Left
    S.L = SPL;
elseif isequal('Right', EXP.AudioChannelsUsed), 
    SPL = SPL(:,2); % Right
    S.R = SPL;
else,
    S.L = SPL(:,1);
    S.R = SPL(:,2);
end







