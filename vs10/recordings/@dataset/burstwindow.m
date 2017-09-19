function Tw=burstwindow(D, iCond, Chan);
% Dataset/burstwindow - time window of stimulus burst
%   Tw=burstwindow(D, iCond) returns a timewindow Tw = [Tstart Tend]
%   indicating the start and end (in ms) of the stimulus "burst" of
%   stimulus condition iCond. Tw is determined from the following 
%   GenericStimParams:
%      OnsetDelay
%        BurstDur.
%   If iCond is an array, Tw is a Nx2 array with Tw(k,:) being the time
%   window of the k-th stimulus condition. Default value for k is 0, mening
%   "all conditions";
%   If the two DAC channels are both active and have have different
%   stimulus timings, the "union" of the two bursts is taken.
%
%   burstwindow(D, iCond, 'L') and burstwindow(D, iCond, 'R') restrict the
%   window to the stimulus segemnt of requested DAC channel. If a 
%   non-active DAC channel is requested, Tw will be filled with NaNs.
%
%   See also GenericStimParams.

[iCond, Chan] = arginDefaults('iCond, Chan', 0, 'B');

GSP = GenericStimparams(D);
DAC = GSP.DAC; 
Ncond = GSP.Ncond;
if isequal(0, iCond),
    iCond = 1:Ncond;
end
OnsetDelay = GSP.OnsetDelay(iCond,:);
BurstDur = GSP.BurstDur(iCond,:);
OnsetDelay = SameSize(OnsetDelay, BurstDur);

% set non-active DAC values to NaN
if isequal('L', DAC(1)),
    OnsetDelay(:,2) = nan;
    BurstDur(:,2) = nan;
elseif isequal('R', DAC(1)),
    OnsetDelay(:,1) = nan;
    BurstDur(:,1) = nan;
end

% select requested channels
OnsetDelay = channelSelect(Chan, OnsetDelay);
BurstDur = channelSelect(Chan, BurstDur);
T0 = OnsetDelay;
T1 = T0 + BurstDur;

T0 = min(T0,[], 2); % start of burst: min across channels
T1 = max(T1,[], 2); % end of burst: max across channels

Tw = [T0, T1];











