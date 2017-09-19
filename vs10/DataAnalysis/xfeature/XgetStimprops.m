function M = XgetStimprops(ExpID, RecID, icond);
% XgetStimprops - stimulus properties in a struct
%   SP = XgetStimprops(ExpID, RecID, icond) returns struct SP containing
%   stimulus properties.
%
%   XgetStimprops(S) is the same as XgetStimprops(S.ExpID, S.RecID, S.icond)
%   
%   See also XcycleHist.

if nargin==1 && isstruct(ExpID),
    D = ExpID;
    [ExpID, RecID, icond] = deal(D.ExpID, D.RecID, D.icond);
end

[D, DS] = readTKABF(ExpID, RecID, icond);
M.DS = DS;
M.icond = icond;
M.StimType = DS.StimType;
[dum M.StimInfo] = DSinfoString(DS);
M.Nsweep = D.Nsweep;
M.SweepDur = D.sweepDur_ms;
M.OrigrepDur = DS.repdur;
M.burstDur = local_reduceifyoucan(DS.burstDur, DS.dachan, M.icond);
M.SPL = local_reduceifyoucan(getFieldOrDefault(DS, 'SPL', nan), DS.dachan, icond);
M.Fcar = local_reduceifyoucan(DS.carfreq, DS.dachan, icond);
M.Xname = DS.X.Name;
M.delay = local_reduceifyoucan(getFieldOrDefault(DS, 'delay', 0), DS.dachan, icond);
M.Xval = local_reduceifyoucan(DS.Xval, DS.dachan, icond);
M.Xunit = DS.Xunit;
if isnan(M.SPL),
    if isequal('SPL', M.StimType),
        M.SPL = M.Xval;
    else,
        warning('Cannot determine SPL of stimulus');
    end
end
M.onset = (1:D.Nsweep-1)*M.SweepDur;
M.offset = M.onset + M.burstDur;


function X = local_reduceifyoucan(X, dachan, icond);
% reduces stumulus param to value of current condition & channel, if applicable
if dachan>0 && size(X,2)>1, % monaural stim, reduce X to that of active DA chan
    X = X(:,dachan);
end
if nargin>2 && size(X,1)>1, % X was varied; reduce to that of current condition
    X = X(icond,:);
end
% if X exists of 2 identical colms (2 channels), reduce to one
if isequal(X(:,1), X(:,end)),
    X = X(:,1);
end





