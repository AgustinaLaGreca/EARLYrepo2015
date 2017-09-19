function [WS, effdur] = anwin(S, AW, irep, isub)
% anwin - Apply analysis window.
%    anwin(X, [W0 W1]) returns those elements x of X which are within the analysis
%    window: W0<=x<W1.
%    anwin(X, W) is the same as anwin(A, [0 W]).
%
%    When X is a cell matrix, the window is applied to the individual cells.
%
%    anwin(DS, W), where DS is a dataset variable, returns the spike times
%    of DS restricted to the analysis window. That is, 
%    anwin(DS, W) = anwin(SPT, W), where SPT = DS.SPT. In this case, the
%    anwin(DS, []) or anwin(DS) use the default window [0 min(DS.burstdur)]
%
%    anwin(X,W,irep) or anwin(DS,W,irep) only selects elements of
%    the "repetitions" irep, that is, if SPT = anwin(X,W), then
%    anwin(X,W,irep) equals SPT(:,irep). Of course, this only works 
%    when SPT is a cell-matrix. 
%    anwin(X,W, ':') and anwin(X,W, 0) are both the same as anwin(X,W)
%
%    anwin(X,W,irep,isub) or anwin(DS,W,irep,isub) also selects one
%    or more subsequences (conditions) of SPT, that is SPT is
%    restricted to SPT(isub,irep). Default is isub=0, meaning all
%    subsequences. As before, W=[] uses the default window 
%    [0 min(DS.burstdur)].
%
%    [S, effDur] = anwin(S,AW,..) also returns the effective duration of the
%    spike trains, i.e. effDur = diff(AW)
%
%    See also DATASET.

if nargin<2
    AW = [];
end
if nargin<3
    irep = ':';
end
if nargin<4
    isub = ':';
end
if isequal(0,irep)
    irep = ':'; % see help text
end
if isequal(0,isub)
    isub = ':'; % see help text
end

if iscell(S)
   % select subs & reps and use recursive call to anwin
   S = S(isub,irep);
   Sz = size(S);
   WS = cell(Sz);
   for isub=1:Sz(1)
      for irep=1:Sz(2)
         [WS{isub,irep}, effdur] = anwin(S{isub,irep}, AW);
      end
   end
elseif isa(S,'dataset') % get default AW if needed, extract spiketimes & use recursive call
   if  isempty(AW)
       AW = [0, min(S.Stim.BurstDur)];
   end
   [WS, effdur] = anwin(spiketimes(S), AW, irep, isub);
elseif isnumeric(S)
   if ~isequal(':', irep) || ~isequal(':', isub)
      error(['Selections of repetitions & subsequences cannot be' ...
          'performed on numerical format of spike times.']);  
   end
   if length(AW)<2
       AW = [0 AW];
   end
   WS = S((S>=AW(1))&(S<AW(2)));
   effdur = diff(AW);
else
    error(['Invalid data class "' class(S) '" for anwin.']);
end
