function [WS, effdur] = AnWin(S, AW, irep, isub);
% AnWin - Apply analysis window.
%    AnWin(X, [W0 W1]) returns those elements X which are within the analysis
%    window: W0<=x<W1.
%    AnWin(X, W) is the same as AnWin(A, [0 W]).
%
%    When X is a cell matrix, the window is applied to the individual cells.
%
%    AnWin(DS, W), where DS is a SGSRdataset variable, returns the spike times
%    of DS restricted to the analysis window. That is, 
%    AnWin(DS, W) = AnWin(DS.SPT, W). In this case, 
%    AnWin(DS, []) or AnWin(DS) uses the default window [0 min(DS.burstdur)]
%
%    AnWin(X,W,irep) or AnWin(DS,W,irep) only selects elements of
%    the "repetitions" irep, that is, if SPT = AnWin(X,W), then
%    AnWin(X,W,irep) equals SPT(:,irep). Of course, this only works 
%    when SPT is a cell-matrix. 
%    AnWin(X,W, ':') and AnWin(X,W, 0) are both the same as AnWin(X,W)
%
%    AnWin(SPT, W ,'pool'), where SPT is a [Nsub x Nrep] cell matrix, pools 
%    the spike times across the "repetitions" of SPT.
%
%    AnWin(DS, W, 'pool'), is the same as AnWin(DS.SPT, W, 'pool').
%
%    AnWin(X,W,irep,isub) or AnWin(DS,W,irep,isub) also selects one
%    or more subsequences (conditions) of SPT, that is SPT is
%    restricted to SPT(isub,irep). Default is isub=0, meaning all
%    subsequences. As before, W=[] uses the default window 
%    [0 min(DS.burstdur)].
%
%    [S, effDur] = AnWin(S,AW,..) also returns the effective duration of the
%    spike trains, i.e. effDur = diff(AW)
%
%    See also SGSRdataset.

if nargin<2, AW = []; end
if nargin<3, irep = ':'; end
if nargin<4, isub = ':'; end

if isequal(0,irep), irep = ':'; end % see help text
if isequal(0,isub), isub = ':'; end % see help text

if iscell(S),
   % select subs & reps and use recursive call to AnWin
   if isequal('pool', irep),
       S = localPool(S);
       irep=1;
   end
   S = S(isub,irep);
   [Nsub, Nrep] = size(S);
   WS = cell(Nsub,Nrep);
   for isub=1:Nsub,
      for irep=1:Nrep,
         [WS{isub,irep}, effdur] = AnWin(S{isub,irep}, AW);
      end
   end
elseif isa(S,'sgsrdataset'), % get default AW if needed, extract spiketimes, & use recursive call
   if  isempty(AW), AW = [0, min(S.burstdur)]; end
   [WS, effdur] = AnWin(S.SPT, AW, irep, isub);
elseif isnumeric(S),
   if ~isequal(':', irep) | ~isequal(':', isub), 
      error('Selecting repetitions or subsequences is not allowed for numerical format of spike times.');  
   end
   if length(AW)<2, AW = [0 AW]; end;
   WS = S(find((S>=AW(1))&(S<AW(2))));
   effdur = diff(AW);
else, error(['Invalid data class "' class(S) '" for AnWin.']);
end
%======================================
function T = localPool(S);
% pool S across across reps
[Ncond, Nrep] = size(S);
for icond=1:Ncond,
    T{icond,1} = [S{icond,:}];
end














