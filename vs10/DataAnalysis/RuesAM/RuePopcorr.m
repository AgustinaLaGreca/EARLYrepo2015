function S = RuePopcorr(FN1, FN2, excludeZerodB);
% RuePopcorr - xcor analysis for population test. Whole audiogram is
% treated as one stimulus, i.e., the recordings are not split into the
% responses to different tones.
%   S = RuePopcorr(FN1, FN2, excludeZerodB) returns a single S with fields
%        FN1: cell #1
%        FN2: cell #2
%  excludeZerodB: true -> zero dB conditions are excluded
%     isPair: true for a cell pair
%     isSame: true when FN1 & FN2 are the same
%
%      Nrep1: # repetitions for cell 1
%      Nrep2: # repetitions for cell 2
%
%      mean1: mean values of waveforms of reps # 1..Nrep1 of cell 1
%      mean2: mean values of waveforms of reps # 1..Nrep2 of cell 2
%      var1: variances of waveforms of reps # 1..Nrep1 of cell 1
%      var2: variances of waveforms of reps # 1..Nrep2 of cell 2
%
%       rho: normalized Xcorr values in Nrep1xNrep2 matrix
%      urho: unnormalized  Xcorr values in Nrep1xNrep2 matrix
%    meanrho: mean across all rep-pairs of rho
%    meanrho_nd: mean across all nondiagonal rep-pairs of rho
%
%     MWvar1: variance of waveform average across reps of cell 1
%     MWvar2: variance of waveform average across reps of cell 2
%      MWrho: normalized Xcorr between mean waveforms
%     MWurho: unnormalized Xcorr between mean waveforms

if nargin>1 && isnumeric(FN2), 
    excludeZerodB = FN2; 
elseif nargin<3,
    excludeZerodB = 0; 
end

if (nargin<2) || (nargin==2 && isnumeric(FN2)),
    FN2 = FN1;
elseif isequal('B', FN2),
    FN2 = strrep(FN1, 'A', 'B');
end
% sort F1 and F2 
FN = sort({FN1 FN2});
FN1 = FN{1}; FN2 = FN{2};

% try to retrieve S from cache
CFN = mfilename;
CPar = {FN1, FN2, excludeZerodB, 2};
S = getcache(CFN, CPar);
if ~isempty(S), return; end

disp(['=== ' FN1, '  ' FN2 '===']);
% compute all fields of S 
isPair = ~isequal(FN1, FN2) && (...
    isequal(strrep(FN1,'A', 'B'), FN2) ...
    || isequal(strrep(FN1,'B', 'A'), FN2)   );
isSame = isequal(FN1, FN2);

Nrep1 = RueNrep(FN1);
Nrep2 = RueNrep(FN2);

S = CollectInStruct(FN1, FN2, excludeZerodB, isPair, isSame, '-', ...
    Nrep1, Nrep2, '-');

S.rho = [];
for irep1=1:Nrep1,
    [Y1, dum, dt] = readRueRep(FN1, irep1);
    if excludeZerodB, Y1 = Y1(:, 1:464); end
    for irep2=1:Nrep2,
        Y2 = readRueRep(FN2, irep2);
        if excludeZerodB, Y2 = Y2(:, 1:464); end
        qq = cov(Y1(:),Y2(:));
        S.urho(irep1,irep2) = qq(1,2);
        S.var1(irep1,irep2) = qq(1,1);
        S.var2(irep1,irep2) = qq(2,2);
    end
end
S.rho = S.urho./sqrt(S.var1.*S.var2);


putcache(CFN, 1e4, CPar, S);


% ================================







