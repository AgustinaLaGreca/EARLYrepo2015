function S = RueTriplot(FN1, FN2);
% function RueTriplot(FN1, FN2);
%   Combine FRAs and cross corr in one figure
% Example
%   S = RueTriplot('p90206Ac', 'p90206Bc')
%   S contains all kind of incomprehensible similarity metrics ;)

hf = figure; 
set(hf,'units', 'normalized', 'position', [0.28 0.259 0.634 0.641])
subplot(2,2,1);
FACA1 = RueceptiveField(FN1,1);
subplot(2,2,2);
FACA2 = RueceptiveField(FN2, 1);
subplot(2,2,3);

subplot(2,2,4);
FCCA = RueCrossShuf(FN1, FN2,1);

subplot(2,2,3);
%RueCross(FN1, FN2);
FACA1(FACA1<0)=0;
FACA2(FACA2<0)=0;
Pred = sqrt(FACA1.*FACA2);
%Pred(Pred<0.1) = inf;
Rueceplot(FCCA-Pred,-1, 1,1,1); % 1,1,1 = draw xlabels, ylabels, colorbars


% similarity metrics ..
SimpleSim = mean(FCCA(:));
TotalWeight1 = sqrt(sum(FACA1(:).^2));
TotalWeight2 = sqrt(sum(FACA2(:).^2));
Overlap = sqrt(sum( abs(FACA1(:).*FACA2(:)) ));
NormOverlap = Overlap/sqrt(TotalWeight1*TotalWeight2);

TotalWeight = sqrt(sum(FACA1(:).^2 + FACA2(:).^2)); % total weights for FACA1 and FACA2 combined
Sweight = sqrt(max(FACA1.^2, FACA2.^2)); % high weight where *either* cell respond
Sweight = Sweight/sum(Sweight(:));

Pweight = abs(sqrt(FACA1.*FACA2)); % high weight only where *both* cells respond
Pweight = Pweight/sum(Pweight(:));

Ssimilarity = sum(Sweight(:).*FCCA(:));
Psimilarity = sum(Pweight(:).*FCCA(:));
Qsimilarity = sum(FCCA(:))/TotalWeight;

isPair = isRuePair(FN1, FN2);
S = CollectInStruct(FN1, FN2, SimpleSim, NormOverlap, TotalWeight, TotalWeight1, TotalWeight2, Overlap, NormOverlap, '-', ...
    Ssimilarity, Psimilarity, Qsimilarity, isPair);

