function DispRueCorr7(Id1,Id2, PostFix);
% DispRueCorr7 - display previosuly computed across-unit analysis
%   DispRueCorr PostFix
%  
if nargin<3, 
    PostFix = '';
end

DPostfix = [strrep(PostFix,'_', ' conditions, ') 'ms'];
Ddir = 'D:\Data\RueData';
load(fullfile(Ddir,['RueCorNorm' PostFix])); % this loads RueCorr struct
RueCorNorm.AZ = RueCorNorm.ZA; %#ok<NODEF>

% field names, e.g. AB, AA, BB
FX = sort(upper([Id1 Id2]));
F1 = FX([1 1]); F2 = FX([2 2]);
C1 = median(RueCorNorm.(F1));
C2 = median(RueCorNorm.(F2));
CX = median(RueCorNorm.(FX));
allCXn = RueCorNorm.(FX)/sqrt(C1*C2);
CXn = median(allCXn);

CorBins = linspace(-1,1.5,60); % bin centers for histograms
figure; set(gcf,'units', 'normalized', 'position', [0.399 0.257 0.383 0.66])
subplot(3,1,1);
hist(RueCorNorm.(F1),CorBins); xlim([-0.5 1]); title([F1 ': ' num2str(C1,2) ' (' DPostfix ')']);
subplot(3,1,2);
hist(RueCorNorm.(F2),CorBins); xlim([-0.5 1]); title([F2 ': ' num2str(C2,2)]);
subplot(3,1,3);
hist(RueCorNorm.(FX),CorBins); 
hh = findobj(gcf,'type','patch'); set(hh,'facecolor',[1 0.7 0.7]);
hold on
hist(allCXn,CorBins); xlim([-0.5 1]); 
hh = findobj(gca,'type','patch'); set(hh(1),'facecolor',[0.7 1 1]);
title([FX ': ' num2str(CX,2) '; normalized: ' num2str(CXn,2)]);
xlabel('Normalized crosscorrelation');














