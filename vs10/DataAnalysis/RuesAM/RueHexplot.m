function S = RueHexplot(FN1, FN2, doPrint);
% function RueHexplot(FN1, FN2, doPrint);
%   Combine FRAs and cross corr, ipsi & contra in one figure
% Example
%   RueHexplot('p90206A', 'p90206B')

if nargin<3, doPrint=0; end

FN1c = [FN1 'c']; FN1i = [FN1 'i'];
FN2c = [FN2 'c']; FN2i = [FN2 'i'];
hf = figure; 
set(hf,'units', 'normalized', 'position', [0.28 0.138 0.523 0.762])
set(hf,'paperpos',[0.25    0.25    8.0000    10.0000])

subplot(3,2,1);
FACA1c = RueceptiveField(FN1c,1);
subplot(3,2,2);
FACA1i = RueceptiveField(FN1i,1);
subplot(3,2,3);
FACA1c = RueceptiveField(FN2c,1);
subplot(3,2,4);
FACA1i = RueceptiveField(FN2i,1);

subplot(3,2,5);
FCCA = RueCrossShuf(FN1c, FN2c,1);
subplot(3,2,6);
FCCA = RueCrossShuf(FN1i, FN2i,1);

if doPrint,
    print(hf, '-djpeg', fullfile('D:\Data\RueData\IBW2\FCCA', [FN1 '_' FN2]));
end


% 
% subplot(2,2,3);
% %RueCross(FN1, FN2);
% FACA1(FACA1<0)=0;
% FACA2(FACA2<0)=0;
% Pred = sqrt(FACA1.*FACA2);
% %Pred(Pred<0.1) = inf;
% Rueceplot(FCCA-Pred,-1, 1,1,1); % 1,1,1 = draw xlabels, ylabels, colorbars



