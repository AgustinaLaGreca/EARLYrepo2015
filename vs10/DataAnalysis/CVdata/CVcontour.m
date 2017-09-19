function CVcontour(DS, Chan, icond, Cdelay, Ref)
% CVcontour - contour plot of gain and phase
%    CVcontour(DS, Chan, icond, Cdelay, Ref)
%    Contourf graphs of gain and phase data. 
%    Input args
%        DS: array of datasets.
%      Chan: AD channel
%     icond: selection of stim conditions. default is 0 == all conditions
%    Cdelay: compensating delay (ms) for phase plots.
%       Ref: transfer object used as reference (e.g. ME trf obtained w ds2trf)
%
%   Examples

[icond, Cdelay, Ref] = arginDefaults('icond/Cdelay/Ref',0,0,[]); 

AA = apple(DS, Chan, icond, Ref); 
AA = sortAccord(AA, -[AA.baseSPL]); % high to low SPL
Freq = AA(1).Fprim/1000; % stim freq in kHz
if isequal('-', DS(1).ProfileType),
    SPL = [AA.baseSPL]; % stim intensity per component in dB SPL
    YL = 'Intensity/cmp (dB SPL)';
    Apdx = 'flat';
elseif isequal('single cmp', DS(1).ProfileType),
    SPL = [AA.baseSPL] + [AA.SPLjump]; % stim intensity of varied component in dB SPL
    YL = 'Suppressor intensity (dB SPL)';
    Fprofile = unique([AA.Fprofile])/1e3;  %kHz profile freq
    Apdx = sprintf(' base SPL = %d dB, %0.1f-Hz comp amplified', AA(1).baseSPL, Fprofile);
else,
    error('NYI');
end
Alpha = cat(1,AA.Alpha);
GG = cat(1, AA.Gain)+pmask(Alpha<=0.05);
PH = cat(1, AA.Phase)+pmask(Alpha<=0.05);
for ispl=1:numel(SPL),
    PH(ispl, :) = delayPhase(PH(ispl,:), Freq, Cdelay);
    PH(ispl, :) = mod(PH(ispl, :), 1);
end

set(figure,'units', 'normalized', 'position', [0.302 0.256 0.641 0.66]);
figh = gcf;
colormap(circcolormap);
% ===Gain=
subplot(2,1,1);
Gmin = 5*floor(min(GG(:))/5)
Gmax = 5*ceil(max(GG(:))/5)
Bdr = Gmin:5:Gmax;
contourf(Freq, SPL, GG, Bdr); 
caxis(Gmin + [0 0.9*(Gmax-Gmin)]); 
xlog125; 
hca = colorbar;
ylabel(YL);
title([AA(1).ExpName '  rec ' sprintf('%d ', [AA.irec]) '  ' Apdx]);
% ===Phase=
subplot(2,1,2);
%contourf(Freq, SPL, PH, 0:0.05:1); 
[fFreq, fSPL, fPH] = local_finerphase(Freq, SPL, PH);
hp = pcolor(fFreq, fSPL, fPH); 
set(hp, 'edgecolor', 'none')
caxis([0 1]); 
xlog125; 
hca = colorbar;
ylabel(YL);
title(['cdelay=' num2str(Cdelay) ' ms']);

% ===========================================
function [fFreq, fSPL, fPH] = local_finerphase(Freq, SPL, PH);
PH = exp(2*pi*i*PH); % convert to complex unit vector to enable lin interp
fFreq = logispace(min(Freq), max(Freq), 200);
fSPL = linspace(min(SPL), max(SPL), 100);
fPH = interp2(Freq(:)', SPL(:), PH, fFreq(:).', fSPL(:));
fPH = mod(cangle(fPH),1);


