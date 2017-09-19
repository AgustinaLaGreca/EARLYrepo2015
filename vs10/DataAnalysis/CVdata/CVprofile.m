function CVP = CVprofile(DS, Chan, iCond, Ref, GainBounds)
% CVprofile - analysis of zwuis data with varied spectral profile 
%    CVprofile(DS, Chan, iCond, Ref, GainBounds)
%    Inputs are 
%           DS: array of ZW datasets
%         Chan: DA channel specification (see dataset/Anachan)
%        iCond: selection of stimulus condition. Default is 0, meaning all
%               conditions. Specify a cell array to select a different stim
%               subset for each element of the aray DS.
%          Ref: transfer object acting as "calibration". Typically obtained
%               from middle-ear measurements and (See dataset/ds2trf).
%   GainBounds: array holding boundaries for gain contour plot. Default
%               5-dB steps from min to max Gain value.

Nfreq = 100; % # freq components for uniform freq spacing

[iCond, Ref, GainBounds] = arginDefaults('iCond/Ref/GainBounds', 0, [], []);

AA = apple(DS, Chan, iCond, Ref);
Na = numel(AA);
% boundary freqs
allFprim = cat(1, AA.Fprim);
Fmin = max(allFprim(:,1));
Fmax = min(allFprim(:,end));
Freq = linspace(Fmin, Fmax, Nfreq);
for ia=1:Na,
    Gain(ia,:) = interp1(AA(ia).Fprim, AA(ia).Gain, Freq);
    phasor = exp(2*pi*i*AA(ia).Phase);
    phasor = interp1(AA(ia).Fprim, phasor, Freq);
    Phase(ia,:) = cunwrap(cangle(phasor));
end
if isempty(GainBounds),
    minGain = 5*floor(min(Gain(:))/5);
    maxGain = 5*ceil(max(Gain(:))/5);
    GainBounds = minGain:5:maxGain;
end
Fprofile = cat(1, AA.Fprofile);
baseSPL = cat(1, AA.baseSPL);
SPLjump = cat(1, AA.SPLjump);


CVP = CollectInStruct(Freq, Gain, Phase, Fprofile, baseSPL, SPLjump);

if all(Fprofile==0), % prob SPL was varied
    contourf(CVP.Freq/1e3, CVP.baseSPL, CVP.Gain, GainBounds);
    ylabel('Base SPL (dB)');
elseif ~isscalar(unique(SPLjump)), % SPLjump  was varied
    contourf(CVP.Freq/1e3, CVP.SPLjump, CVP.Gain, GainBounds);
    ylabel('Profile frequency (kHz)');
else, % Fprofile was varied
    contourf(CVP.Freq/1e3, CVP.Fprofile/1e3, CVP.Gain, GainBounds);
    ylabel('Profile frequency (kHz)');
end
xlabel('Stimulus frequency (kHz)'); 
title([AA(1).ExpName '  rec ' sprintf('%d ', unique([AA.irec])) ' ' unique([AA.ProfileType])]);





