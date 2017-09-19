function Mess = ZWOAEtest(Fz, Fs);
% ZWOAEtest - test zwuisness of candidate stimulus frequencies
%   ZWOAEtest(Fz, Fs) tests whether the array Fz (zwuis freqs) and 
%   the single number Fs (single freq) yield non-overlapping ZWOAE
%   components. If the frequencies pass the test, '' is returned, otherwise
%   ZWOAEtest returns a message indicating the problem.
%
%  See also ZWOAEmatrices.

% put all primary components in single col vector Fprim (see ZWOAEmatrices)
Fz = sort(Fz(:));
Fprim = [Fz; Fs];

% compute the complete aggregate of ZWOAEs
Nzw = numel(Fz); % # zwuis comp
MS = ZWOAEmatrices(Nzw);
Fdist = MS.Mall*Fprim; % all ZWOAE distortion products
% numerical errors might spoil equality, so round the values
if localAreUnique(Fdist),
    Mess = ''; % okay
    return;
end
% not good - find details
Enearfar = ~localAreUnique([MS.Mnear; MS.Mfar]*Fprim);
Enearsup = ~localAreUnique([MS.Mnear; MS.Msupall]*Fprim);
Efarsup = ~localAreUnique([MS.Mfar; MS.Msupall]*Fprim);
Mess = ['Primary frequencies yield non-unique ZWOAE components. Clash between '];
if Enearfar, Mess = [Mess, 'near/far ']; end
if Enearsup, Mess = [Mess, 'near/supp ']; end
if Efarsup, Mess = [Mess, 'far/supp ']; end
if ~(Enearfar | Enearsup | Efarsup), Mess = [Mess '??/??']; end
Mess = [Mess ' type components.']; 

% =========================
function au=localAreUnique(Fdist);
% test uniqueness of components, correcting for numerical "noise"
[fmin, fmax] = minmax(Fdist);
RR = fmax-fmin;
df = max(1e-10,RR/1e6); % smallest freq difference considered nonzero
Fdist = df*round(Fdist/df); % round to multiple of df
Ndist = numel(Fdist); % # ZWOAE components
Ndiff = numel(unique(Fdist)); % # different ZWOAE components
au = isequal(Ndist, Ndiff);

