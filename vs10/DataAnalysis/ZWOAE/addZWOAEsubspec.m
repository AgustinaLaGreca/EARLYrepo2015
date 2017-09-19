function S = addZWOAEsubspec(iGerbil, idataset, subType, Ref);
% addZWOAEsubspec - add multiple partial spectra of ZWOAE data
%    S = addZWOAEsubspec(iGerbil, iRec, SubType, Ref) combines the data
%    from recordings iRec(1) to iRec(end) into a single spectrum, by
%    applying a weighted average of the individual spectra. The weights are
%    proportional to the durations of the corresponding recordings.
%
%    See ZWOAEsubspec for more details on the inputs and outputs.
%
%  See also ZWOAEsubspec.

if nargin<4, Ref = 'stim'; end

Freq = [];
Cspec = 0; % initialize complex spectrum
TotWeight = 0;
for irec=idataset(:).',
    D = getZWOAEdata(iGerbil, irec);
    Weight = D.Stimparam.recdur; 
    S = ZWOAEsubspec(iGerbil, irec, subType, Ref);
    if ~isempty(Freq) && ~isequal(Freq, S.Fr),
        error('Incompatible frequency spectra.');
    end
    Freq = S.Fr;
    Cspec = Cspec + Weight*dB2A(S.MG).*exp(2*pi*i*S.PH);
    TotWeight = TotWeight + Weight;
end
Cspec = Cspec/TotWeight;
% convert complex spec back to magn & phase, and place combined spectrum in
% last read S
S.MG = A2dB(abs(Cspec));
S.PH = cangle(Cspec);
S.iRec = idataset; % for the record ;)







