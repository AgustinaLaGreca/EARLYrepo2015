function IR= infreq_range(EC, Freq, ichan);
% earcalib/infreq_range - true if frequency is in range of Earcalib object
%    infreq_range(EC, Freq) returns Freq if Freq is within the frequency
%    range of Earcalib object EC. Freq may be an array.
%
%    infreq_range(EC, Freq,ichan) uses DA channel ichan. Default ichan=1.
%
%    See Earcalib/frequency.

if nargin<3, ichan=1; end

if ~isfilled(EC),
    error('Transfer object T is not filled.');
end

Trf = EC.Transfer(ichan);
IR = (Freq>=fmin(Trf)) & (Freq<=fmax(Trf));

