function T= transfer(EC, ichan);
% Earcalib/Transfer - transfer object of Earcalib object
%    Transfer(EC) returns the D/A Volt -> dB SPL transfer function of Earcalib
%    object EC. This is a 1-by-2 Transfer object corresponding to channels 
%    1 and 2 of the D/A converter.
%
%    Transfer(EC,ichan) only returns the component corresponding to DA
%    channel ichan.
%
%    See Transfer, Transfer/measure, Transfer/resp_amp.

if nargin<2, ichan=1:numel(EC.Transfer); end
T = EC.Transfer(ichan);



