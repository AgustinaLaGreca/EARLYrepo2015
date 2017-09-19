function W = Waveform(P,j)
% playlist/waveform - waveforms stored in playlist object.
%
%   Waveform(P) or P.Waveform returns all waveforms stored in P in a cell
%   array of numerical columns vectors.
%
%   Waveform(P,j) or P.Waveform{j} returns only the jth waveform.
%
%   See playlist/double, "methods playlist".

W = P.Waveform;
if nargin>1, % select single jth waveform
    if ~isscalar(j),
        error('Index j must be integer scalar.');
    end
    error(numericTest(j, 'posint/noninf', 'Index j is '))
    W = W{j};
end






