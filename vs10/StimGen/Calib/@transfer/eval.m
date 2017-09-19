function M = eval(T, Freq, BeFreqTolerant);
% transfer/eval - evaluate transfer fcn at given frequencies
%    eval(T, Freq) returns an array containing the complex transfer
%    function of Transfer object T evaluated at frequencies Freq [Hz].
%    Evaluation outside T's frequency range produces an error.
%    Note that the evaluation does NOT include any phase compensation to
%    incorporate the T's wideband delay (see below).
%
%    eval(T, Freq, 1) suppresses out-of-frequency-range errors and returns
%    an arbitrary value of one outside the range.
%
%    eval(T, Freq, TolerateOutOfRange) suppresses out-of-frequency-range errors and returns
%    an arbitrary value of one outside the range.
%
%    See Transfer, Transfer/ztrf.

BeFreqTolerant = arginDefaults('BeFreqTolerant', 0);

if ~isfilled(T),
    error('Transfer object T is not filled.');
end
% linear interpolation
Sz = size(Freq);
if BeFreqTolerant,
    M = interp1(T.Freq, T.Ztrf, Freq(:), 'linear', 1); % 1 is the "extrap" value, i.e. default outside T.Freq range
else,
    M = interp1(T.Freq, T.Ztrf, Freq(:));
    if any(isnan(M)),
        error('Stimulus frequencies exceed calibration range.');
    end
end
M = reshape(M(:), Sz); % same shape as Freq input arg




