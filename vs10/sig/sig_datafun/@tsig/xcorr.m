function Y = xcorr(S,T,tmax,SCL);
% tsig/xcorr - XCORR for tsig objects: crosscorrelation function
%    X=xcorr(S,T) a tsig object whose kth channel contains the
%    crosscorrelation function between the kth channels of S and T, i.e.,
%    S(k) and T(k). The onsets of S and T are accounted for as follows.
%    First the time domains of each channel S(k) and T(k) are equalized by
%    padding zeros. Second these zer-padded waveforms are crosscorrelated.
%
%    xcorr(S,T,[tmin tmax]) restricts the lag from tmin ms to tmax ms.
%    xcorr(S,T,tmax]) is the same as xcorr(S,T,[-tmax tmax]).
% 
%    xcorr(S,T,tmax,scaleopt) normalizes the correlation according to
%    scaleopt, which is one of 'biased' 'unbiased' 'coeff' 'none'.
%    See Matlab's XCORR for details. Note: unlike XCORR, the default
%    normalization for tsig/xcorr is 'coeff'.
% 
%    See also XCORR, tsig/align.

if nargin<3,
    tmax = [];
end
if nargin<4,
    SCL = 'coeff'; % default normalization. See help text
end

if ~isempty(tmax),
    error(numericTest(tmax,'rreal', 'Time lag argument is '))
    if isscalar(tmax), tmax=[-1 1]*abs(tmax); end % see help text
    if ~isequal([1,2], size(tmax)),
        error('Time lag argument must be single number or 1 x 2 row vector [tmin tmax].');
    end
    tmax = sort(tmax);
end

% channelwise time alignment
[S,T, Mess] = align(S,T);
switch Mess,
    case 'void', error('Void tsig objects may not be crosscorrelated.');
    case 'nchan', error('Incompatible channel counts between tsig objects.')
    case 'fsam', error('Tsig object having different sample frequencies cannot be be crosscorrelated.');
end

[SCL, Mess] = keywordMatch(SCL,{'biased' 'unbiased' 'coeff' 'none'},'scaling option');
error(Mess);

DT = DT(T);
Nmaxlag = round(max(abs(tmax))/DT);
if isempty(Nmaxlag),
    Nmaxlag = 2*Nsam(S)-1; % full range; see XCORR
end
Nmaxlag = SameSize(Nmaxlag,S.Waveform); % cellfun needs uniform input sizes
my_xcorr = @(x,y,ml)xcorr(x,y,ml,SCL); % portable xcorr with additional args
Y.Waveform = cellfun(my_xcorr,S.Waveform ,T.Waveform , num2cell(Nmaxlag), 'UniformOutput', false);
Yonset = -Nmaxlag*DT;
Y = tsig(S.Fsam,Y.Waveform,Yonset);
% provide restriction of lags if requested
if ~isempty(tmax),
    Y = cut(Y,tmax(1),tmax(2), tmax(1));
end




