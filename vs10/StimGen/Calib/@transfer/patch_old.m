function P = patch(S,T);
% transfer/patch - combine transfer functions increase spectral range
%    patch(S,T) is the Transfer object characterizing S and T "in parallel."
%    In order to combine S and T in this order, the following
%    properties of S and T must be identical:
%         Q_resp and Q_stim
%         Ref_resp and Ref_stim
%         dBRef_resp dBRef_stim
%    In the region of overlap of the frequency bands of S and T, the mean
%    of the two complex transfer functions is taken.
%
%    See Transfer, Transfer/mtimes.

if ~isequal(S.Q_resp, T.Q_resp),
    error('Q_resp properties of S and T must be identical.');
end
if ~isequal(S.Q_stim, T.Q_stim),
    error('Q_stim properties of S and T must be identical.');
end

if ~isequal(S.Ref_stim, T.Ref_stim),
    error('Ref_stim properties of S and T must be identical.');
end
if ~isequal(S.Ref_resp, T.Ref_resp),
    error('Ref_resp properties of S and T must be identical.');
end

if ~isequal(S.dBref_stim, T.dBref_stim),
    error('dBref_stim properties of S and T must be identical.');
end
if ~isequal(S.dBref_resp, T.dBref_resp),
    error('dBref_resp properties of S and T must be identical.');
end

% make sure S is the one starting at lower freqs. 
if fmin(S)>fmin(T),
    [S,T] = swap(S,T);
end
if fmax(S)<fmin(T),
    error('S and T may not have a spectral gap between them.');
end
% disp(['S: ' num2str(fmin(S)) '-' num2str(fmax(S)) ' Hz']);
% disp(['T: ' num2str(fmin(T)) '-' num2str(fmax(T)) ' Hz']);

P = S;
Ns = numel(S.Ztrf); Nt = numel(T.Ztrf);
N = Ns+Nt;
[Fmin, Fmax] = minmax([S.Freq; T.Freq]);
P.Freq = logispace(Fmin, Fmax, N).'; % new freq axis
P.Ztrf = nan(size(P.Freq));
P.WB_delay_ms = max([S.WB_delay_ms T.WB_delay_ms]);
S = setWBdelay(S, P.WB_delay_ms);
T = setWBdelay(T, P.WB_delay_ms);

% assess the contributions of S and T to P over new freq range by simple
% interpolation. NaNs will them indicate the validity of the span.
Sz = interp1(S.Freq, S.Ztrf, P.Freq);
Tz = interp1(T.Freq, T.Ztrf, P.Freq);

inS = find(~isnan(Sz));
inT = find(~isnan(Tz));
inBoth = sort(intersect(inS,inT));
inS = setdiff(inS,inBoth);
inT = setdiff(inT,inBoth);
P.Ztrf(inS) = Sz(inS);
P.Ztrf(inT) = Tz(inT);
NsamBoth = numel(inBoth);
if NsamBoth>0,
    w = linspace(0,1,NsamBoth+2); w = w(2:end-1).';
    P.Ztrf(inBoth) = w.*Sz(inBoth) + (1-w).*Tz(inBoth); % fade S out, fade T in
end

if any(isnan(ztrf(P))),
    error('NaNs in Patch???!');
end













