function P = mtimes(S,T);
% transfer/times - cascade of two transfer functions
%    S*T is the Transfer object characterizing S and T "in series": first
%    S then T. In order to combine S and T in this order, the following
%    pairs of properties must be identical:
%         S.Q_resp and T.Q_stim
%         S.Ref_resp and T.Ref_stim
%         S.dBRef_resp and T.dBRef_stim
%    The Transfer function is only defined the region of overlap of the
%    frequency bands of S and T.
%    The Description of S*T is simply the concatenation of S.Descr &
%    T.Descr.
%
%    c*S, where c is a numeric scalar and T a transfer object, multiplies 
%    the complex transfer function of T by c. S*c is the same thing.
%
%    See Transfer, Transfer/inverse.

if isnumeric(S), % simple scalar multiplication
    P = T;
    P.Ztrf = S*P.Ztrf;
    return
elseif isnumeric(T), % simple scalar multiplication
    P = S;
    P.Ztrf = T*P.Ztrf;
    return
end

if ~isequal(S.Q_resp, T.Q_stim),
    error('In a product S*T, S.Q_resp and T.Q_stim must be identical.');
end
if ~isequal(S.Ref_resp, T.Ref_stim),
    error('In a product S*T, S.Ref_resp and T.Ref_stim must be identical.');
end
if ~isequal(S.dBref_resp, T.dBref_stim),
    error('In a product S*T, S.dBref_resp and T.dBref_stim must be identical.');
end

P = S;
P.Q_resp = T.Q_resp;
P.Ref_resp = T.Ref_resp;
P.dBref_resp = T.dBref_resp;
P.Description = [S.Description ' * ' T.Description];
P.CalibParam = [S.CalibParam T.CalibParam];
Fmin = max(min(S.Freq), min(T.Freq));
Fmax = min(max(S.Freq), max(T.Freq));
NfS = sum((S.Freq>=Fmin) & (S.Freq<=Fmax));
NfT = sum((T.Freq>=Fmin) & (T.Freq<=Fmax));
P.Freq = logispace(Fmin, Fmax, min(NfS,NfT));
Zs = interp1(S.Freq, S.Ztrf, P.Freq);
Zt = interp1(T.Freq, T.Ztrf, P.Freq);
P.Ztrf = Zs.*Zt;
P.WB_delay_ms = S.WB_delay_ms + T.WB_delay_ms;

















