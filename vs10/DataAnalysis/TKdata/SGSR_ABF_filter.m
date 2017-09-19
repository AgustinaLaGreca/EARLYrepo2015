function s = SGSR_ABF_filter(S, ExpID, RecID);
% SGSR_ABF_filter - select element of SGSR_ABF link data base
%   SGSR2ABF(Sall, 'RG09110', '3-1-FS') selects that element of Sall that has
%   the specified ExpID and RecID. Sall is the SGSR_ABF database produced
%   by TKpool. RecID must be unique prefix.
%
%   See also compile_all_SGSR_ABFs.

% select exp
iexp = strmatch(ExpID, {S.SGSRExpID}, 'exact');
S = S(iexp);
% select rec
irec = strmatch(upper(RecID), {S.SGSRidentifier});
if isempty(irec),
    error(['No ABF recordings matching ''' ExpID, '/' RecID '''.']);
elseif numel(irec)>1,
    error(['Multiple ABF recordings matching ''' ExpID, '/' RecID '''.']);
end
s = S(irec);







