function V = JL_APparam(Uidx, ParamName, ParamVal);
% JL_APparam - get/set AP parameter in database
%    Pval = JL_APparam(Uidx, Pname) returns the AP parameter Pname of the
%    recording having unique recording index Uidx. Pname is one of
%        APthrSlope, CutoutThrSlope, APwindow_start, APwindow_end.
%    The value is retrieved from a database. 
%
%    JL_APparam(Uidx, '-all') retrieves all the parameter values in a
%    struct.
%
%    JL_APthrSlope(Uidx, Pname, Pval) overrides the value in the database.
%
%    See also JLgetrec.

persistent PPP

CFN = fullfile(processed_datadir, 'JL', 'JL_APparam', mfilename);
if isempty(PPP), % initialize
    PPP = getcache(CFN,1);
end


V = PPP(ismember([PPP.UniqueRecordingIndex], Uidx));
if isempty(V),
    error('Requested UniqueRecordingIndex not found in JLdbase.');
end
if nargin==2, % get specific param
    if ~isfield(V, ParamName),
        error(['Unknown AP param '''  ParamName '''.']);
    else,
        V = [V.(ParamName)];
    end
elseif nargin==3, % set
    if ~isscalar(Uidx),
        error('Uidx must be single number when updating a AP parameter.');        
    end
    if ~isfield(V, ParamName),
        error(['Unknown AP param '''  ParamName '''.']);
    end
    PPP([PPP.UniqueRecordingIndex]==Uidx).(ParamName) = ParamVal;
    putcache(CFN, 1, 1 , PPP); 
    V = ParamVal;
end










