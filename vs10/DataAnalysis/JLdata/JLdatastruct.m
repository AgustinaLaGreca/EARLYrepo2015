function S = JLdatastruct(X, RecID, icond)
% JLdatastruct - return struct that fully describes a single JLbeat recording
%    JLdatastruct(X) returns a struct describing the details of a single
%    JL recording, i.e., a single frequency condition of a beat sweep.
%    This standard recording descriptor which serves as a start for 
%    specific analysis programs.
%    X can have one of the following formats:
%        - a struct as returned by JLdbase
%        - a struct containing fields: ExpID, RecID, icond.
%        - a single number equal to the UniqueRecordingIndex of the recording
%    When X is an array, S is also an array.
%
%    Alternative input formats are
%        S = JLdatastruct(ExpID, RecID, icond)
%        S = JLdatastruct({ExpID, RecID, icond})
%    Here icond (but not the another input args) may be an array.


%clear JLdbase % force JLdbase to read from disk @ next call

if nargin>1, % convert to struct format
    ExpID = X;
    clear X
    IC = icond;
    for ii=1:numel(IC),
        icond = IC(ii);
        X(ii) = CollectInStruct(ExpID, RecID, icond);
    end
end

if iscell(X), % unpack cells
    S = JLdatastruct(X{:});
    return;
end

if numel(X)>1, % recursive handling of multiple recordings
    for ii=1:numel(X),
        S(ii) = JLdatastruct(X(ii));
    end
    return;
end


%======single recording from here =============
D = JLdbase;
if isstruct(X) && isfield(X,'UniqueRecordingIndex'), % e.g., JLdbase output
    d = D([D.UniqueRecordingIndex]==X.UniqueRecordingIndex);
elseif isstruct(X) && isfield(X,'ExpID'), % use ExpID etc
    [iexp, icell, iseq] = local_parse(X);
    d = D((iexp==[D.iexp]) & (icell==[D.icell]) & (iseq==[D.iseq]) & (X.icond==[D.icond]));
elseif isnumeric(X),
    d = D([D.UniqueRecordingIndex]==X);
else,
    error('Invalid format of input arg X.');
end

if isempty(d),
    display(X);
    error('Specified Recording not found');
elseif numel(d)>1,
    for ii=1:numel(d),
        disp(d(ii));
    end
    error('Multiple recordings match specification.');
end

DBdir = fullfile(processed_datadir,  '\JL\JLdbase');
FN = fullfile(DBdir, ['Spart_' num2str(d.i_SPfile)]);
load(FN, 'Spart'); % this loads variable Spart (see JLbuild_database)
S = Spart([Spart.UniqueRecordingIndex]==d.UniqueRecordingIndex);

%=================================================
function [iexp, icell, iseq] = local_parse(X);
R = X.RecID;
idash = [find(R=='-') numel(R)+1];
icell = str2num(R(1:idash(1)-1));
iseq = str2num(R(idash(1)+1:idash(2)-1));
E = X.ExpID;
if isletter(E(end)),
    E = E(1:end-1); % 216B -> 216
end
iexp = str2num(E(end-2:end)); % three last digits




