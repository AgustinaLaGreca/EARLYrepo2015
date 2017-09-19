function DD = JLdbase(Uidx, X);
% JLdbase - database for JLbeat data
%    JLdbase() returns whole database (struct array).
%
%    JLdbase(Uidx) returns only the entry/entries having the specied
%    recording index as given by the UniqueRecordingIndex field. Uidx may
%    be an array.
%
%    JLdbase('Foo', X) only returns those entries whose Foo field equals X.
%    X must be a single value.
%
%    Example
%      JLdbase(178035411) % unique recording index equals 178035411
%      JLdbase('icell_run', 13) % all recs from this cell

persistent D

if isempty(D), % load dabse previously saved by JLbuild_database
    DBdir = fullfile(processed_datadir,  '\JL\JLdbase');
    SFN = fullfile(DBdir, 'JLdbase');
    load(SFN, 'D');
    D = sortAccord(D, [D.UniqueRecordingIndex]);
end

if nargin<1,
    DD = D;
elseif isnumeric(Uidx),
    for ii=1:numel(Uidx),
        DD(ii) = D([D.UniqueRecordingIndex]==Uidx(ii));
    end
elseif ischar(Uidx),
    DD = D([D.(Uidx)]==X);
else, error('Invalid input args.')
end





