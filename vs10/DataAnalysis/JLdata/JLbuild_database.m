function D = JLbuild_database(Sp);
% JLbuild_database.m - build database for JLbeat data
if nargin<1,
    load(fullfile(processed_datadir, '\JL\JLvarStats\varstats_with_partners'), 'Sp');
end

DBdir = fullfile(processed_datadir,  '\JL\JLdbase');

% generate lightweight, ID-only, struct array D
D = [];
[prev_sp, prev_d] = deal([]);
for ii=1:numel(Sp),
    ii
    sp = Sp(ii);
    d.iexp = local_iexp(sp.ExpID);
    d.icell = sp.icell;  % cell count within exp
    newcell = 1; % default assumption
    if isempty(prev_sp), % first cell
        d.icell_run = 1;
    elseif ~isequal(prev_sp.UniqueCellIndex, sp.UniqueCellIndex), % next cell
        d.icell_run = prev_d.icell_run+1;
    else, % same cell
        d.icell_run = prev_d.icell_run;
        newcell=0;
    end
    d.iseq = local_iseq(sp);
    d.SPL = max([sp.SPL1 sp.SPL2]);
    d.chan = sp.StimType;
    if isempty(prev_sp), % first cell
        d.iseries = 1;
        d.iseries_run = 1;
    elseif newcell, 
        d.iseries = 1;
        d.iseries_run = prev_d.iseries_run + 1;
    elseif ~isequal(prev_sp.UniqueSeriesIndex , sp.UniqueSeriesIndex), % new series in same cell series
        d.iseries = prev_d.iseries + 1;
        d.iseries_run = prev_d.iseries_run + 1;
    else, % same series
        d.iseries = prev_d.iseries;
        d.iseries_run = prev_d.iseries_run;
    end
    d.icond = sp.icond;
    d.freq = sp.Freq1;
    d.UniqueRecordingIndex = sp.UniqueRecordingIndex;
    d.MinutesSinceRecStart = JLtimeSinceRecstart(sp.UniqueRecordingIndex);
    d.i_SPfile = nan;
    prev_sp = sp;
    prev_d = d;
    D = [D, d];
end

% cut Sp in blocks and save those
BlockSize = 100;
iblock = 0;
for ii=0:numel(Sp),
    if rem(ii,BlockSize)==0 || ii>=numel(Sp),
        if iblock>0, % save previous Spart
            FN = fullfile(DBdir, ['Spart_' num2str(iblock)]);
            save(FN, 'Spart', '-V6');
        end
        iblock = iblock + 1
        Spart = [];
    end
    if ii>=numel(Sp), break; end
    s = local_extract(D(ii+1),Sp(ii+1));
    D(ii+1).i_SPfile = iblock;
    Spart = [Spart, s];
end

% save D, including the pointgers to the files containg the block of Sp
SFN = fullfile(DBdir, 'JLdbase');
save(SFN, 'D');


%=========================
function iexp = local_iexp(ExpID);
if isletter(ExpID(end)),
    ExpID = ExpID(1:end-1); % 216B -> 216
end
iexp = str2num(ExpID(end-2:end)); % three last digits

function d = local_extract(d0,sp);
sp1 = structpart(sp, {'ExpID' 'RecID' 'seriesID' 'icond'});
sp2 = structpart(sp, {'UniqueCellIndex' 'UniqueSeriesIndex' 'UniqueRecordingIndex', ...
    'ABFname' 'ABFdir' 'abfdate' 'abfdatenum'});
sp3 = structpart(sp, {'APthrSlope' 'APwindow_start' 'APwindow_end' 'CutoutThrSlope'  'Ttrans' 'JLcomment' 'MHcomment'});
sp4 = structpart(sp, {'StimType' 'Freq1' 'Freq2' 'Fbeat' 'burstDur' 'SPL1' 'SPL2' 'preDur'});
sp5 = structpart(sp, {'varName' 'TTT'});
sp6 = structpart(sp, {'I_partners' 'C_partners' 'B_partners' 'hasMonPartners'});
d = structJoin(d0, local_sep('ID'), sp1, sp2, sp5, ...
    local_sep('AnaParam'), sp3,  local_sep('StimParam'), sp4,  ...
    local_sep('Partners'), sp6);

function s = local_sep(Str);
Nsc = 25-numel(Str);
s.([Str repmat('_', [1 Nsc])]) = '_____________________';

function iseq = local_iseq(sp);
idash = find(sp.RecID=='-', 1 ,'first');
iseq = str2num(sp.RecID(idash+1:end-4));






