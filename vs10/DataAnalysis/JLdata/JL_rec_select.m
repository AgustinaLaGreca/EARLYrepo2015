function S = JL_rec_select;
% JL_rec_select - recording selection
%    S = JL_rec_select returns struct array whose elements correspond to 
%    those of JLdbase. Its fields are the logicals CanBeUsed and
%    ShouldBeUsed.
%
% See also JLdbase, JLcellview.

[S, CFN, CP] = getcache(mfilename, dir('D:\USR\Gerard\MSO experiments list 20110923 GB.xls')); % cache validity depends on XLS file 
if ~isempty(S), return; end

DB = JLdbase;
Allicell = unique([DB.icell_run]);

S = [];
for icell_run = Allicell(:).',
    db = DB([DB.icell_run]==icell_run);
    db = sortAccord(db, [db.MinutesSinceRecStart]); % sort according to abs recording time
    [db.irec_per_cell] = DealElements(1:numel(db)); % add simple counting index (the one used by GB to judge admission)
    db = local_select(db);
    S = [S, structpart(db, {'UniqueRecordingIndex' 'irec_per_cell' 'CanBeUsed' 'ShouldBeUsed'})];
end

putcache(CFN, 10, CP, S);
%=====================================
function T = local_xlsfield(ic,C, Pref);
FFN = 'D:\USR\Gerard\MSO experiments list 20110923 GB.xls';
Xr = @(C)[C num2str(ic+1) ':' C num2str(ic+1)];
[dum, T] = xlsread(FFN,1, Xr(C));
if isempty(T), 
    T = '-';
else,
    T = T{1};
end
if nargin>2,
    T = [Pref ': ' T];
end

function DB = local_select(DB, F, Prf);
icell_run = DB(1).icell_run;
F = local_xlsfield(icell_run,'N');
Prf = local_xlsfield(icell_run,'P');
N = numel(DB);
if isequal('-', Prf),
    ipref = [];
else,
    Prf = strrep(Prf,' ','');
    Prf = strrep(Prf,'-',':');
    ipref = eval(Prf);
end
qpref = ismember([DB.irec_per_cell], ipref);
switch lower(F),
    case {'all' '-'}; % no restriction
        qsel = true(size(DB));
    otherwise,
        F = strrep(F,' ','');
        F = Words2cell(F,',');
        isel = [];
        for ii=1:numel(F),
            f = F{ii};
            f = strrep(f,'-',':');
            f = strrep(f,'end',num2str(N));
            isel = [isel eval(f)];
        end
        qsel = ismember([DB.irec_per_cell], isel);
end
qsel = qsel & [DB.freq]~=50 & local_cell_OK(icell_run);
qpref = qpref & qsel;
[DB.CanBeUsed] = DealElements(qsel);
[DB.ShouldBeUsed] = DealElements(qpref);

function  ok = local_cell_OK(icell_run);
%ok = ismember(icell_run, [1 4 6 7 13 14 15 16 17 18 19 20 22 23 25 30 31 32 33 35 37 38 40 41 42]);
ok = ismember(icell_run, [1 4 6 7 11 12 13 14 15 16 17 18 19 20 22 23 25 28 30 31 32 33 35 37 38 40 41 42]);
% TC = [0 0 0]; % black text
% switch icell_run
%     case {1 4 6 7 13 14 15 16 17 18 19 20 22 23 25 30 31 32 33 35 37 38 40 41 42},
%         CLR = [0 1 0];
%     case {2 3 5 9 10 26 29},
%         CLR = [0 0 1];
%         TC = [1 1 1];
%     case {8 21 24 34 36 39 43 44 45},
%         CLR = [1 0 0];
%     case {11 12 28},
%         CLR = [1 0.7 0.7];
% end


function [Hd, Str, Uidx] = local_listbox(DB, icell_run, CanBeUsed, ShouldBeUsed);
% restrict database to this cell and select relevant fields
DB = DB([DB.icell_run]==icell_run);
DB = sortAccord(DB, [DB.MinutesSinceRecStart]); % sort according to abs recording time
[DB.irec] = DealElements(1:numel(DB)); % add simple counting index
DB = DB([DB.freq] ~= 50);
[DB, cpref] = local_restrict(DB, CanBeUsed, ShouldBeUsed);
Sdisp = structpart(DB, {'irec' 'chan', 'freq', 'SPL'});
[Sdisp.Minutes] = DealElements(round(10*[DB.MinutesSinceRecStart])/10);
[Sdisp.series] = DealElements([DB.iseries]);
[Sdisp.series_run] = DealElements([DB.iseries_run]);
% store uniq rec indices
Uidx = [DB.UniqueRecordingIndex];
%structview(Sdisp);

% construct char matrix for listbox view
[Str, Hd] = struct2char(Sdisp, 8);
Str = [cpref Str];
%Str = [blanks(size(Str,1)).', Str]; % prepend spaces to disable navigation by key presses
Hd = [blanks(size(Hd,1)).', Hd];
% insert separators in display string & in rec indices, so they stay matched 
qsep = diff([Sdisp.series])~=0; % where to insert 
Str = local_sep(Str, qsep, '-');
Uidx = local_sep(Uidx(:), qsep, nan).';



function S = local_cellQuality(ic);
UseTxt = local_xlsfield(ic, 'M', 'Useful');
StabTxt = local_xlsfield(ic, 'J', 'Stability');
BaseTxt = local_xlsfield(ic, 'K', 'Baseline var');
BinTxt = local_xlsfield(ic, 'R', 'Binaurality');
DrivTxt = local_xlsfield(ic, 'Z', 'Driven');
[BC, TC] = local_cell_color(ic);
S = struct('style', 'text', ...
    'fontsize', 10, ...
    'horizontalalign', 'left', ...
    'foregroundcolor', TC, ...
    'backgroundcolor', BC, ...
    'string', [UseTxt char(10) StabTxt char(10) BaseTxt char(10) BinTxt char(10) DrivTxt char(10)]);

