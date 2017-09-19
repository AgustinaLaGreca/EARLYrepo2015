function S = JLfindpartners(S);
% JLfindpartners - find partner recordings that share stimulus freq & SPL
%    S = JLfindpartners(S) adds the following fields the the output S of
%    JLvarstats
%      I_partners: array of UniqueRecordingIndex values of all members of
%          S containing ipsi-only recordings from the same cell and
%          obtained with identical Freq1 and SPL1 values.
%      C_partners: array of UniqueRecordingIndex values of all members of
%          S containing contra-only recordings from the same cell and
%          obtained with identical Freq2 and SPL2 values.
%      B_partners:  array of UniqueRecordingIndex values of all members of
%          S containing binaural recordings from the same cell and
%          obtained with identical Freq1 and SPL1 (or Freq2 and SPL2 for
%          contra recordings) values.

S = local_addfields(S);

MatchingFields = {'Freq1' 'Freq2' };

CellIndex = [S.UniqueCellIndex];
unCellIndex = unique(CellIndex);
Ncell = numel(unCellIndex);
for icell=1:Ncell,
    icell
    cell_index = unCellIndex(icell);
    RecIndex = find([S.UniqueCellIndex]==cell_index);
    for irec=RecIndex(:).',
        S(irec).I_partners = local_match_I(S, S(irec));
        S(irec).C_partners = local_match_C(S, S(irec));
        S(irec).B_partners = local_match_B(S, S(irec));
    end
end

%==================================================
function S = local_addfields(S);
[S.partners________________] = deal('_____________');
[S.I_partners] = deal([]);
[S.C_partners] = deal([]);
[S.B_partners] = deal([]);

function idx = local_match_I(S, Sref);
% hits must be from same cell, have ipsi stim and matching Freq1
qhit =  ([S.UniqueCellIndex]==[Sref.UniqueCellIndex]) ...
    & ([S.StimType]=='I') & ([S.Freq1]==[Sref.Freq1]);
% which freq and SPL to compare depends on ref stim
switch Sref.StimType,
    case 'I', qhit = qhit & ([S.SPL1]==[Sref.SPL1]);
    case 'C', qhit = qhit & ([S.SPL1]==[Sref.SPL2]);
    case 'B', qhit = qhit & ([S.SPL1]==[Sref.SPL1]);
end
idx = [S(qhit).UniqueRecordingIndex];

function idx = local_match_C(S, Sref);
% hits must be from same cell, have contra stim and matching Freq2
qhit =  ([S.UniqueCellIndex]==[Sref.UniqueCellIndex]) ...
    & ([S.StimType]=='C') & ([S.Freq2]==[Sref.Freq2]);
% which freq and SPL to compare depends on ref stim
switch Sref.StimType,
    case 'I', qhit = qhit & ([S.SPL2]==[Sref.SPL1]);
    case 'C', qhit = qhit & ([S.SPL2]==[Sref.SPL2]);
    case 'B', qhit = qhit & ([S.SPL2]==[Sref.SPL2]);
end
idx = [S(qhit).UniqueRecordingIndex];

function idx = local_match_B(S, Sref);
% hits must be from same cell, have bin stim and matching Freq1, Freq2
qhit =  ([S.UniqueCellIndex]==[Sref.UniqueCellIndex]) ...
    & ([S.StimType]=='B') & ([S.Freq1]==[Sref.Freq1]) & ([S.Freq2]==[Sref.Freq2]);
% which freq and SPL to compare depends on ref stim
switch Sref.StimType,
    case 'I', qhit = qhit & ([S.SPL1]==[Sref.SPL1]);
    case 'C', qhit = qhit & ([S.SPL2]==[Sref.SPL2]);
    case 'B', qhit = qhit & ([S.SPL1]==[Sref.SPL1]) & ([S.SPL2]==[Sref.SPL2]);
end
idx = [S(qhit).UniqueRecordingIndex];








