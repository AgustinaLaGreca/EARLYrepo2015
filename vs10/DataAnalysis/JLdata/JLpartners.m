function Y = JLpartners(S, ST);
% JLpartners - find monaural or binaural partners of JL recording
%    Y = JLpartners(S, ST); returns struct with fields
%    Ipsi, Contra, Bin. Each is an array of ST indices pointing to
%    Ipsi, Contra and Binaural recordings, respectively, from the same cell
%    which share with S the stimulation frequencies and the SPL.
% 
%    If the ST  argument is omitted, JLpartners looks ST in mat file
%    STall.mat in the JL\JLbeatStats\ subdir of processed_datadir.
%
%    See also JLbeatStats, processed_datadir.

persistent last_loaded
if nargin<2,
    if isempty(last_loaded), % load from file
        FN = fullfile(processed_datadir, 'JL', 'JLbeatStats', 'STall.mat');
        if ~exist(FN, 'file');
            error('Missing input arg ST and missing STall.mat  file.');
        end
        load(FN);
        last_loaded.ST = ST;
    end
    ST = last_loaded.ST;
end

UCI = S.UniqueCellIndex;
switch S.StimType,
    case 'I',
        Y.Ipsi = local_match_fields(ST, UCI, 'I', S.Freq1, [], S.SPL1, []);
        Y.Contra = local_match_fields(ST, UCI, 'C', [], S.Freq2, [], S.SPL1);
        Y.Bin = local_match_fields(ST, UCI, 'B', S.Freq1, S.Freq2, S.SPL1, S.SPL1);
    case 'C',
        Y.Ipsi = local_match_fields(ST, UCI, 'I', S.Freq1, [], S.SPL2, []);
        Y.Contra = local_match_fields(ST, UCI, 'C', [], S.Freq2, [], S.SPL2);
        Y.Bin = local_match_fields(ST, UCI, 'B', S.Freq1, S.Freq2, S.SPL2, S.SPL2);
    case 'B',
        Y.Ipsi = local_match_fields(ST, UCI, 'I', S.Freq1, [], S.SPL1, []);
        Y.Contra = local_match_fields(ST, UCI, 'C', [], S.Freq2, [], S.SPL2);
        Y.Bin = local_match_fields(ST, UCI, 'B', S.Freq1, S.Freq2, S.SPL1, S.SPL2);
end


function iok = local_match_fields(ST, UCI, StimType, Freq1, Freq2, SPL1, SPL2);
iok = [ST.StimType]==StimType;
iok = iok & [ST.UniqueCellIndex]==UCI;
if ~isempty(Freq1),
    iok = iok & [ST.Freq1]==Freq1;
end
if ~isempty(Freq2),
    iok = iok & [ST.Freq2]==Freq2;
end
if ~isempty(SPL1),
    iok = iok & [ST.SPL1]==SPL1;
end
if ~isempty(SPL2),
    iok = iok & [ST.SPL2]==SPL2;
end
iok = find(iok);








