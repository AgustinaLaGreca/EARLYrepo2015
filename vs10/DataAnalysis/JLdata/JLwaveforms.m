function W = JLwaveforms(Uidx, flag);
% JLwaveforms - retrieve waveforms previously computed by JLviewrec
%    W = JLwaveforms(Uidx)
%    See JL_NNTP for storage conventions.
%
%    W = JLwaveforms(Uidx, 'rmcache')
%    removes cached waveforms and waveform stats. Needed when changing AP
%    threshold criterion etc. Will be restored at next standard call to
%    JLwaveforms.
%
%    See also JLviewrec, JL_NNTP.

if ~isscalar(Uidx), error('Uidx arg must be single number.'); end
GBdir = fullfile(processed_datadir,  '\JL\NNTP');
FFN = fullfile(GBdir, ['waveForms_' num2str(Uidx) '.mat']);
if nargin==2 && isequal('rmcache', flag),
    delete(FFN);
elseif exist(FFN, 'file'),
    load(FFN, 'W');
else, % recompute & save
    load(fullfile(GBdir, 'cycleStats'), 'T');
    ihit = find([T.UniqueRecordingIndex]==Uidx);
    [T(ihit), W] = JLviewrec(Uidx);
    save(FFN, 'W');
    save(fullfile(GBdir, 'cycleStats'), 'T', '-V6');
end









