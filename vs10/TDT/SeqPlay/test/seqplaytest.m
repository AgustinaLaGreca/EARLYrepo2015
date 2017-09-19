function testSeqPlay(wav);
% test function for seq play

Nwav = length(wav);
Lwav = cellfun(@length, wav);
iswitch = cumsum([2 Lwav]); % switch events are determined by start/end of waveforms
Nsam = max(iswitch+1);
ioffset = [1:Nwav -Nsam]; % switch events are determined by start/end of waveforms

sys3loadCircuit('testseqplay');

% collect waveforms in single column vector
Grand_wav = [zeros(2,1)];
for iwav=1:Nwav,
    Grand_wav = [Grand_wav; wav{iwav}(:)]; % 
end

local_loadList(iswitch, ioffset);

sys3write(Grand_wav, 'Waveforms', '', 0, 'F32');

sys3setpar(-1,'EndSample'); 
NextSwitch = sys3getpar('NextSwitch')


% start the ticking clock by setting Endsample to positive value
sys3run;
NextSwitch = sys3getpar('NextSwitch')
% reset all buffers
sys3trig(1);
NextSwitch = sys3getpar('NextSwitch')
inow = sys3getpar('isampnow')

% start playing
sys3setpar(Nsam,'EndSample'); sys3getpar('Tick');
pause(1);
TickMon = sys3read('TickMon', 'Ntick','', 0, 'I32'); iSamp = sys3read('iSamp', 'Nsamp','', 0, 'I32'); plot(TickMon, iSamp-TickMon, '.'); ff
% plot(TickMon, iSamp-TickMon, '.');

%================
function local_loadList(iswitch, ioffset);
iswitch = iswitch(:).';
ioffset = ioffset(:).';
if any(diff([0 iswitch])<2), 
    error('Switch events must be at least 2 samples apart.');
end
if ~isequal(length(iswitch), length(ioffset)),
    error('Switch list and Offset list must have equal length');
end
sys3write([0 -1  iswitch-1], 'SwitchList', '', 0, 'I32'); 
sys3write([0 0 0 ioffset], 'OffsetList', '', 0, 'I32'); 



