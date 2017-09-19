function R=prepare(R);
% playsound/prepare - prepare playsound object for action
%    prepare(actP) prepares playsound object actP for doing its job,
%    delivering a sound stimulus.
%
%    This amounts to preparing the hardware for D/A conversion. The info 
%    needed to do this has been provided at constructor time. Note that
%    playsound/start is overloaded, because the action is not really
%    carried out by the timer, but is started once by calling a seqplaygo
%    which triggers the hardware. The timer still has the important task of
%    regularly testing whether actP is ready (see action/isready), and if
%    so, of changing the status of actP to 'finished'.
%
%  See also makestimFS, sortConditions, action, dynamic.

persistent LastWavID 

eval(IamAt);
if ~isequal('initialized', status(R)),
    error('Can only prepare object whose status is ''initialized''.');
end

R = download(R); % get the latest version

R.Wav = download(R.Wav); 
error(TestNsam(R.Wav));
% due to Matlab's overloaded subsref's failure to return multiple argouts, 
% we need to convert the waveforms to a struct
R.Wav = struct(R.Wav); % must be reversed @ end (before final uploading of R)
if isequal(LastWavID, [R.Wav(:).UniqueID]),
    % exactly the same set of waveforms was passed last time; assume the 
    % previous prepare is still valid; skip it now, except for the
    % attenuator settings
    SetAttenuators(R.Experiment, R.Atten.AnaAtten);
    R = status(R, 'prepared'); 
    return;
end

%=======the following code is too specific (TDT!) to be here. Needs cleaning.
SetAttenuators(R.Experiment,'max'); % set attenuators to their maximum value
% Initialize RPvdS circuit
Dev = R.Experiment.AudioDevice;
[CI, Recycled]=seqplayinit(R.Wav(1).Fsam/1e3, Dev, circpostfix(R.Rec.CircuitInstr.RX6_circuit), 1); % #noload: circuit has been loaded already; 1 means recycle circuit if possible
% if ~Recycled,
%     error('Correct Seqplay circuit was not loaded prior to
%     playsound/prepare call');
% end

% Expand the play list: W has multiple conditions, each containing multiple
% segments. SeqplayList does not care, it only needs to know the segment
% indices and Nreps in the grand pool of waveform segments.
[iWav, Nrep, ConditionIndex, Wav]=local_expand_Playlist(R.StimPres, R.Wav);
[Triggers] = local_get_triggers(iWav, Nrep, ConditionIndex, Wav);

% upload waveforms
[Ncond, Nchan] = size(Wav);
Scale = R.Atten.NumScale; % scale factors per element of W
Scale = local_expand_scaleFactors(Scale, Wav); % now per segment making up the elements of W

if Nchan==1,
    switch Wav(1).DAchan,
        case 'L',
            seqplayupload([Wav(:,1).Samples], {}, Scale{1});
            seqplaylist(iWav{1}, Nrep{1});
        case 'R',
            seqplayupload({}, [Wav(:,1).Samples], 1, Scale{1});
            seqplaylist([],[], iWav{1}, Nrep{1});
    end
else, % dual channel play
    seqplayupload([Wav(:,1).Samples], [Wav(:,2).Samples], Scale{:});
    seqplaylist(iWav{1}, Nrep{1}, iWav{2}, Nrep{2});
end

% Upload trigger signals
% % % % % % sys3setpar(length(Triggers),'BaseLineBuffSize',Dev);
% % % % % % sys3write(Triggers, 'BaselineIndices', Dev,0,'I32');

% set analog attenuators
SetAttenuators(R.Experiment, R.Atten.AnaAtten);
LastWavID = [R.Wav(:).UniqueID]; % store to maybe save tine on replay

R = status(R, 'prepared'); % note: uploading is implicit in this call, too
end

%============================================================
%============================================================
function [iWav, Nrep, iCond, W]=local_expand_Playlist(Pres, W);
Nchan = size(W,2);
RS = cellify(Pres.ReducedStorage);
if numel(RS) < Nchan, RS(1:Nchan) = RS(1); end; 
for ichan = 1:size(W,2),
    iwaveOffset = cumsum([0 W(1:end-1,ichan).Nwav]); % offsets per condition in the order of storage in W (and as uploaded)
    iwaveOffset = iwaveOffset(Pres.iCond); % ditto, but now in the order of presentation
    %disp('========')
    Nwav = [W(Pres.iCond,ichan).Nwav]; % wave counts of each presentation
    iw = arrayfun(@(v,n)(v+(1:n)), iwaveOffset(:), Nwav(:), 'UniformOutput', false);
    ic = arrayfun(@(c,n)(c*ones(1,n)), Pres.iCond(:), Nwav(:), 'UniformOutput', false);
    iWav{ichan} = [iw{:}];
    % added: reduce storage by removing duplicate segments
    [iWav, W] = local_reduce_Playlist(W, iWav, ichan, RS{ichan});
    iCond{ichan} = [ic{:}];
    Nrep{ichan} = [W(Pres.iCond,ichan).Nrep];
end
end

function S = local_expand_scaleFactors(Scale, W);
    for ichan = 1:size(W,2),
        Nwav = [W(:,ichan).Nwav]; % wave counts of each condition
        qq = arrayfun(@(c,n)(c*ones(1,n)), Scale(:,ichan).', Nwav, 'UniformOutput', false); 
        S{ichan} = [qq{:}];
    end
end
function wav = local_struct_to_waveform(wave_struct)
    for i=1:length(wave_struct)
        wav(i) = Waveform(wave_struct(i).Fsam, wave_struct(i).DAchan, wave_struct(i).MaxMagSam, ...
            wave_struct(i).SPL, wave_struct(i).Param, wave_struct(i).Samples, wave_struct(i).Nrep);
    end

end

function [iWav,W] = local_reduce_Playlist(W, iWav, ichan, crit)
switch crit
    case ''
        return;
    case 'nonzero'
        if size(W,1)<4, return; end % no more reduction possible
        new = find(cellfun(@(X)(numel(find(X))),W(1,ichan).Samples)); % use 1st cond as ref
        W(1,ichan).Nonzero = new;
        for icond=2:size(W,1)-2 % no baseline
            findNonzero = find(cellfun(@(X)(numel(find(X))),W(icond,ichan).Samples));
            W(icond,ichan).Samples(:,findNonzero) = [];
            W(icond,ichan).Nonzero = findNonzero;
        end
        NNonzero = cellfun(@(X)(numel(X)), {W(2:end-2,ichan).Nonzero}); % no baseline
        if ~all(NNonzero==NNonzero(1)), error('Playlist could not be reduced. ''nonzero'' criterion is not applicable.')'; end    
        Nonzero = [W(2:end,ichan).Nonzero];
        iwaveOffset = [cumsum([W(1:end-3,ichan).Nwav])]; % no baseline
        iwaveOffset = ones(NNonzero(1),1)*iwaveOffset;
        iwaveOffset = iwaveOffset(:)';
        inz = arrayfun(@(v,n)(v+n), iwaveOffset(:), Nonzero(:), 'UniformOutput', false);
        iNonzero = [inz{:}];
        iw = iWav{1,ichan};
        for i=1:numel(iw)
            if ismember(iw(i),iNonzero)
                iw(i) = new(mod(find(iw(i)==iNonzero),NNonzero(1))+1); % replace copies
            end
        end
        Uind = unique(iw);
        for ii=Uind, iw(iw==ii) = find(Uind==ii); end % lower other indices
        iWav{1,ichan} = iw; % store

end
end

function [Triggers] = local_get_triggers(iWavs, Nreps, iConds, Wav)

    Nchan = size(Wav,2);
    Triggers = [];
    % Execute Loop For Each Channel
    for chan=1:Nchan
        iCond = iConds{1};
        Nrep = Nreps{1};
        iWav = iWavs{1};
        iTrigger = 1;
        waveform_order(1) = iCond(1);
        total_samples = 40;
        
        for i=1:length(iCond)
            if waveform_order(end) ~= iCond(i)
                waveform_order(end + 1) = iCond(i);
            end
        end
        
        for i=1:length(waveform_order)
            W = Wav(waveform_order(i),1);
            rep = W.Annotations.Repetitions;
            for j=1:rep
                if(W.Annotations.BeginOfBaseline == 0)
                    total_samples = total_samples + W.NsamPlay;
                else
                    Triggers(iTrigger) = total_samples + W.Annotations.BeginOfBaseline;
                    iTrigger = iTrigger + 1;
                    for k=1:10
                       Triggers(end+1) = Triggers(end)+5;  
                       iTrigger = iTrigger + 1;
                    end
                    total_samples = total_samples + W.NsamPlay;
                    
                end
            end
        end
    end

end

% % % function [Triggers] = local_get_triggers(iWavs, Nreps, iConds, Wav)
% % % 
% % %     Nchan = size(Wav,2);
% % %     Triggers = [];
% % %     Execute Loop For Each Channel
% % %     for chan=1:Nchan
% % %         iCond = iConds{chan};
% % %         Nrep = Nreps{chan};
% % %         iWav = iWavs{chan};
% % %         iTrigger = 1;
% % %         prev_label = 'stim';
% % %         prev_cond = -1;
% % %         sample_index = 0;
% % %         index_in_waveform = 1;
% % %         nWaves = -1;
% % %         
% % %         for i=1:length(iWav)
% % %             Check if new waveform, if not so check if last Sample array in
% % %             waveform to update the prev_cond
% % %             if iCond(i) == prev_cond
% % %                 index_in_waveform = index_in_waveform + 1;
% % %                 if index_in_waveform > nWaves
% % %                     index_in_waveform = 1;
% % %                 end
% % %                 
% % %                 
% % %                 Check if it is baseline or stim
% % %                 if strcmpi(Wav(iCond(i),chan).Annotations.Label(index_in_waveform),'stim');
% % %                     len = Wav(iCond(i),chan).Annotations.Length(index_in_waveform);
% % %                     rep = Wav(iCond(i),chan).Annotations.Nrep(index_in_waveform);
% % %                     sample_index = sample_index + len;
% % %                     prev_label = 'stim';
% % %                 else
% % %                     if strcmpi(prev_label,'stim')
% % %                         Triggers(iTrigger) = sample_index+1;
% % %                         iTrigger = iTrigger + 1;
% % %                     end
% % %                     len = Wav(iCond(i),chan).Annotations.Length(index_in_waveform);
% % %                     rep = Wav(iCond(i),chan).Annotations.Nrep(index_in_waveform);
% % %                     sample_index = sample_index + len;
% % %                     prev_label = 'baseline';
% % %                 end    
% % %                 
% % %                 prev_cond = iCond(i);
% % %             else
% % %                 index_in_waveform = 1;
% % %                 nWaves = length(Wav(iCond(i),chan).Annotations.Label);
% % %                 Check if it is baseline or stim
% % %                 if strcmpi(Wav(iCond(i),chan).Annotations.Label(index_in_waveform),'stim');
% % %                     len = Wav(iCond(i),chan).Annotations.Length(index_in_waveform);
% % %                     rep = Wav(iCond(i),chan).Annotations.Nrep(index_in_waveform);
% % %                     sample_index = sample_index + len;
% % %                     prev_label = 'stim';
% % %                 else
% % %                     if strcmpi(prev_label,'stim')
% % %                         Triggers(iTrigger) = sample_index+1;
% % %                         iTrigger = iTrigger + 1;
% % %                     end
% % %                     len = Wav(iCond(i),chan).Annotations.Length(index_in_waveform);
% % %                     rep = Wav(iCond(i),chan).Annotations.Nrep(index_in_waveform);
% % %                     sample_index = sample_index + len;
% % %                     prev_label = 'baseline';
% % %                 end
% % %                 
% % %                 prev_cond = iCond(i);
% % %             end
% % %         end
% % %     end
% % % 
% % % end