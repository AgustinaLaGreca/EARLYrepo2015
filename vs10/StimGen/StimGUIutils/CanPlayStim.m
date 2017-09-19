function cp=CanPlayStim;
% CanPlayStim - true if hardware for auditory stimulation is present
%    CanPlayStim returns True (1) if hardware for auditory stimuli is
%    present. Currently, this is equivalent to a non-empty sys3devicelist,
%    but that might change one day ;)

cp = ~isempty(sys3devicelist);


