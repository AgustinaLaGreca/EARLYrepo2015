function okay=datasave(D,del);
% Dataset/datasave - temporarily save Dataset object.
%
%   datasave(D,'delete') removes all files of temporarily saved
%
%   See Dataset, dataset/save.

okay = true; % optimistic default;

if nargin<2, del = []; end

EE = D.ID.Experiment;

% never compete with the final version of D, so make sure...
fp = [filename(EE) '_DS' zeropadint2str(nan,5) '.EarlyDS']; %... iDataset is nan

try
    if ~isempty(del)
        delete(fp);
    else
        warning('off','MATLAB:unassignedOutputs');
        save(fp, 'D', '-mat'); % save D in MAT format
        warning('on','MATLAB:unassignedOutputs');
    end
catch IOE % possible I/O error
    disp(IOE.getReport());
    okay = false;
end


