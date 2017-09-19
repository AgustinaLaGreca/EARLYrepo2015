function I = private_seqPlayInfo(Field, Value);
% private_seqPlayInfo - get/set seqplay info (helper for SeqplayXXX)
persistent PI

if nargin==0, % return all info
    I = PI;
elseif isequal('clear', Field),
    PI = [];
elseif isstruct(Field), % set fields of PI to fields of Field
    FNS = fieldnames(Field);
    FVS = struct2cell(Field);
    for ii=1:length(FNS), 
        private_seqPlayInfo(FNS{ii}, FVS{ii});
    end
elseif nargin==1, % get requested field
    I = PI.(Field);
elseif nargin==2, % set requested field
    PI.(Field) = Value;
end










