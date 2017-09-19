function Y = imagZWOAErecording(X,irec,comp_irec);
% imagZWOAErecording - storage & retrieval of "imaginary" ZWOAE recordings
%   private helper fcn for poolZWOAEdata and ZWOAEimport.
%       iRec = imagZWOAErecording(S) stores S and returns imag index iRec.
%       S = imagZWOAErecording(iGerbil, iRec) retrieves S with imag index iRec.
%       S = imagZWOAErecording(iGerbil, iRec, iComp) sets iRec's companion to iComp

if nargin<3, comp_irec=[]; end

persistent TTT BK
N = numel(TTT);

if isstruct(X), % store
    % check if already stored
    UID = {X.ExpID X.recID}; % unique ID
    for ii=1:N,
        if isequal(BK{ii}, UID), % yes, stored before. return i*index
            Y = i*ii;
            return;
        end
    end
    % if we get here, it's a new one
    TTT = [TTT X];
    BK{N+1} = UID;
    Y = i*(N+1);
elseif ~isempty(comp_irec), % change comp of stored irec
    index = imag(irec);
    TTT(index).companionID = comp_irec;
else, % retrieve
    index = imag(irec);
    if index>N,
        error(['No pooled ZWOAE data available having index i*' num2str(index) '.']);
    end
    Y = TTT(index);
end




