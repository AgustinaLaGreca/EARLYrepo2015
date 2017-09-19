function [irecAC, irecCM]=ZWOAEcompanions(igerbil, irec, flag);
% ZWOAEcompanions - find AC/CM pair of ZWOAE recordings
%    [irecAC, irecCM]=ZWOAEcompanions(iGerbil, irec) returns the data
%    indices [irecAC, irecCM] of a pair AC and CM recordings, one of which
%    equals the input index irec. iGerbil is the experiment index.
%
%    ZWOAEcompanions(iGerbil, irec, 'unique') uniquifies the outputs. Note
%    that this sorts the output arguments.
%
%    See also ZWOAEimport.

if nargin<3, flag = ''; end

doUniq = isequal('unique', flag);
%handle multiple datasets recursively
if numel(irec)>1,
    for ii=1:numel(irec),
        [irecAC(ii), irecCM(ii)]=ZWOAEcompanions(igerbil, irec(ii));
    end
    if doUniq,
        [irecAC, iuniq] = unique(irecAC);
        irecCM = irecCM(iuniq);
    end
    return;
end

%-----single irec from here-----------
qq = ZWOAEimport(igerbil,irec,'-nosig');
icomp = qq.companionID;
if ischar(icomp), error(icomp); end
rType = qq.RecType;
if isempty(icomp) | isequal(0,icomp), % try to find it "manually" (note: only for gerbils 70,71,72)
    if isequal('ACOUSTIC',rType),
        icomp = irec+1;
    elseif isequal('CM',rType),
        icomp = irec-1;
    else,
        error('No companion recording found; recording type is stored in data.');
    end
    qq2 = ZWOAEimport(igerbil,icomp,'-nosig');
    % check if it's indeed a companion
    qq = rmfield(qq,{'RecType' 'recID' });
    qq2 = rmfield(qq2,{'RecType' 'recID'});
    if ~isequal(qq,qq2),
        error(['No companion recording found; candidate companion irec=' num2str(icomp) ' has different stimulus parameters.']);
    end
end
switch rType,
    case 'ACOUSTIC',
        [irecAC, irecCM] = deal(irec, icomp);
    case 'CM',
        [irecAC, irecCM] = deal(icomp, irec);
end









