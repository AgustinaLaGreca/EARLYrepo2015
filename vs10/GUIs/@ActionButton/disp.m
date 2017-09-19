function Str=disp(Q)
% ActionButton/disp - DISP for ActionButton objects.
%
%   See ParamQuery.

if numel(Q)>1,
    for ii=1:numel(Q),
        if nargout>0, Str{ii}=disp(Q(ii));
        else, disp(Q(ii));
        end
    end
    return
end
% single element from here

if isempty(Q.Name),
    Str = 'Void ActionButton object';
elseif ischar(Q.String),
    Str = ['ActionButton  ' Q.Name ' labeled ' Q.String  ';   Callback = ' char(Q.Callback')];
else, % cellstr
    StrStr = ['{' Q.String{1}]; 
    for ii=2:length(Q.String), StrStr=[StrStr '|' Q.String{ii}]; end
    StrStr = [StrStr '} showing ' Q.CurrentString];
    Str = ['ActionButton  ' Q.Name ' labeled ' StrStr ';   Callback = ' char(Q.Callback')];
end

if nargout<1,
    disp(Str);
    clear Str
end






