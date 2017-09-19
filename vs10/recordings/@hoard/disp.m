function Str = disp(H);
% hoard/disp - DISP for hoard object
%    disp(H) displays hoard object H
%
%    See also hoard/set, hoard/cat.

nam = get(H.figh, 'tag');
d = get(H.h, 'userdata');
Sstr = num2str(numel(H));;
Str = ['hoard object named ''' nam ''' containing ' Sstr ' samples.'];
if nargout<1,
    disp(Str);
    clear Str
end







