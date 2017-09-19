function SF = stimlist_strfun(D)
% Dataset/stimlist_strfun - helper functions for stimlist char strings.
%   stimlist_strfun(dataset()) returns a struct with function-handle valued
%   fields 
%       xrange: displays range of values like 100:100:300 Hz 
%     shstring: compact string displaying number(s) in A|B|C format
%     num2kstr: format using k for 1000, eg 2k4 == 2400
%       modstr: display modulation like 30-Hz mod, 10|12-Hz mod, no mod
%
%   The signature of these functions is
%    xrange(X), with X the X or Y field of a stimpresentx object
%    shstring(x,sep), with x numeric and sep the separator defaulting to |
%    num2kstr(x), with x numeric
%    modstr(Stim, Prefix), with Stim the stimulus struct of a dataset, and
%            Prefix a prefix of the modulation-specs in Stim.
%    
%
%   See also Dataset/stimlist.

SF = struct('xrange', @local_xr, 'shstring', @local_str, 'num2kstr', @local_num2kstr, 'modstr', @local_mod);

%=====================================
%=====================================
%=====================================
function Str = local_xr(X);
[Xmin, Xmax] = minmax(X.PlotVal);
U = X.ParUnit;
if isequal('octave', lower(U)),
    U = 'Oct';
end
if isequal('Hz', U),
    Str = [local_num2kstr(Xmin) ':' local_num2kstr(Xmax) ' ' U];
else,
    Str = [local_str(Xmin) ':' local_str(Xmax) ' ' U];
end

function s=local_str(x,sep);
if nargin<2, sep='|'; end
x = unique(x);
if numel(x)>1,
    s = [local_str(x(1),sep) sep local_str(x(2:end),sep)];
    return;
end
if abs(x)>100,
    s = num2str(round(x));
else
    s = num2str(deciRound(x,3));
end

function Str = local_mod(Stim, Prefix);
if nargin<2, Prefix=''; end
if isfield(Stim,[Prefix 'ModFreq'])
    modfreq = Stim.([Prefix 'ModFreq']);
else
    modfreq = 0;
end
if isfield(Stim,[Prefix 'ModDepth'])
    moddepth = Stim.([Prefix 'ModDepth']);
else
   moddepth = 0; 
end
if all(0==moddepth.*modfreq), Str = 'no mod'; 
else, Str = [local_str(modfreq) '-Hz mod'];
end

function Str = local_num2kstr(X);
Str = '';
for ii=1:numel(X),
    x = X(ii);
    if abs(x)>= 1000,
        xa = abs(x);
        xk= floor(xa/1000);
        xh= round(rem(xa,1000)/100);
        if xh==10,
            xk=xk+1;
            xh=0;
        end
        s = [num2str(xk) 'k'];
        if xh>0, s= [s num2str(xh)]; end
        if x<0, s = ['-' s]; end
    else,
        s = num2str(round(x));
    end
    Str = [Str s ' '];
end
Str(end) = [];
if ~isempty(strfind(Str,'k')) && isequal('0', Str(end)),
    Str = Str(1:end-1);
end






