function S = disp(T)
% tsig/disp - DISP for tsig objects.
%
%   See tsig.

if isvoid(T), Str = 'void tsig object.'; 
else,
    Nc = nchan(T);
    if Nc==1, cstr = 'single';
    elseif Nc==2, cstr = 'double';
    else cstr = num2str(Nc);
    end
    Ns = nsam(T);
    samstr = trimspace(num2str(Ns));
    if length(Ns)>1,
        samstr = ['[' samstr ']'];
    end
    fstr = num2str(0.1*round(10*T.Fsam));
    Rstr = 'real';
    if any(~isreal(T)), 
        Rstr='complex'; 
    end
    if all(haslogic(T)), 
        Rstr='logical'; 
    end
    Str = [cstr '-channel ' Rstr ' time signal containing ' samstr ' samples @ ' fstr ' kHz.'];
    dur = Uniquify(duration(T));
    t0 = Uniquify(starttime(T));
    t1 = Uniquify(endtime(T));
    timestr = [vector2str(dur) '-ms long; starting at ' vector2str(t0,2) ' ms;  ' ...
        'ending at '  vector2str(t1,2)     ' ms.'];
    Str2 = ['   ' timestr];
    Str = strvcat(Str, Str2);
end

if nargout<1,
    disp(Str);
else,
    S = Str;
end









