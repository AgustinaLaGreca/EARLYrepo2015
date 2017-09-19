function S = disp(T)
% playlist/disp - DISP for playlist objects.
%
%   See playlist.

if isvoid(T), Str = 'void playlist object.'; 
else,
    Nw = nwave(T);
    if Nw==1, wvstr = 'single';
    else,  wvstr = num2str(nwave(T));
    end
    wvstr = [wvstr '-waveform Playlist'];
    Nsam = nsam(T); % sample count of stored waveforms
    Nsam = Nsam(T.iplay); % sample counts of waveforms in list
    Nrep = T.Nrep;
    Nplay = nsamplay(T);
    ip = T.iplay;
    if isempty(ip),
        samstr = '; no play list specified';
    else,
        samstr = [' playing ' num2str(Nrep(1)) 'x' num2str(Nsam(1))];
        for ii=2:length(Nrep),
            samstr = [samstr ' + ' num2str(Nrep(ii)) 'x' num2str(Nsam(ii)) ];
        end
        samstr = [samstr ' = ' num2str(Nplay) ' samples'];
    end
    Str = [wvstr samstr '.'];
end

if nargout<1,
    disp(Str);
else,
    S = Str;
end









