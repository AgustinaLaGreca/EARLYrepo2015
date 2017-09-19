function Y = subsref(P, S);
% playlist/SUBSREF - SUBSREF for playsig objects
%
%   P(k) returns the waveform to be played in the kth position in P, 
%   counting repetitions.
%
%   P(k,J) returns the Jth sample(s) of P(k). J may be vector.
%
%   P(K), where K is a vector, returns the requested waveforms concatenated
%   into a single column vector.
%   
%   P.Foo is the the same as Foo(P) provided that Foo is a method for playlist
%   objects.
%
%   P.help produces a list of all fields that may be assessed by subscripted
%   reference of playlist objects. These "fields" are really methods, so
%   type help playlist/Foo for help on field Foo.
%
%   See also playlist, "methods playlist".


if length(S)>1, % use recursion from the left
    y = subsref(P,S(1));
    Y = subsref(y,S(2:end));
    return;
end

%----single-level subsref from here (i.e., S is scalar) ---

switch S.type,
    case '()',
        if isempty(P.iplay),
            error('No list specified for playlist object.');
        end
        idx = S.subs;
        if length(idx)>2,
            error('Invalid indexing of playlist object. Too many indices');
        end
        iwa = iLongWave(P, idx{1}); % expanded play indices; select the ones requested in idx.
        Y = P.Waveform(iwa); % requested waveform(s) in cell array
        Y = cat(1, Y{:}); % same, now pasted into single column array
        if length(idx)>1, Y = Y(idx{2}); end
    case '.',
        Field = S.subs; % field name
        [Method, Mess] = keywordMatch(Field, methods(P), 'playlist field name');
        error(Mess);
        Y = feval(Method, P);
    case '{}',
        error('Playlist objects do not allow subscripted reference using braces {}.');
end % switch/case


