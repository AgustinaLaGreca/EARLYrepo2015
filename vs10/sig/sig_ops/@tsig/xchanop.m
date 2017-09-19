function [Y, Mess] = xchanop(S, operator, defS);
% tsig/xchanop - across-channel function evaluation.
%   Helper function for tsig/median, etc. See those mfiles.
%   Syntax:
%    [P, Mess]  = xchanop(S, @Fun, defS)
%    or [P, Mess]  = xchanop(S, @fun, 'existingchannelsonly');
%
%    xchanop computes a sample-wise Fun([s1,s2,.]), where s1, s2 .. are the
%    simultaneous samples of channels S(1), S(2) .. .
%
%    The crucial thing is the time alignment. For those time instances
%    where some channels exist, and others don't, either a default
%    value DefS is substituted, or only those channels are passed to Fun which
%    exists at that time. The latter mode ('existingchannelsonly' option)
%    is done in blockmode using tsig/exist, so it is fast.
%
%    Note: use an anonymous function when you need to pass additional args
%    to Fun: myfun = @(x)fun(x,par1,..), and pass myfun to xchanop.
%
%   See also tsig/binop, tsig/relop.

if ~isa(operator, 'function_handle'),
    error('Operator argument must be function handle.');
end
if ~isTsig(S),
    Mess = 'nontsig';
    return;
end
Y = []; % allow premature exit

% time alignment: construct zero-valued signal having max time range
t0=min(onset(S)); t1=max(offset(S));
S0 = zeroTsig(Fsam(S),t1-t0,t0); % S0 has full time range
if ~ischar(defS), % all channels participate, even outside their region of existence
    [S, Mess] = binop(S,S0,@plus,defS,0); % this time-aligns all channels, providing defS values where needed
    if ~isempty(Mess), return; end
    Y = operator([S.Waveform{:}]); % dump waveforms in big matrix; compute median
else, 
    [defS, Mess] = keywordMatch(defS,{'existingchannelsonly'}, 'defS string');
    error(Mess); % this is an error from the calling function, not the "end user"
    switch  defS,
        case 'existingchannelsonly', % what else ;)
            % more tricky; we need to track the regions of existence of each of the signals
            % evaluate the subsequent intervals ("Episodes") in which ...
            % ... different subsets of the channels exist.
            Episodes = exist(S); %#ok<EXIST> MLint does not understand that exist is overloaded
            W = Waveform(S+S0); % adding S0 fills the gaps with zeros, so we can access the ... 
            W = [W{:}]; % ... samples in a matrix
            Y = zeros(1,Nsam(S0));
            for epi=Episodes, % run median on this episode using the subset of channels that cover it (see tsig/exist)
                Y(1,epi.istart:epi.iend) = operator(W(epi.istart:epi.iend, epi.ichan));
            end
        otherwise, error(['Invalid defS char value ''' defS '''.']);
    end
end
S0.Waveform = Y;
[Y, Mess] = test(S0);
error(Mess);



