function [P, Mess, varargout] = binop(S,T, operator, defS, defT);
% tsig/binop - binary operator for tsig objects.
%   Helper function for tsig/plus, etc. See those mfiles.
%   Syntax:
%     [P, Mess]  = binop(S,T, operator, defS, defT);
%   defS and defT are the values appended to S and T in time intervals 
%   that are present in either of hem, but not both.
%   List of possible messages and their meanings:
%                 void: a void tsig was passed to binop
%      numnonrowvector: numeric S or T specified that is not a row vector
%                nchan: incompatible channel counts between S and T
%              nontsig: S or T specified that is nor tsig nor numeric.
%                 fsam: S and T have different sample frequencies.
% 
%   See also tsig/plus.

if ~isa(operator, 'function_handle'),
    error('Operator argument must be function handle.');
end
Nxtra = nargout-2; % # additional output args
P =[]; % default output to allow premature exit w/o warnings on missing outarg
if isnumeric(S), 
	Mess = local_testnum(S,T);
    if ~isempty(Mess),
        return; % the calling function will produce a readable error message
    end
    [T.Waveform,S] = SameSize(T.Waveform,S);
    for ichan=1:length(S),
        Xout = cell(1,Nxtra);
        [T.Waveform{ichan}, Xout{:}] = operator(S(ichan),T.Waveform{ichan});
        for iargout=1:Nxtra,
            varargout{iargout}{ichan} = Xout{airgout};
        end
    end
    P = tsig(T.Fsam, T.Waveform, T.t0);
elseif isnumeric(T),
	Mess = local_testnum(T,S);
    if ~isempty(Mess),
        return; % the calling function will produce a readable error message
    end
    [S.Waveform,T] = SameSize(S.Waveform,T);
    for ichan=1:length(T),
        Xout = cell(1,Nxtra);
        [S.Waveform{ichan}, Xout{:}] = operator(S.Waveform{ichan},T(ichan));
        for iargout=1:Nxtra,
            varargout{iargout}{ichan} = Xout{airgout};
        end
    end
    P = tsig(S.Fsam, S.Waveform, S.t0);
elseif ~isTsig(S) || ~isTsig(T),
    Mess = 'nontsig';
    return;
else, % two S & T are both tsig objects
    % first aligh S and T, i.e. equalize channel counts & provide time 
    % domain in either that is missing in the other.
    [S,T,Mess] = align(S,T, defS, defT); % newly added S has defS value; dito T & defT
    if ~isempty(Mess),
        return; % the calling function will produce a readable error message
    end
    % now that S&T have the same format, and the default values have been ..
    % put in place, we can simply apply the operator channelwise.
    P = S;
    for ichan=1:nchan(S),
        Xout = cell(1,Nxtra);
        [P.Waveform{ichan}, Xout{:}] = operator(S.Waveform{ichan},T.Waveform{ichan});
        for iargout=1:Nxtra,
            varargout{iargout}{ichan} = Xout{airgout};
        end
    end
    [P, Mess] = test(P); % check consistency
    error(Mess);
end

%-------------------------------
function Mess = local_testnum(X,T);
% test compatibility of numerical X and tsig T
Mess = '';
if isvoid(T),
    Mess = 'void';
elseif ~isvector(X),
    Mess = 'numnonrowvector';
elseif ~isequal(1,size(X,1)),
    Mess = 'numnonrowvector';
elseif length(unique([1 length(X) nchan(T)]))>2, % incompatible channel counts
    Mess = 'nchan';
end







