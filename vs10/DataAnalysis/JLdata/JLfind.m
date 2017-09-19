function [st, ihit] = JLfind(Jb, ST);
% JLfind - find JLbeatStats entry corresponding to single JLbeat output
%    JLfind(Jb, ST);
%
%    See also JLgetBeatStruct.

if nargin<2, ST=[]; end
if isempty(ST),
    % import from caller space
    qqname = ['qq' num2str(RandomInt(1e9))];
    evalin('caller', ['global ' qqname ';']);
    evalin('caller', [qqname '= ST;']);
    eval(['global ' qqname ';']);
    eval(['ST = ' qqname ';']);
    evalin('caller', ['clear global ' qqname ';']);
    eval(['clear global ' qqname ';']);
end

ihit = find([Jb.UniqueRecordingIndex] == [ST.UniqueRecordingIndex]);
st = ST(ihit);

