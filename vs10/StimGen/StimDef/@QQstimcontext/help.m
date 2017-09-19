function y=help(Q,Nindent)
% stimcontext/help - help for StimContext objects
%    help(S) or S.help displays the help text for StimContext objects.
%    S is a dummy stimContext object allowing help to be a sytimContext
%    method.
%
%    help(S,N) uses N blanks for indentation.
%
%    See also StimContext, methods.
if nargin<2, Nindent=0; end
B = blanks(Nindent+3);
y=invisible; % dummy  output
disp([B '        DAchan: active DA channel(s) [char]    Left|Right|Both']);
disp([B '       minFreq: minimum stimulus frequency [Hz]']);
disp([B '       maxFreq: maximum stimulus frequency [Hz]']);
disp([B '      maxNcond: max # conditions in a recording']);
disp([B '       recSide: recording side [char]   Left|Right'])
disp([B ' ITDconvention: meaning of ITD>0 [char] IpsiLeading|ContraLeading|LeftLeading|RightLeading']);
disp([B '      Hardware: hardwareSpec object']); %help(HardwareSpec, Nindent+3);
disp([B '   Preferences: EarlyPref object']);    %help(EarlyPref, Nindent+3);






