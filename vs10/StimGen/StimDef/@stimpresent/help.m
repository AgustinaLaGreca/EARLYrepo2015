function y=help(Q,Nindent) %#ok<INUSL>
% stimpresent/help - help for stimpresent objects
%    help(S) or S.help displays the help text for stimpresent objects.
%    S is a dummy stimpresent object allowing help to be a stimpresent
%    method.
%
%    help(S,N) uses N blanks for indentation.
%
%    See also stimpresent, methods.

if nargin<2, Nindent=0; end
B = blanks(Nindent+3);
y=invisible; % dummy  output
disp([B '       Fsam:  sample rate [Hz]']);
disp([B '      Npres: number of stimulus presentations, reps counted']);
disp([B '      Ncond: number of different conditions']);
disp([B '       Nrep: number of repetitions per condition']);
disp([B 'TotNsamPlay: total # samples in stimulus']);
disp([B '     TotDur: total stim duration [ms]']);
disp([B '-----------------------------------']);
disp([B '      iCond: index array; iCond(k)==j -> k-th pres plays j-th condition']);
disp([B '       iRep: index array; iRep(k)==j -> k-th pres is the j-th rep of that condition']);
disp([B '   NsamPres: sample counts; NsamPres(k)=n -> k-th pres plays n samples']);
disp([B '  SamOffset: sample offsets of pres since start of stim']);
disp([B '  PresOnset: time offsets [ms] of pres since start of stim']);
disp([B '    PresDur: pres durations [ms]']);
disp([B '-----------------------------------']);
disp([B '          X: details on 1st varied stimulus param [struct]']);
disp([B '          Y: details on 2nd varied stimulus param [struct]']);






