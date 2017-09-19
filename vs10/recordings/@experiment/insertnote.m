function okay = insertnote(E); %#ok<INUSL>
% Experiment/insertnote - insert user note in experiment log file.
%    D = insertnote(E) opens notepad and allows the user to type a note
%    that will be inserted in the log file of experiment E.
%
%    See Datadir, Experiment, experiment/addtolog.

okay = 0;
Instruction = '====insert note below, save and exit====';
FN = tempname;
textwrite(FN, Instruction);
CMD = ['!notepad ' FN ];
eval(CMD);
Note = mytextread(FN);
delete(FN);
if isequal(Instruction, Note{1}),
    Note = Note(2:end);
end
if isempty(Note),
    return;
end
Note = char(Note);
Nline = size(Note, 1);
Note = [repmat('% ', Nline, 1), Note];
Head = ['% =========== note added ' datestr(now) '========='];
Tail = repmat('=', size(Head)); Tail(1:2) = '% ';
Note = strvcat(Head, Note, Tail);
addtolog(E, Note);
okay = 1;





