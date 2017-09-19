function Str = disp(D);
% Dataset/disp - DISP for Dataset objects.
%
%   See Dataset.

if numel(D)>1, % recursion
    for ii=1:numel(D),
        disp(D(ii));
    end
    return
end

% ======single element from here=======
% if ~D.isstatic,
%     try, D = download(D); end;
% end
if isempty(D),
    [M,N] = size(D);
    disp([num2str(M) ' x ' num2str(N) ' Dataset object']);
    return
elseif isvoid(D),
    disp('Void Dataset object');
    return;
elseif isdummy(D),
    disp(['Dummy dataset ' expname(D) ' Rec ' num2str(irec(D))]);
    return;
end
%=====non-empty, non-void D from here=======
id = D.ID;
ID = [IDstring(D, 'full') ' Pen ' num2str(id.iPen) ', ' num2str(id.PenDepth) ' um  ' id.Location '/' id.Computer ' ' id.modified];
Pres = D.Stim.Presentation;

[xmin, xmax] = minmax(Pres.X.PlotVal(:)); 
Ndec = 2+abs(log10(2*abs(xmax-xmin)/abs(xmax+xmin))); Ndec = max(3,round(Ndec));
[xmin, xmax] = deal(deciRound(xmin, Ndec), deciRound(xmax,Ndec));
if isequal('Linear', Pres.X.PlotScale), StepStr = ' (lin steps)'; else, StepStr = ' (log steps)'; end;
ParStr = ['   ' Pres.X.ParName ' = ' num2str(xmin) ':' num2str(xmax) ' ' Pres.X.ParUnit StepStr];
if has2varparams(D),
    [ymin, ymax] = minmax(Pres.Y.PlotVal(:)); [ymin, ymax] = deal(deciRound(ymin), deciRound(ymax));
    if isequal('Linear', Pres.Y.PlotScale), StepStr = ' (lin steps)'; else, StepStr = ' (log steps)'; end;
    ParStr = strvcat(ParStr, ...
        ['   ' Pres.Y.ParName ' = ' num2str(ymin) ':' num2str(ymax) ' ' Pres.Y.ParUnit StepStr]);
end

StimStr = ['   ' D.Stim.StimType ' stim   ' durationstring(D)];
if isempty(D.Stopped),
    StimStr = [StimStr '  (' D.Status ')'];
else,
    Npr = NpresRecorded(D);
    Ncond = numel(unique(Pres.iCond(1:Npr)));
    Nrep = numel(unique(Pres.iRep(1:Npr)));
    StimStr = [StimStr '  *** STOPPED at pres. # ' num2str(D.Stopped.ipres) ' '];
    StimStr = [StimStr '(' num2str(Ncond) ' conditions, '  num2str(Nrep) ' reps) ***'];
end

RecStr = '';

DFNS = fieldnames(D.Data);
for ii=1:numel(DFNS),
    dfn = DFNS{ii};
    d = D.Data.(dfn);
    if isstruct(d)
        RecStr = 'Obsolete data format.'; 
        break; 
    end
    
    %by abel, not implemented yet for THR data
    if strcmpi(dfn, 'Thr')
        RecStr = 'Thr not inplemented yet in disp function'; 
        break;
    end
    
    rstr = sprintf('   %15s: ' , dataType(d));
    switch class(D.Data.(dfn)),
        case 'adc_data', 
            adc = D.Data.(dfn);
            Dur = 1e3*Nsam(adc)/Fsam(adc);
            rstr = [rstr  sprintf('%6d ms @ %0.1f kHz (%i samples)',round(Dur), deciRound(Fsam(adc)/1e3, 4), Nsam(adc))];
        case 'eventdata',
            ed = D.Data.(dfn);
            rstr = [rstr sprintf('%i spikes', Nevent(ed))];
        otherwise,
    end
    RecStr = strvcat(RecStr, rstr);
    
end

% [Pres.X.ParName ' = ' num2str([ max(Pres.X.PlotVal)]))]
%DD.DataStatus = grabstatus(D);
disp(ID);
disp(StimStr);
disp(ParStr);
disp(RecStr);





