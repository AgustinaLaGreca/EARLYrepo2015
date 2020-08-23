function iRec = ExpandIsub(iRec, ds)
% ExpandIsub - convert Isub specification to vector containing the Isub

if isequal(0,iRec) 
   iRec=1:Ncond(ds); 
elseif isequal('up',iRec)
   [~, iRec] = sort(ds.Stim.Presentation.X.PlotVal(:).'); % row vector
elseif isequal('down',iRec)
   [~, iRec] = sort(ds.Stim.Presentation.X.PlotVal(:).'); % row vector
   iRec = iRec(end:-1:1);
elseif ~isnumeric(iRec), error('Invalid iSub specification');
end
