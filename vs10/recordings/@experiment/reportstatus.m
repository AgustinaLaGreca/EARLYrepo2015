function S = reportstatus(E, M);
% experiment/reportstatus - report status of experiment
%    reportstatus(E) returns a cell string describing the current status 
%    of experiment E. This is used by dashboard.
%    
%    reportstatus(E,M) displays this info to messenger M.
%
%    See also experiment/status, dashboard, messenger/report.

if nargin<2, M=[]; end

if isvoid(E),
    S = 'No experiment defined.';
else,
    St = status(E);
    S = {[upper(name(E)) ' (' state(E) '); ' num2str(St.Ndataset) ' saved.'], ...
        ['Last saved DS = <' St.IDlastSaved '>'], ...
        ['(Unit,k) = (' num2str(St.iCell) ', ' num2str(St.iRecOfCell) ')'], ...
        ['Electrode # ' num2str(St.iPen) ';  depth = ' num2str(St.PenDepth) ' um'] ...
        };
end
if ~isempty(M),
    report(M,S);
end










