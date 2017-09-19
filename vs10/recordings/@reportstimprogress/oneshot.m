function oneshot(R, MessMode);
% reportstimprogress/oneshot - report progress one time
%    oneshot(R) performs the elementary action of R by reporting the
%    stimulus progress once. The status of R must be 'prepared' or
%    'started'.
%
%    The progress itself is queried by invoking the DAstatus on the
%    StimPresent information passed to R at constructor time.
%
%    oneshot(R, MessMode) uses a specific GUImessage mode for the display.
%    The default mode is 'shy'.
%
%    See also sortConditions, GUImessage, stimpresent/DAstatus.


if nargin<2, , MessMode = 'shy' ; end % default: don't bring GUI to the fore when displaying the info
R = download(R);
if ~isequal('prepared', status(R)) && ~isequal('started', status(R)),
    error('Status of R must be ''prepared'' or ''started''.');
end

St = DAstatus(R.SP);
ds_name = R.ds_name;
if isequal(-inf, St.icond), % not started yet; don't report anything
    Disp_Str = '';
elseif isequal(inf, St.icond), % ready
    Disp_Str = 'D/A finished';
else, % build string to be displayed on GUI 
    fillChar = '  ';
    if ~(R.Stutter)
        switch St.icond,
            case R.SP.Ncond+1, % pre-stim baseline
                Disp_Str = 'pre-stimulus baseline';
            case R.SP.Ncond+2, % post-stim baseline
                Disp_Str = 'post-stimulus baseline';
            otherwise, % regular stim conditions
                Disp_Str = ['pres ' MoutofNstr(St.ipres, R.SP.Npres, fillChar), '  --- ', ...
                    sprintf(R.SP.X.FormatString, R.SP.X.PlotVal(St.icond)) ];
                if ~isempty(R.SP.Y),
                    Disp_Str = [Disp_Str ' & ' sprintf(R.SP.Y.FormatString, R.SP.Y.PlotVal(St.icond))];
                end
                Disp_Str = [Disp_Str '  --- rep ' MoutofNstr(St.irep, R.SP.Nrep, fillChar)];
        end
    else % If the stutter is on
        switch St.icond,
            case R.SP.Ncond+1, % pre-stim baseline
                Disp_Str = 'pre-stimulus baseline';
            case R.SP.Ncond+2, % post-stim baseline
                Disp_Str = 'post-stimulus baseline';
            otherwise, % regular stim conditions
                Disp_Str = ['pres ' MoutofNstr(St.ipres, R.SP.Npres, fillChar), '  --- ', ...
                    sprintf(R.SP.X.FormatString, R.SP.X.PlotVal(St.icond)) ];
                if ~isempty(R.SP.Y),
                    Disp_Str = [Disp_Str ' & ' sprintf(R.SP.Y.FormatString, R.SP.Y.PlotVal(St.icond))];
                end
                Disp_Str = [Disp_Str '  --- rep ' MoutofNstr(St.irep, R.SP.Nrep, fillChar)];
    
        end
    end
    if strcmp(ds_name, '') 
        Disp_Str = {Disp_Str; ['Remaining D/A time: ' DurString(St.remtime)]}; % 2nd line: remaining time
    else
        Disp_Str = {ds_name; Disp_Str; ['Remaining D/A time: ' DurString(St.remtime)]}; % 2nd line: remaining time
    end
end
% display it
GUImessage(R.figh, Disp_Str, MessMode);






