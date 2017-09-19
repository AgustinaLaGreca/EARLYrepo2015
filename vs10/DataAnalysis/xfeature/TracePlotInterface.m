function TracePlotInterface(figh,link);
% TracePlotInterface - enable keyboard navigation for plots of recording traces
%   TracePlotInterface(Figh) enables navigation keys for plots of recording
%   traces. Figh is the figure handle.
%   Here is a list of the fast keys and their meaning
%
%            leftarrow: backward in time, small step
%      shift-leftarrow: backward in time, large step
%           rightarrow: forward in time, small step
%     shift-rightarrow: forward in time, large step
%             spacebar: forward in time, large step
%
%              uparrow: slow zoom in
%        shift uparrow: fast zoom in
%            downarrow: slow zoom out
%      shift downarrow: fast zoom out
%
%         home: jump to start of trace
%          end: jump to end of trace
%            0: jump to start of trace, set X range to [0 30] ms
%            A: auto X-range (complete trace)
%            s: jump to next stimulus onset
%            S: jump to previous stimulus onset
%
%       delete: close figure
%       shift-delete: close all figures
%            d: display stimulus properties to command window
%            h: help TracePlotInterface (display this text)
%            ?: same as h
%            k: set focus to command window
%            P: impose position/size convenient for time analysis

if nargin < 2, link = true; end
ha = findobj(figh, 'type', 'axes', '-not', 'tag', 'legend');
for ii=1:numel(ha), 
    ylim(ha(ii), ylim(ha(ii))); 
end; % fix ylimits
set(figh, 'keypressfcn', @local_CallBack);
if ~isempty(ha) && link, linkaxes(ha,'x'); end

function local_CallBack(Src,Ev,varargin);
Shifted = ismember('shift', Ev.Modifier);
%Ev.Key, Ev.Modifier
shh = get(0,'showhidden');
set(0,'showhidden','on'); % axes to be affected by be hidden in a GUI; temporarlily force their visibility
switch Ev.Key,
    case 'leftarrow', % move to left
        if Shifted, Factor = 1; else, Factor = 0.1; end
        xlim(xlim-0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim-0.5*Factor*diff(xlim)); 
    case 'rightarrow', % move to right
        if Shifted, Factor = 1; else, Factor = 0.1; end
        xlim(xlim+0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim+0.5*Factor*diff(xlim)); 
    case 'space', % move to right
        Factor = 1; 
        xlim(xlim+0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim+0.5*Factor*diff(xlim)); 
    case 'backspace', % move to right
        Factor = 1; 
        xlim(xlim-0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim-0.5*Factor*diff(xlim)); 
    case 'uparrow',  % zoom in
        if Shifted, Factor = 0.5; else, Factor = 0.9; end
        xlim(mean(xlim)+0.5*sqrt(Factor)*[-1 1]*diff(xlim)); pause(0.025);
        xlim(mean(xlim)+0.5*sqrt(Factor)*[-1 1]*diff(xlim));
    case 'downarrow', % zoom out
        if Shifted, Factor = 0.5; else, Factor = 0.9; end
        xlim(mean(xlim)+0.5/sqrt(Factor)*[-1 1]*diff(xlim)); pause(0.025);
        xlim(mean(xlim)+0.5/sqrt(Factor)*[-1 1]*diff(xlim));
    case 'home', % to start of trace
        xlim(xlim-min(xlim));
    case 'end', % to end of trace
        DX = diff(xlim);
        xlim('auto');
        xlim(max(xlim)+[-DX 0]);
    case 'a', % default view
        if Shifted, xlim('auto'); end
    case 'delete', % close figure
        delete(gcbf);
        if Shifted, aa; end
    case '0',
        xlim([0 30]);
    case {'k' 'escape'},
        commandwindow;
    case {'h', 'slash'},
        commandwindow;
        help(mfilename);
    case 'g',
        TimeTo = str2double(inputdlg('Specify X value: '));
        if ~isempty(TimeTo),
            xlim(xlim-mean(xlim)+TimeTo);
        end
    case 's',
        M=getGUIdata(gcbf,'StimProps',[]); % retrieve stimulus info
        if ~isempty(M),
            XL = xlim; DX=diff(XL); X0=XL(1);
            if ~Shifted, % jump to next onset
                X0 = min(M.onset(M.onset>=X0+1.1));
            else, % previous
                X0 = max(M.onset(M.onset<=X0-1.1));
            end
            if ~isempty(X0),
                xlim((X0-1)+DX*[0 1]);
            end
        end
    case 'd',
        M=getGUIdata(gcbf,'StimProps'); % retrieve stimulus info
        if ~isempty(M),
            disp(M);
        end
    case 'p', % default view
        if Shifted, set(gcbf,'units', 'normalized', 'position', [0.00625 0.534 0.987 0.376]); end
    otherwise, 
        %disp(Ev)
end % swicth/case
set(0,'showhidden',shh);
