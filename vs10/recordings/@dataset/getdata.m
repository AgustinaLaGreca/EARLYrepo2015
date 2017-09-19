function D=getdata(D);
% Dataset/getdata - get data from the suppliers of dataset object
%    CD = getdata(DS) visits the suppliers S, .. of DS one by one, and invokes
%    getdata(S1), .. to collect their data. The the return value CD is a 
%    copy of DS with the data in place. DS itself is not affected 
%    (no "uploading" of DS). This copy can be used for realtime data
%    analysis during data collection. 
%
%    getdata(DS) with no output args *does* upload DS with the data
%    collected. This call changes DS irreversibly, so it can only be done
%    once - typically at wrapup time. A typical construction is to create a
%    "dataset finalizer" object as is done in PlayRecordStop:
%
%         actabit(figh, 'finalizeDS', 'wrapup', {@getdata DS}, 'clear', {@save DS});
%
%    This finalizer should be placed at the end of the action list of the
%    GUI, after the suppliers of DS. This way the suppliers are wrapped up
%    first, so dataset/getdata gets the wrapped-up versions of the data.
%
%   See Dataset/addsupplier, grabevents/getdata, PlayRecordStop.

if ~isbeingcollected(D), % getdata has already done its work
    return;
end
D = download(D);
Dtypes = fieldnames(D.Data); SupStatus = {}; StopSamples = [];
for ii=1:numel(Dtypes), % visit the D.Data fields & get them data
    dtype = Dtypes{ii};
    datafield = D.Data.(dtype); % current value; may either hold Grabbeddata or Datagrabber
    if isa(datafield, 'address'), datafield = get(datafield); end % pointer -> content
    if isa(datafield,'datagrabber'), % replace it by its current data
        D.Data.(dtype) = getdata(datafield);
    elseif isa(D.Data.(dtype),'grabbeddata'), % already filled w data - get fresh data from supplier
        SuSt = supplierstatus(datafield);
        % is there anything new? Then get it
        if ~isequal('finished',SuSt) && ~isequal('stopped',SuSt),
            D.Data.(dtype) = getdata(supplier(datafield));
        end
    else, 
        datafield, itsclass = class(datafield),
        error(['''' dtype ''' field of D.Data is neither a datagrabber nor a grabbeddata object.'])
    end
    StopSamples = [StopSamples, stopstatus(datafield)];
    SupStatus{ii} = supplierstatus(D.Data.(dtype));
end
if isequal({'finished'}, unique(SupStatus)), D.Status = 'complete'; % all suppliers have finished
elseif ismember('started', SupStatus), D.Status = 'expecting_more'; % at least one supplier is still collecting data
elseif ismember('stopped', SupStatus), D.Status = 'interrupted'; % at least one supplier was interrupted, none is running anymore
else, D.Status = 'premature'; % not all suppliers have started
end
% store progress of D/A when stopped - if stopped
istop = min(StopSamples);
if ~isempty(istop),
    D.Stopped = DAstatus(D.Stim.Presentation, istop);
end

D.ID.modified = datestr(now);

if nargout<1,
    upload(D);
end


