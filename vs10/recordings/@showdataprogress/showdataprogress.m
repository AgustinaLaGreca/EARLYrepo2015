function S = showdataprogress(adDS, DT, dataviewer, figh, Name, You);
% showdataprogress - data analysis plot during data collection
%    S=showdataprogress(DS, DT, Viewer, figh, 'XYZ', 'You') constructs an 
%    Action object S whose task it is to present online data analysis 
%    during data collection. S is added to the list of Action of the figure
%    with handle figh under the name XYZ, and it will be finished as soon
%    as the Action object named You is finished. The dataset DS to be 
%    analyzed is assumed to have an address (in fact, it is okay to pass 
%    the address instead of the dataset). Every DT ms, an updated copy of 
%    DS is passed to function Dataviewer, which performs the desired data 
%    analysis. At wrapup time, the data are once more updated and passed 
%    to the Dataviewer. Dataviewer is typically a function handle; if it is
%    a char string, it will be converted to a function handle.
%
%    A special value for the Dataviewer is '-', in which case
%    showdataprogress calls Nope (=do nothing).
%
%    S=showdataprogress(DS, DT, {@Foo Arg1 ..}) will pass additional 
%    arguments Arg1 .. to Foo, in addition the its obligatory first 
%    argument DS.
%
%    See also Action, dataset/getdata.


if ~isa(adDS, 'address'), adDS = address(adDS); end

if isequal('-', dataviewer), % no dataviewer specified; launch none
    return;
elseif ischar(dataviewer),
    dataviewer = fhandle(dataviewer);
elseif iscell(dataviewer) && ischar(dataviewer{1})
    dataviewer{1} = fhandle(dataviewer{1});
end

dataviewer = cellify(dataviewer); % easier passing to feval
Tlast = clock;
S = CollectInStruct(adDS, Tlast, DT, dataviewer); % note that DT is not passed to Action constructor. Instead, the timer ..
S = class(S,mfilename, action('initialized',[200 2000])); % ... interval is a fixed 200 ms; see ....
S = upload(S, figh, Name); % ...showdataprogress/oneshot why. The long start delay gives competing actions a chance 
afteryou(S, figh, You);






