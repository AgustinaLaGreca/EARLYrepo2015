function D = getdatabuf(G, varargin);
% datagrabber/getdata - get data buffer from datagrabber object
%   getdatabuf(G) returns the databuffer filled by datagrabber G. 
%   This fcn is typically use by the getdata methods of datagrabber
%   subclasses.
%
%   getdatabuf(G, 'Foo') only returns field Foo of the databuffer.
%
%   getdatabuf(G, 'Foo', 'Goo') only returns subfield Foo.Goo, etc.
%
%   See also datagrabber/setdatabuf, Grabevents/getdata, Dataset/getdata.

D = G.DataBuf;
for ii=1:nargin-1,
    D = D.(varargin{ii});
end









