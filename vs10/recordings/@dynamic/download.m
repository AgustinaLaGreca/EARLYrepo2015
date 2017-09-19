function D=download(D,flag);
% Dynamic/download - download updated version of Dynamic object
%   D=download(D) retrieves the value of D that was previously stored
%   using upload(D). This is a hack to use pointer-like programming
%   constructs. If D has not previously been uploaded, attempting to
%   download generates an error.
%
%   D=download(Ax) where Ax is an address as returned by GUIaccess does the
%   same (it's handled by address/download). It is equivalent to D=Ax.get().
%
%   D=download(D, 'sloppy') suppresses the error of non-uploaded S and
%   returns D unchanged.
%
%   D=download(D, 'homeless') downloads D and returns D with teh address
%   removed (leaving the uploaded version unaffected).
%
%   See Dynamic/upload, dynamic/hasbeenuploaded.

if nargin<2, flag=''; end

if isempty(D.access),
    if ~isequal('sloppy',flag),
        error('Cannot download dynamic object that has never been uploaded. Use the ''sloppy'' flag if you don''t care.');
    end
else,
    D = get(D.access); 
    if isequal('homeless', flag), % remove address
        D.access = [];
    end
end


