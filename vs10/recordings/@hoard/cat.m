function cat(H,X, Dim);
% hoard/cat - concatenate data to hoard object
%    cat(H,X) concatenates data X to hoard X.
%
%    cat(H,X,Dim) concatenates along dimension Dim. Default Dim=1;
%
%    See also hoard/set.

%qq=dbstack; {qq.name}
%H.Name
if nargin<3, Dim=1; end
[Xold, Nch] = get(H);
Xnew = cat(Dim, Xold, X);
% dsize(Xold, X, Xnew)
if numel(Xnew)>=Nch,
    [Chunk, Xnew] = deal(Xnew(1:Nch), Xnew(Nch+1:end));
else, Chunk=[];
end

if isempty(Chunk), % simply replace H's data
    set(H, Xnew);
else, % save Xnew to H, Chunk to file
    % update userdata fields & save them quickly (before the next call to cat)
    UD = get(H.h, 'userdata');
    UD.data = Xnew;
    UD.NsamSaved = [UD.NsamSaved numel(Chunk)];
    ichunk = UD.NchunkSaved+1;
    Xmax = max(abs(Chunk));
    if Xmax==0, Xmax=1; end; % avoid dividing by zero
    UD.ScaleFactor = [UD.ScaleFactor, Xmax/(2^31-1)];
    UD.NchunkSaved = ichunk;
    set(H.h, 'userdata', UD);
    % save samples
    Chunk = int32(Chunk/UD.ScaleFactor(end));
    Filename = [tempfileprefix(H) num2str(ichunk) '.bin'];
    fid = fopen(Filename, 'wb');
    fwrite(fid, Chunk, 'int32');
    fclose(fid);
end







