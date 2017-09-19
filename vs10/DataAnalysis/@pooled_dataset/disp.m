function disp(D);
% Pooled_dataset/disp - DISP for Pooled_dataset objects.
%
%   See Pooled_dataset.


% ======single element from here=======
% if ~D.isstatic,
%     try, D = download(D); end;
% end
%=====non-empty, non-void D from here=======
for ii=1:numel(D),
    disp('--Pooled dataset--')
    disp(D(ii).DS);
    disp('------------------')
end




