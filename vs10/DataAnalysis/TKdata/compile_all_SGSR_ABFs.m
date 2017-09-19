function S = compile_all_SGSR_ABFs(Experimenter);
% compile_all_SGSR_ABFs - combine all SGSR_ABF files and save.
%   compile_all_SGSR_ABFs(Experimenter) looks for all experiments in
%   TKdatadir(Experimenter) and compiles a database linking the ABF files
%   to their SGSR counterparts by calling SGSR_ABF_link.
% 
%  See also SGSR_ABF_link, ABFupdata, SGSRupdata.

if nargin<1, Experimenter='TK'; end % Thomas the pioneer.

Ddir = TKdatadir(Experimenter);
% find out which experiments are in Ddir.
qq = dir(Ddir);
qq = qq([qq.isdir]);
qq = qq(3:end); % remove .\ and ..\

S = [];
for ii=1:numel(qq),
    xp = qq(ii).name;
    if ~isequal('SGSR_ABF_LINK', upper(xp)),
        try,
            [s, Mess] = SGSR_ABF_link(xp);
            warning(Mess);
        catch,
            warning(['SGSR_ABF_link crashed for ''' xp '''.'])
        end
        if ~isempty(s),
            if isempty(S), S=s; else, S=[S,s]; end
        end
    end
end

% combine all .SGSR_ABF files in pool
PoolFile = fullfile(TKdatadir(Experimenter), 'SGSR_ABF_link' , 'All_Exps.SGSR_ABF');
save(PoolFile, 'S');













