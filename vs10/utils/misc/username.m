function Name = username(Name,varargin);
% username - name of current user
%   Username returns a character string identifying the user account.
%
%   Username('John') sets the username to 'John'. Usernames must be valid
%   Matlab identifiers (see ISVARNAME). Once specified, the username is
%   remembered across Matlab sessions. Different accounts on the same PC
%   have their own username.
%
%   vargargin{1} can be to ask now for a username, usefull for the startup
%   function.
%
%   Note: these names are meaningful only within EARLY.
%
%   See also Compuname, Where.


FileName = fullfile(prefdir(1), 'userID.EARLY'); % setup file must be kept outside EARLY dirs

    
if nargin<1, % get
    Name = MyFlag('Username___'); % try to retrieve from persistent flag (no file I/O)
    if isempty(Name),
        try,
            qq = load(FileName, '-mat'); Name = qq.Name;
            MyFlag('Username___', Name); % store as persistent flag
        catch, % not known - set it now
            while 1,
                UN = input('Specify username: ', 's');
                if isvarname(UN), break; end
                disp('Username must be valid Matlab identifier.')
            end
            username(UN);
            Name = UN;
        end
    end
else, % set
    if isnumeric(Name)
        while 1,
            UN = input('Specify username: ', 's');
            if isvarname(UN), break; end
            disp('Username must be valid Matlab identifier.')
        end
        username(UN);
        Name = UN;
        return;
    end
    if ~isvarname(Name),
        error('Username must be valid Matlab identifier.');
    end
    save(lower(FileName), 'Name', '-mat');
    MyFlag('Username___', Name); % store as persistent flag
end

