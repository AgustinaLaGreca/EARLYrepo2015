function AllRueCorrsNorm(Postfix, varargin);
% AllRueCorrsNorm - crosscorr analysis, Normalized, among Rue's three test sets
%  

if nargin<1, Postfix = ''; end

Ddir = 'D:\Data\RueData';

FileName = fullfile(Ddir,'RueCorNorm', Postfix);

RueCorNorm.AA = RueCorrNorm('A', 'A', varargin{:});
RueCorNorm.BB = RueCorrNorm('B', 'B', varargin{:});
RueCorNorm.ZZ = RueCorrNorm('Z', 'Z', varargin{:});
RueCorNorm.AB = RueCorrNorm('A', 'B', varargin{:});
RueCorNorm.BZ = RueCorrNorm('B', 'Z', varargin{:});
RueCorNorm.ZA = RueCorrNorm('Z', 'A', varargin{:});

save(FileName, 'RueCorNorm');















