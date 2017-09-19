function [S, T] = ZWOAEmatrices(Nzw, DPtype);
% ZWOAEmatrices - matrices for the evaluation of 3rd order DPs
%   ZWOAEmatrices(Nzw) returns a struct whose fields are matrices for 
%   computing frequencies and phases of 3rd-order DP evoked  
%   by presenting a group of Nzw primaries ("F1 group") together with a 
%   single F2 primary tone (note that F2 may be either below OR above the
%   F1's). 
%
%   Specifically, let F be a (Nzw+1)x1 column vector containing all
%   primary frequencies: Fprim == [F1(1); F1(2); ..., F1(Nzw); F2].
%   The elements of F1 must be in ascending order. Then the DP frequencies 
%   of type Foo are given by
%
%           Fdp = Mfoo*Fprim
%
%   The different types of DPs are indicated by the fieldname of the
%   returned structure:
%       Mnear: DPs near F1 of type F1(i)+F1(j)-F2
%        Mfar: DPs near F2 of type 2*F2-F1(i)
%      Msuplo: DPs just below F2 ("suppression") of type F2+F(i)-F(j), i<j
%      Msuphi: DPs just above F2 ("suppression") of type F2-F(i)+F(j), i<j
%     Msupall: the 2 previous type together, i.e., Msupall = [Msuplo; Msuphi]
%        Mall: all previous types, i.e. Mall = [Mnear; Mfar; Msuplo; Msuphi]
%
%   [S, T] = ZWOAEmatrices(Nzw) also returns a struct T whose fieldnames
%   are identical to those of S. Each field Foo is a Nx3 matrix MI whose
%   jth row MI(j,:) == [k l m] are the indices of the interacting primaries
%   that give rise to the DP described by the jth row of the Mfoo matrix
%   described above. Specifically, the DP of the jth row is evoked by the
%   interaction of Fprim([k l m]), with Fdp(j) = Fprim(k)+Fprim(l)-Fprim(m).
%   Note that k & l occur symmetrically in this formula. Per convention, we
%   choose j<=k. Also note that j and k may coincide, in which case we have
%   the "classical" case 2*Fprim(k)-Fprim(m).
%
%   [S, T] = ZWOAEmatrices(Nzw, Foo) only returns the matrices of type Foo, 
%   where Foo is one of the fieldnames in the above list (abbreviations 
%   and case mismatches are allowed, as long as there is no ambiguity).
%
%   Note that the Fdp computed by Mfoo*Fprim are generally not sorted, so
%   care must be taken when unwrapping phases, etc.
%
%   See also ZWOAEplotStruct.

if nargin<2, 
    DPtype = ''; % return whole struct
end

% aux matrices
Md = -EvaDifMat(Nzw); % differences across Nzw components; if F1 are excending, Md*F1 are positive
Ms = [abs(Md); 2*eye(Nzw)];  % sums across Nzw components, including same-term sums

% DP matrices themselves
Mnear = [Ms, -ones(size(Ms,1),1)]; % F1(i)+F1(j)-F2 cmp
Mfar =  [-eye(Nzw), 2*ones(Nzw,1)]; % 2*F2-F1 cmp
Msuplo = [-Md, ones(size(Md,1),1)]; % F2 +F1b-F1a components (upper suppr band)
Msuphi = [Md, ones(size(Md,1),1)]; % F2 -F1b+F1a components (lower suppr band)
Msupall = [Msuplo; Msuphi];
Mall = [Mnear; Mfar; Msupall];
S = CollectInStruct(Mnear, Mfar, Msuplo, Msuphi, Msupall, Mall);

% re-analyze the rows of the matrices Mfoo to find the 3 interacting
% primaries
FNS = fieldnames(S);
for ii=1:length(FNS),
    fn = FNS{ii};
    M = S.(fn);
    Nrow = size(M,1);
    F = zeros(Nrow,3);
    for irow=1:Nrow,
        row = M(irow,:);
        itwo = find(row==2, 1);
        iones = find(row==1);
        iminus = find(row==-1);
        if ~isempty(itwo),
            F(irow,1:2) = itwo;
        else,
            F(irow,1:2) = iones;
        end
        F(irow,3) = iminus;
    end
    T.(fn) = F;
end

if ~isempty(DPtype), % select specified field
    [Fieldname, Mess] = keywordMatch(DPtype, fieldnames(S), 'distortion product type');
    error(Mess);
    S = S.(Fieldname);
    T = T.(Fieldname);
end





