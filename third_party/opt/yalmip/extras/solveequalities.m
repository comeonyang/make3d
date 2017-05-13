function [x_equ,H,A_equ,b_equ,factors] = solveequalities(F_struc,K,unitary)
%SOLVEEQUALITIES Internal function remove equality constraints

% Author Johan Löfberg
% $Id: solveequalities.m,v 1.12 2006/05/14 17:21:09 joloef Exp $

% Extract the inequalities
A_equ = F_struc(1:K.f,2:end);
b_equ = -F_struc(1:K.f,1);

if nargin<3
    unitary = 1;
end

factors = [];

% Improve numerics by removing obviously
% redundant constraints
hash = 1+rand(1+size(A_equ,2),1);
[i,j] = unique([A_equ b_equ]*hash);
remove =  setdiff(1:size(A_equ,1),j);
if ~isempty(remove)
    A_equ(remove,:) = [];
    b_equ(remove,:) = [];
end

if ~unitary
    % Just use a crappy basis derived from A
    % FIX : fix dependent row case
    if 0
        [L,U,P] = lu(A_equ');
        s = find(diag(U));
        r = setdiff(1:length(U),s);
        n = size(s);
        [i,j] = find(P');
        H1 = A_equ(:,i(s));
        H2 = A_equ(:,i(r));
        x_equ = P'*(L'\(U'\b_equ));
        % FIX : use L and U stupid!
        H = P'*[-H1\H2;eye(size(H2,2))];

    else
        [L,U,P] = lu(A_equ');
        n = max(find(diag(U)));
        [i,j] = find(P');
        H1 = A_equ(:,i(1:n));
        H2 = A_equ(:,i(n+1:end));
        try
            x_equ = P'*linsolve(full(L'),linsolve(full(U'),full(b_equ),struct('LT',1==1)),struct('UT',1==1));
        catch
            x_equ = A_equ\b_equ;
        end
        % FIX : use L and U stupid!
        H = P'*[-H1\H2;eye(size(H2,2))];
    end
else
    % Use unitary basis
    try
        [Q,R,E] = qr(full(A_equ)');
    catch
        [Q,R,E] = qr(A_equ'); % Ouch, that big!
    end
    n = max(find(sum(abs(R),2)>1e-14*size(R,2)));

    Q1 = Q(:,1:n);
    R = R(1:n,:);
    x_equ = Q1*(R'\E'*b_equ);
    H = Q(:,n+1:end); % New basis  
end
