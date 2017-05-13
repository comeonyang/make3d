function [F_struc,F_blksz,G_struc,G_blksz]  = sedumi2maxdet(F_struc,K)
%SEDUMI2MAXDET Internal function to convert SeDuMi structure to format needed in MAXDET

% Author Johan Löfberg
% $Id: sedumi2maxdet.m,v 1.2 2004/07/02 08:17:32 johanl Exp $

switch K.m
    case 0
        G_struc = [];
        G_blksz = [];        
        F_struc = F_struc;
        F_blksz = [repmat(1,1,K.l) K.s];
    case 1
        G_struc = F_struc(K.l,:);
        G_blksz = [1];        
        F_blksz = [repmat(1,1,K.l-1) K.s];
        F_struc = [F_struc(1:1:K.l-1,:);F_struc(K.l+1:1:end,:)];
    otherwise
        G_struc = F_struc(end-K.s(end)^2+1:end,:);
        G_blksz = K.s(end);        
        F_blksz = [repmat(1,1,K.l) K.s(1:end-1)];
        F_struc = F_struc(1:end-K.s(end)^2,:);
end
     