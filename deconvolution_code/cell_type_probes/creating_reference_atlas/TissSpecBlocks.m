%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file was adapted from Joshua Moss and Tommy Kaplan from the paper                                  %
% Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease %
% doi: https://doi.org/10.1038/s41467-018-07466-6                                                                    %
% Any reuse of this code should include that attribution                                                             %               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% For Research Use only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,tabu,indic,indic2] = TissSpecBlocks(X,n,chrs,chr,pos,rng)
[N,d]=size(X);
X1 = bsxfun(@rdivide, X, nansum(X,2));
X2 = 1-X; X2 = bsxfun(@rdivide, X2, nansum(X2,2));

S=[]; indic =[];
for i=1:d,
    [~,J]=sort(X1(:,i),'descend');
    S = [S; J(1:n)];
    indic = [indic; repmat(i,n,1)];
end

for i=1:d,
    [~,J]=sort(X2(:,i),'descend');
    S = [S; J(1:n)];
    indic = [indic; repmat(-i,n,1)];
end

% include neigh
seenchr = unique(chr(S));
for i=1:length(seenchr),
    c = seenchr(i);
    I2{c} = find(chr==c); % all CpGs in chr
end

tabu = []; indic2 = [];
for i=1:length(S),
    si = S(i);
    c = chr(si);
    I = I2{c}(abs(pos(I2{c}) - pos(si))<=rng);
    tabu = [tabu;I];
    indic2 = [indic2; repmat(indic(i),length(I),1)];
end
