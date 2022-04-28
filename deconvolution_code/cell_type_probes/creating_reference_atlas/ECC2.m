%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file was adapted from Joshua Moss and Tommy Kaplan from the paper                                  %
% Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease %
% doi: https://doi.org/10.1038/s41467-018-07466-6                                                                    %
% Any reuse of this code should include that attribution                                                             %               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% For Research Use only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,L,Sij] = ECC(X,NN,tabu,sat)
% argmax_S argmin_{i,j} ||X_i-X_j||^2

[N,d]=size(X);
X1 = bsxfun(@rdivide, X, nansum(X,2));
X2 = 1-X; X2 = bsxfun(@rdivide, X2, nansum(X2,2));

% find all pairs of tissues
alpha = triu(ones(d,d),1);

% compute squared-differences between Xi and Xj
[I,J]=find(alpha); IJ=sub2ind([d,d],I,J); Lij = length(IJ);
Xij = ones(N,Lij,'single');
for i=1:length(IJ)
    Xij(:,i) = (X(:,I(i))-X(:,J(i))).^2;
end
clear I J

% initialize S: set of selected features
S=ones(0,1); Sij=S;

% I=setdiff(find(var(X,[],2)>0.01),tabu);
% S1=CompSensing(X(I,:),[]); S=[S;I(S1)];
% S2=CompSensing(X(I,:),S1); S=[S;I(S2)];
L=NaN*zeros(1,NN);
% df = var(X,[],2);

for iter=1:NN,
    % find two nearest cell-types i and j, given current selection S
    % pairwise distance upon S
    dis = sqrt(sum(Xij(S,:),1));
    saturated = find(histc(Sij(:,1),1:length(IJ))>=sat);
    dis(saturated)=NaN;
    [mindis,ij]=min(dis);
    if isnan(mindis), break; end;
    L(iter)=mindis;
    [i,j]=ind2sub([d d],IJ(ij));
    % find CpG maximizing the distance between two types
    df = Xij(:,ij);
    % ignore tabu or already-selected CpGs
    df([S;tabu])=0;
    % add k into set S
    [~,k]=max(df);
    S=[S;k];

    % record why k was selected
    Sij=[Sij;[ij i j]]; % ij index in Xij, convert to IJ(ij) to get pair-index
    % fprintf('k:%d ij %d i:j %d:%d %.3f\n',k,ij,i,j,sqrt(df(k)));

    if 0,
	I1 = find(mean(X(S,:),2)<0.5); [~,II]=max(X1(S(I1),:),[],2); [~,J1]=sort(II); J1=I1(J1);
	I2 = find(mean(X(S,:),2)>=0.5);[~,II]=min(X1(S(I2),:),[],2); [~,J2]=sort(II); J2=I2(J2);
	S=S([J1;J2]);
	clear I1 I2 I J1 J2

	if mod(iter,10)==0,
	    subplot(4,1,1:3);imagesc(X(S,:),[0 1]);
	    subplot(4,1,4); bar(1:NN,L); axis([0 NN 0 max([5,L])]); drawnow;
	end
    end
end

if 0,
    X = importdata('Atlas.csv',',',1);
    hdrs = X.textdata(1,2:end);
    cg = X.textdata(2:end,1);
    X = X.data;

    I=find(var(X,[],2)>0.01 & ~any(isnan(X),2));
    S=ECC(X(I,:),1000,[]);
    S = I(S);
end
