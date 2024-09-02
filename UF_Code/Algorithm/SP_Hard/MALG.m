function [results, y,S,W ] = SP( A, gnd,isrobust, islocal,reguType,maxIter)
% code for multi-view self-pace graph-based clustering
% A: cell contains the given similarity matrix from each view
% c: number of clusters
% isrobust = 1 means we use the 1-norm
% islocal = 1 means only nonzero values in A are updated
if ~exist('reguType','var');     reguType = 'exp';    end
if ~exist('isrobust','var');     isrobust = 1;        end
if ~exist('islocal','var');      islocal = 1;          end
if ~exist('maxIter','var');      maxIter = 50;         end
ncluster = length(unique(gnd));  viewNum = length(A); nSmp = size(A{1},1);

%%  %%%%%%%%%%%%% parameter
%  lambda1 = 8;   
lambda1 = randperm(30,1);
gamma =1; % gamma >= 1 controls the increase of lambda2 
zr = 10e-11;


%% %%%%%%%%%% initilize F 
A0 = zeros(nSmp);
for v = 1:viewNum
    A0 = A0 + A{v};
end
A0 = A0/viewNum;
A0 = A0-diag(diag(A0));
A10 = (A0+A0')/2;
D10 = diag(sum(A10));
L0 = D10 - A10;
[F0, ~, evs]=eig1(L0, nSmp, 0);
if sum(evs(1:ncluster+1)) < zr
    error('The original graph has more than %d connected component', ncluster);
end;
if sum(evs(1:ncluster)) < zr
    [clusternum, y]=graphconncomp(sparse(A10)); y = y';
    S = A0;
    return;
end;

F = F0(:,1:ncluster);

%% %%%% prepare for U (find the index of non-zero ones)
u = cell(viewNum,nSmp);
for i=1:nSmp
    for v = 1:viewNum
        u{v,i} = ones(1,nSmp);
    end
    if islocal == 1
        idxa0 = [];
        for v = 1:viewNum
            A0 = A{v};
            a0 = A0(i,:);
            idxa0 = union(idxa0,find(a0>0));
        end
        nullidx = setdiff(1:nSmp,idxa0);
        for v = 1:viewNum
            u{v,i}(nullidx) = 0;
        end
    end
end

%%     %%%%% initialize W
W = cell(viewNum,1);
for v = 1:viewNum
    W{v} = ones(nSmp,nSmp);
end


%% %%%%%%%%%%%%%%% begin iteration
for iter = 1:maxIter
     
    %%%%%%% update S
    dist = EuDist2(F,F,0);
    S = zeros(nSmp);
    for i=1:nSmp
        if islocal == 1
            idxa0 = find(u{1,i}>0);
        else
            idxa0 = 1:nSmp;
        end
        
        di = dist(i,idxa0);
        ua = zeros(1,length(idxa0));
        usum = ua;
        for v = 1:viewNum
            A0 = A{v};
            a0 = A0(i,:);
            if isrobust
                u{v,i} = u{v,i}.*abs(W{v}(i,:)); %%该步为什么这样计算？
            else
                u{v,i} = u{v,i}.*(W{v}(i,:).^2); 
            end
            ua = ua + (u{v,i}(idxa0)).*a0(idxa0); 
            usum = usum + u{v,i}(idxa0);
        end
        p = ua - lambda1*di/2;
        si = EProjSimplexdiag(p, usum);
        if isrobust
            for v = 1:viewNum
                ai = A{v}(i,idxa0);
                u{v,i}(idxa0) = 1./(2*sqrt((si-ai).^2+eps));
            end
        end
        S(i,idxa0) = si;
    end
    
    %%%%%%%% calculate loss
    Loss = cell(viewNum,1);
    for v = 1:viewNum
        if isrobust
            Loss{v} = abs(S - A{v});
        else
            Loss{v} = (S - A{v}).^2;
        end
    end
    lambda2 = updateLambda(Loss,gamma,iter); %% 在该文中，lambda设置基于样本的均值
    %%%%%%%% update W (choose sample pairs)
    W = updateW(Loss,reguType,lambda2);
    
    AA = S;
    AA = (AA+AA')/2;
    D = diag(sum(AA));
    L = D-AA;
    F_old = F;
    [F, ~, ev]=eig1(L, ncluster, 0);
    obj = 0;
    for v = 1:viewNum
        if isrobust
            obj = obj +sum(sum( abs(W{v}.*(S - A{v}))));
        else
            obj = obj + sum(sum((W{v}.*(S - A{v})).^2));
        end
    end
    obj = obj + lambda1*trace(F'*L*F);
    objs(iter) = obj;
    evs(:,iter+1) = ev;
    
    fn1 = sum(ev(1:ncluster));
    fn2 = sum(ev(1:ncluster+1));
    if fn1 > zr
        lambda1 = 2*lambda1;
    elseif fn2 < zr
        lambda1 = lambda1/2;  F = F_old;
    else
        break;
    end;
    
end;

[clusternum, y]=graphconncomp(sparse(AA)); y = y';
results = ClusteringMeasure(gnd,y);
if clusternum ~= ncluster
    sprintf('Can not find the correct cluster number: %d', ncluster)
end;


end





