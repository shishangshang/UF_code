function [res,a]= runpaceTest_update( Ls,Y,k,lambda,reguType)

[label,a] = solveMultiGraphpace_update(Ls, k, lambda,reguType);


% F=F./repmat(sqrt(sum(F.^2,2)),1,k);
% 
% % Final kmeans
% label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label); 

end