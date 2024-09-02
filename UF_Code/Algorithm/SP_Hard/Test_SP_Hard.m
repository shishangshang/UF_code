function [res,lambda2]= runpaceTest( Ls,k, Y,lambda,reguType)

[F,lambda2] = solveMultiGraphpace(Ls, k, lambda,reguType);


F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end