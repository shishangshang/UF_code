function [res,label] = runNETest(Ls, k, Y, lambda)

F = solveMultiGraphNE(Ls, k, lambda);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end