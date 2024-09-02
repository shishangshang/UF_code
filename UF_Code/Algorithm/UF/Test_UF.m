function [res,label] = runParameterizedTest(Ls, k, Y, lambda, q)

F = solveMultiGraphUP(Ls, k, lambda, q);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end