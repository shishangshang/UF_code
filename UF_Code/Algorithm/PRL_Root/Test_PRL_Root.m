function [res,label, fs_star] = runPRLTest(Ls, k, Y, p,manner)

[F, fs_star] = solveMultiGraphPRL(Ls, k, p,manner);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end