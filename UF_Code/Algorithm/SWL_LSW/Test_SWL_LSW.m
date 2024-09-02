function [res,label,a, fs_star] = Test_SWL_LSW(Ls, k, Y, manner)

[F,a,fs_star] = solveMultiGraphPRL(Ls, k,manner);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end