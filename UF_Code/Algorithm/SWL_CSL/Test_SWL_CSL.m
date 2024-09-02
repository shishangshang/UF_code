function [res,label, fs_star] = Test_SWL_CSL(Ls, k, Y,manner,sata)

[F, fs_star] = solveMultiGraphPRL(Ls, k,manner,sata);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end