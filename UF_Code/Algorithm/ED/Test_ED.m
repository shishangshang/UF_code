function [res,label,fs_star] = Test_ED(Ls, k, Y, gamma)

[F, fs_star] = solveMultiGraphED(Ls, k, gamma);

F=F./repmat(sqrt(sum(F.^2,2)),1,k);

% Final kmeans
label=litekmeans(F,k, 'Replicates', 10);
res = ClusteringMeasure(Y, label);

end