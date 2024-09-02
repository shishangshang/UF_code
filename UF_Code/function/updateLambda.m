function lambda = updateLambda(Loss,gamma,iter)
viewNum = length(Loss);
lambda = zeros(viewNum,1);
for v = 1:viewNum
    idx = Loss{v} > 0;
    temp = median(Loss{v}(idx));
    lambda(v) = temp + log(iter*(temp^2 + gamma));
end
end