function [F,lambda2] = solveMultiGraphpace(Ls, k,gamma,reguType)
if ~exist('reguType','var');     reguType = 'exp';    end
NITER = 20;
c = length(Ls);
%[~, n] = size(Ls{1});

a=cell(1,c);
for v=1:c
    a{v}=1/c;
end

for iter = 1:NITER
    % solve eign system
    S = zeros(size(Ls{1}));
    for v = 1:c
        S = S + a{v}.*Ls{v};
    end
    [F, ~] = eig1(S, k, 0);
    
    obj_v = cell(1, c);
    for v = 1:c
        obj_v{v} = trace(F'*Ls{v}*F);
    end
    
    % solve alphas via self-paced
     lambda2 = updateLambda(obj_v,gamma,iter);
         %%%%%%%% update W (choose sample pairs)
     a = updatea(obj_v,reguType,lambda2);
     
%    a = obj_v.^(1./(1-gamma));
%    a = a ./ sum(a);
end

fs_star = obj_v;

end