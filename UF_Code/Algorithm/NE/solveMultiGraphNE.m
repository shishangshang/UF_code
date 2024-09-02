function F = solveMultiGraphNE(Ls, k, lambda)

NITER = 20;
c = length(Ls);
[~, n] = size(Ls{1});

a = ones(1, c) ./ c;         

for iter = 1:NITER                        
    % solve eign system
    S = zeros(size(Ls{1}));
    for v = 1:c
        S = S + a(v).*Ls{v};
    end
    [F, val] = eig1(S, k, 0);
    
    obj_v = zeros(1, c);
    for v = 1:c
        obj_v(v) = trace(F'*Ls{v}*F);
    end
    
%     % solve alphas via cvx
%     cvx_begin quiet
% %        cvx_precision low
%         variable a(c)
%         minimize( obj_v*a + lambda*sum(+log(a).*a))
%         subject to
%         a >= 0
%         sum(a) == 1
%     cvx_end
     %% solve alpha via exp
      a=exp(-obj_v/lambda);
      a = a ./ sum(a);
     %% solve alpha via equal weight
%      a = ones(1, c) ./ c;  
end

end