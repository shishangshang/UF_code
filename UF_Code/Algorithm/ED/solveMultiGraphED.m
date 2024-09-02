function [F, fs_star] = solveMultiGraphED(Ls, k, gamma)

NITER = 20;
c = length(Ls);
%[~, n] = size(Ls{1});

a = ones(1, c) ./ c;

for iter = 1:NITER
    % solve eign system
    S = zeros(size(Ls{1}));
    for v = 1:c
        S = S + a(v).*Ls{v};
    end
    [F, ~] = eig1(S, k, 0);
    
    obj_v = zeros(1, c);
    for v = 1:c
        obj_v(v) = trace(F'*Ls{v}*F);
    end
    
    % solve alphas via ED
%   a = obj_v.^(gamma./(1-gamma));
    a = obj_v.^(1./(1-gamma)); %%correcting by shi
    a = a ./ sum(a);
end

fs_star = obj_v;

end