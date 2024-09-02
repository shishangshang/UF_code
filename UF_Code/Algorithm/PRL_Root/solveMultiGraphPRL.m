function [F,a,fs_star] = solveMultiGraphPRL(Ls, k,manner,sata)

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
    switch manner
      case 'root'
    % solve alphas via root
    a = obj_v.^(p-1);
      case 'log'
    % solve alphas via log
    a= 1./obj_v.*log(obj_v);
      case 'capped'  
    % solve alphas via capped
      for v=1:c
          if obj_v(v)<sata
              a(v)=1;
          else
              a(v)=0;
          end
      end
    end
    
end

fs_star = obj_v;

end