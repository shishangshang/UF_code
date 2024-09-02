function [W] = updateW(Loss, reguType,lambda)
% % Loss: cell, contains the simlarity differences of each view
% reguType: the regularization type imposed on W

% lambda: parameter
viewNum = length(Loss);
[m n] = size(Loss{1});
W = cell(viewNum,1);
for v = 1:viewNum
    Wv = zeros(m,n);
    switch reguType
        case 'hard'
            Wv(Loss{v} < lambda(v)) = 1;
            W{v} = Wv;
        case 'linear'
            idx = Loss{v} < lambda(v);
            Wv(idx) = 1 - Loss{v}(idx)./lambda(v);
            W{v} = Wv;
        case 'log'
            % lambda must be in (0,1);
            if min(lambda <= 0) || max(lambda) >= 1
                error('when the regularization type is logarithmic, the parameter lambda must be within (0,1)!\n');
            end
            eta = 1 - lambda;
            idx = Loss{v} < lambda(v);
            Wv(idx) = log(Loss{v}(idx)+eta(v))/log(eta(v));
            W{v} = Wv;           
        case 'exp'
            W{v} = (1 + exp(-lambda(v)))./(1 + exp(Loss{v} - lambda(v)));
        otherwise
            error('unknown regularization type of W!\n');
    end
end
end