function [result,idx]=test_SPLMVC(L,k,Y,reguType)

 V=length(L);
 [n,~]=size(L(1));
%% 初始化
for v=1:V
    alpha(v)=1/V;
end

%% 循环过程
lambda=0.5;
for iter=1:30
%% 更新U
LL=zeros(n,n);   
for v=1:V
    LL=LL+alpha(v)*L{v};
end
[Vector,D,W]=eig(LL);
U=Vector(:,2:k+1);
U=U./repmat(sqrt(sum(U.^2,2)),1,k);

%% 更新alpha
for v=1:V
    F=U'*L{v}*U;
    loss(v)=trace(F) ; 
end
Loss=sum(loss);
loss1=loss/Loss;

 switch reguType
case 'equal'  %new add
         for v=1:V 
             alpha(v)=1/V; 
         end 
case  'hard' 
     for v=1:V
         if loss1(v)<=1/lambda
             alpha(v)=1;
         else
             alpha=0;
         end
     end

case 'linear'
    for v=1:V
        if   loss1(v)<=1/lambda
            alpha(v)=1-lambda*loss1(v);
        else
            alpha=0;
        end
    end
    
case 'exponential'
    for v=1:V
       alpha(v) = (1 + exp(-lambda))./(1 + exp(loss1(v) - lambda));  
    end
                                      
case 'mixture'      
  for v=1:V
    if loss1(v)<lambda/(lambda+exp(1));
        alpha(v)=1;
    elseif loss1(v)>lambda/(lambda+1);
        alpha(v)=0;
    else
        alpha(v)=log(lambda*(1/loss1(v)-1));
    end      
  end
end


end
[idx,C]=kmeans(U,k);
result = ClusteringMeasure(Y, idx);