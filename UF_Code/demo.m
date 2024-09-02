clc
clear all
close all
%% 添加路径
addpath(genpath('Algorithm'))  
addpath(genpath('data_set'))
addpath(genpath('function'))
%% 载入数据
Files = dir(fullfile('..\Code\data_set', '*.mat'));   
Max_datanum = length(Files);   % 数据集的个数
t = 1;                % 防止加入新的算法
 for re_num = 1:1 % 第一个循环，控制重复次数   
 for data_num = 1:1 % 第二个循环，控制测试的数据集的序号   
 Dname = Files(data_num).name; % data name 
  disp(['***********The test data name is: ***' num2str(data_num) '***'  Dname '****************'])
 load(Dname);
 V=length(X);           
 for v=1:V
    X{v}=double(full(X{v}));               
 end
%% 计算相似度矩阵 
v = length(X);
k = max(Y);
Fs = cell(v, 1);
n = length(Y); 
tic;
for i = 1:v
    for  j = 1:n
        X{i}(j,:) = ( X{i}(j,:) - mean( X{i}(j,:) ) ) / std( X{i}(j,:) ) ;
    end
end 
Ls = cell(1,v);
for idx = 1:v
    A0 = constructW_PKN(X{idx}',10);
    A0 = A0-diag(diag(A0));
    A10 = (A0+A0')/2;
    D10 = diag(sum(A10));
    Ls{idx} = D10 - A10;
    A{idx}=A10;
end
 %% test algorithms 
 do_ED=0;  %% 对应文章MSE
 do_PRL_Root=0; % 对应文章AMGL  
 do_UF=0;  %% 联合框架
 do_MALG_Hard=0;
 do_MALG_LSW=0; 
 do_MALG_EF=0; 
 do_SPLMVC=0;
 do_NE=0;
 do_SWL_LSW=0;    
 do_SWL_CSL=1;
 %% Test Exponential Decay
  if do_ED == 1
     addpath(genpath('Algorithm\ED'))  
     q = .01:.1:1;
     gamma = 1./q;
     results = cell(length(gamma),1);
     lbd_prime = zeros(length(q), 1);
     for ii  = 1:length(gamma)
     rng('default');
     [results{ii},label,fs_star] = Test_ED(Ls, k, Y, gamma(ii));
     lbd_prime(ii) = -sum( (fs_star./q(ii)).^(1/(q(ii)-1)) )^(q(ii)-1);
     end  
  end
  %% Test Unified Learning 
    if do_UF == 1
       addpath(genpath('Algorithm\UF'))  
%           q = -(.01:.1:1);
%           lambda = -[-(1:5) -(5:2:10) -(10:5:30)];
%           % q = (.01:.1:1);
%           % lambda = [-(1:5) -(5:2:10) -(10:5:30)];
            q=-0.910;
            lambda=0.073; 
            results = cell(length(lambda), length(q)); 
            for ii  = 1:length(lambda)
            for jj = 1:length(q)
            rng('default');
            [results{ii,jj},label] = Test_UF(Ls, k, Y, lambda(ii), q(jj));  
            end
         end
    end    
  %% Test MALG_Hard (lambda>0)(里面涉及随机变量，需要运行多次，求平均结果)  
   if do_MALG_Hard == 1
      addpath(genpath('Algorithm\SP_Hard'))
     [results_hard(re_num,:),y,S,W ] = MALG(A,Y,1,1,'hard');     
   end 
%      hard_mean=mean(results_hard);   
   
   %% Test MALG_LSW(lambda>0) (里面涉及随机变量，需要运行多次，求平均结果)  
   if do_MALG_LSW == 1
      addpath(genpath('Algorithm\SP_Linear'))
      [results_linear(re_num,:),y,S,W ] = MALG(A,Y,1,1,'linear'); 
   end  
%    linear_mean=mean(results_linear); 
   
   %% Test MALG_EF (lambda>0)(里面涉及随机变量，需要运行多次，求平均结果)  
    if do_MALG_EF == 1
       addpath(genpath('Algorithm\SP_exp'))
       [results_exp(re_num,:),y,S,W ] = MALG(A,Y,1,1,'exp');
    end 
    %  Exp_mean=mean(results_exp);         
    
    %% Test_SPL_MVC; (里面涉及随机变量，需要运行多次，求平均结果)     
    if do_SPLMVC == 1
       addpath(genpath('Algorithm\SPLMVC'))
       [results(re_num,:),idx] = test_SPLMVC(Ls,k,Y,'linear'); %%'equal' 'hard','linear','exponential','mixture'
    end     
%      SPL_mean=mean(results);

   %% Test Negative-entropy(lambda>0)对应于文章中的MALG算法中lambda=1,3,5,7,9;
    if do_NE == 1 
       addpath(genpath('Algorithm\NE'))
       lambda = 1:2:10;  %% lambda调节
       results = cell(length(lambda),1);  
       for ii = 1:length(lambda)
       rng('default');
       [results{ii},idx] = Test_NE(Ls, k, Y, lambda(ii));
       end
    end  
    
    %% Test PRL_Root(文章中的AMGL)   
    if do_PRL_Root == 1
       addpath(genpath('Algorithm\PRL_Root'))
      %  q = -(.01:.2:1); 
      %  p = q./(q-1); 
         q= -0.11;
         p=q./(q-1); 
       results = cell(length(p),1); 
       lbd_prime = zeros(length(p), 1);
       for ii  =1:length(p)
       rng('default');
       [results{ii}, label,fs_star] = Test_PRL_Root(Ls, k, Y, p(ii),'root'); %% manner includes 'root','log','capped',sata=3 represenst capped degree,it may tune.
       lbd_prime(ii) = sum( (-fs_star./q(ii)).^(1/(q(ii)-1)) )^(q(ii)-1);
       end
    end
    
    %% Test SWL_LSW (p趋近于0) 
    if do_SWL_LSW ==1; 
        addpath(genpath('Algrithm\SWL_LSW'))
        rng('default');
        [results,label,fs_star] = Test_SWL_LSW(Ls, k, Y,'log'); %% manner includes 'root','log','capped',sata=3 represenst capped degree,it may tune.
    end
    
    %% Test SWL_CSL （与p无关）
    if  do_SWL_CSL == 1
        addpath(genpath('Algorithm\SWL_CSL'))
        rng('default');
        [results,label, fs_star] = Test_SWL_CSL(Ls, k, Y,'capped',3); %% manner includes 'root','log','capped',sata=3 represenst capped degree,it may tune.
    end   
    results
 end
%       Result(re_num,:)=results;

 end

        
  
      
      