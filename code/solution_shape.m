%% visualize q_i^{bounds}
clear;
startup;

q_i=linspace(0,1,100);
quad=q_i.*(1-q_i);
w_i_price=0.0675;
w_i_inc=0.25;
affine=w_i_price+w_i_inc.*q_i;

plot(q_i,quad);
hold on;
plot(q_i,affine);
[q_i_lb,q_i_ub] = get_bounds(w_i_price,w_i_inc);
xlabel('q_i');
legend('RHS','LHS');
vline(q_i_lb,'--r','q_i^{lb}');
vline(q_i_ub,'--r','q_i^{ub}');
hold off;

%% Compare optimization methods on subset of data
clear all;
startup;
load('../data/data_BHP2.mat');
data=[share log_p log_y];
ind=find(log_y<log(12000)); % locally giffen consumers. budget share stats
data2=data(ind,:);
hist(data2(:,1),20);
mean(data2(:,1))

n_obs=100;
rng(3);
sample=datasample(data,n_obs,'Replace',false);
y=sample(:,1);
x1=sample(:,2);
x2=sample(:,3);

len=2;
quantiles = [0.50]; %[0.25 0.50 0.75];
optimizers = {'lbfgs'};  %{'sd','bb','pr','lbfgs'};
inf_string='infEP_shape';

for idq = 1:length(quantiles)
    
    quantileToPredict=quantiles(idq);
    
    tic;
    
    figure;
    for ido = 1:length(optimizers)
        
        optimizer=optimizers{ido};
        fprintf('quantile: %3.2f optimizer: %s\n',quantileToPredict,optimizer);
        
        [~,~,~,~,~,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
        
        param_path=output.trace.param_path;
        n_iter=output.iterations;

        subplot(2,2,ido);
        scatter3(param_path(1,:),param_path(2,:),param_path(4,:),'filled');
        line(param_path(1,:),param_path(2,:),param_path(4,:));
        title(strcat(optimizers{ido},': ',num2str(n_iter)));
        xlabel('log(\gamma_{price})');
        ylabel('log(\gamma_{inc})');
        zlabel('log(\sigma)');
        view(30,30);
        
    end
    
    toc;
    
end

%% run infEP then infEP_shape
clear all;
startup;
load('../data/data_BHP2.mat');
data=[share log_p log_y];

n_obs=100;
rng(3);
sample=datasample(data,n_obs,'Replace',false);
y=sample(:,1);
x1=sample(:,2);
x2=sample(:,3);

len=2;
quantiles = [0.50]; %[0.25 0.50 0.75];
optimizers = {'lbfgs'};
inf_string='infEP';
inf_string2='infEP_shape';

for idq = 1:length(quantiles)
    
    quantileToPredict=quantiles(idq);
    
    tic;
    
    figure;
    for ido = 1:length(optimizers)
        
        optimizer=optimizers{ido};
        fprintf('quantile: %3.2f optimizer: %s\n',quantileToPredict,optimizer);
        
        [~,~,~,~,~,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
        
        param_path=output.trace.param_path;
        n_iter=output.iterations;

        subplot(1,2,ido);
        scatter3(param_path(1,:),param_path(2,:),param_path(4,:),'filled');
        line(param_path(1,:),param_path(2,:),param_path(4,:));
        title(strcat(optimizers{ido},': ',num2str(n_iter)));
        xlabel('log(\gamma_{price})');
        ylabel('log(\gamma_{inc})');
        zlabel('log(\sigma)');
        view(30,30);
        
        [~,~,~,~,~,output]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string2);
        
        param_path=output.trace.param_path;
        n_iter=output.iterations;

        subplot(1,2,1+ido);
        scatter3(param_path(1,:),param_path(2,:),param_path(4,:),'filled');
        line(param_path(1,:),param_path(2,:),param_path(4,:));
        title(strcat(optimizers{ido},': ',num2str(n_iter)));
        xlabel('log(\gamma_{price})');
        ylabel('log(\gamma_{inc})');
        zlabel('log(\sigma)');
        view(30,30);
        
    end
    
    toc;
    
end
%% Implement chosen method on full data set

clear;
load('../data/data_BHP2.mat');
y=share;
x1=log_p;
x2=log_y;

len=10;
optimizer='lbfgs';
inf_string='infEP_shape';
quantiles = [0.50]; %run one by one for some reason

%note: this takes many hours to run
train=false;


if train

    tic;
    %ticBytes(gcp);
    
    for ido = 1:numel(quantiles)
        quantileToPredict = quantiles(ido);
        [~,~,~,~,~,~]=save_quantile(quantileToPredict,len,optimizer,x1,x2,y,inf_string);
    end
    
    %tocBytes(gcp);
    toc;

end

%% Visualize quantiles of gasoline demand

clear;
load('../data/data_BHP2.mat');
y=log_q;
x1=log_p;
x2=log_y;

%debugging plots
load('../data/results_0.5_100_lbfgs_restricted.mat');
plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block);
title('tau=0.50');

%%
% =============================
%  Subfunctions
% =============================
%%

%% loss_vis.m
% visualize different loss functions
%
% <include>loss_vis.m</include>
%
%%

%% prelim.m
% visualize mean gasoline demand
%
% <include>prelim.m</include>
%
%%

%% save_quantile.m
% train site parameters, tune hyperparameters, and save output
%
% <include>save_quantile.m</include>
%
%%

%% plot_quantile.m
% visualize quantile gasoline demand with credible interval
%
% <include>plot_quantile.m</include>
%
%%