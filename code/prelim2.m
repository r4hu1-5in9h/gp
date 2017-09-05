clear;
load('../data/data_BHP2.mat');

x1 = log_p;
x2 = log_y;    % Contains NaN data
y = log_q+log_p-log_y; % demand as G (budget share) instead of q like Blundell et al!

X = [ones(size(x1)) x1 x2];
b = regress(y,X);    % Removes NaN data

figure;
scatter3(x1,x2,y,'filled')
hold on
x1fit = linspace(min(x1),max(x1),10);
x2fit = linspace(min(x2),max(x2),10);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('log(Price)')
ylabel('log(Income)')
zlabel('log(Demand)')
view(50,10)
hold off

%%%%%%%%%%%%%%%%%%

X = [ones(size(x1)) x1];
b = regress(y,X);    % Removes NaN data

figure;
scatter(x1,y,'filled')
hold on
y_hat = b(1) + b(2)*x1fit;
plot(x1fit,y_hat)
xlabel('log(Price)')
ylabel('log(Demand)')
hold off

%%%%%%%%%%%%%%

X = [ones(size(x2)) x2];
b = regress(y,X);    % Removes NaN data

figure;
scatter(x2,y,'filled')
hold on
y_hat = b(1) + b(2)*x2fit;
plot(x2fit,y_hat)
xlabel('log(Income)')
ylabel('log(Demand)')
hold off