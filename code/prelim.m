clear;
load('../data/data_BHP2.mat');

x1 = log_p;
x2 = log_y;    % Contains NaN data
y = log_q;

x1fit = linspace(min(x1),max(x1),100);
x2fit = linspace(min(x2),max(x2),100);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);

%%%%%%%%%%%%%%%%%%%

X = [ones(size(x1)) x1 x2];
[b,bint] = regress(y,X);    % Removes NaN data

YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
Y_LOW = bint(1,1)+bint(2,1)*X1FIT + bint(3,1)*X2FIT;
Y_UPP = bint(1,2)+bint(2,2)*X1FIT + bint(3,2)*X2FIT;

figure;
scatter3(x1,x2,y,'filled');
hold on;
mesh(X1FIT,X2FIT,YFIT);
mesh(X1FIT,X2FIT,Y_LOW);
mesh(X1FIT,X2FIT,Y_UPP);
xlabel('log(Price)');
ylabel('log(Income)');
zlabel('log(Demand)');
view(50,10);
hold off;

%%%%%%%%%%%%%%%%%%

X = [ones(size(x1)) x1];
[b bint] = regress(y,X);    % Removes NaN data

figure;
scatter(x1,y,'filled');
hold on;
y_hat = b(1) + b(2)*x1fit;
y_low = bint(1,1)+bint(2,1)*x1fit;
y_upp = bint(1,2)+bint(2,2)*x1fit;
plot(x1fit,y_hat);
plot(x1fit,y_low,'r--');
plot(x1fit,y_upp,'r--');
xlabel('log(Price)');
ylabel('log(Demand)');
hold off;

%%%%%%%%%%%%%%

X = [ones(size(x2)) x2];
[b bint] = regress(y,X);    % Removes NaN data

figure;
scatter(x2,y,'filled');
hold on;
y_hat = b(1) + b(2)*x2fit;
y_low = bint(1,1)+bint(2,1)*x2fit;
y_upp = bint(1,2)+bint(2,2)*x2fit;
plot(x2fit,y_hat);
plot(x2fit,y_low,'r--');
plot(x2fit,y_upp,'r--');
xlabel('log(Income)');
ylabel('log(Demand)');
hold off;