function plot_quantile5(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block,x2_obs)
%PLOT_QUANTILE generates figures for section 5.2 without using the delta
%method

[X1_TEST,X2_TEST]=meshgrid(x1_test,x2_test);

CI_top=fmu_block+2.*sqrt(fs2_block);
CI_bottom=fmu_block-2.*sqrt(fs2_block);

fmu_block=log(fmu_block)-x1_test'+x2_test;
CI_top=log(CI_top)-x1_test'+x2_test;
CI_bottom=log(CI_bottom)-x1_test'+x2_test;

imagesc(fmu_block);

li=x2_obs-1;
ui=x2_obs+1;

figure;
scatter3(x1,x2,y,'filled');
hold on;
mesh(X1_TEST(li:ui,:),X2_TEST(li:ui,:),fmu_block(li:ui,:));
mesh(X1_TEST(li:ui,:),X2_TEST(li:ui,:),CI_top(li:ui,:));
mesh(X1_TEST(li:ui,:),X2_TEST(li:ui,:),CI_bottom(li:ui,:));
xlabel('log(Price)');
ylabel('log(Income)');
zlabel('log(Demand)');
view(-50,10);
hold off;

test=fmu_block(x2_obs,:);
test_top=CI_top(x2_obs,:);
test_bottom=CI_bottom(x2_obs,:);

figure;
plot(x1_test,test,'-','Color','b');
hold on;
plot(x1_test,test_top,'--','Color','r');
plot(x1_test,test_bottom,'--','Color','r');
xlabel('log(Price)');
ylabel('log(Demand)');
hold off;
income_slice=exp(x2_test(x2_obs));
disp(income_slice)






end