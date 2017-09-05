function plot_quantile3(x1,x2,y,x1_test,x2_test,fmu_block,fs2_block)
%PLOT_QUANTILE generates figures for section 5.2 without using the delta
%method

CI_top=fmu_block+2.*sqrt(fs2_block);
CI_bottom=fmu_block-2.*sqrt(fs2_block);

fmu_block=log(fmu_block)-x1_test'+x2_test;
CI_top=log(CI_top)-x1_test'+x2_test;
CI_bottom=log(CI_bottom)-x1_test'+x2_test;

figure;
scatter3(x1,x2,y,'filled');
hold on;
mesh(x1_test,x2_test,fmu_block);
mesh(x1_test,x2_test,CI_top);
mesh(x1_test,x2_test,CI_bottom);
xlabel('log(Price)');
ylabel('log(Income)');
zlabel('log(Demand)');
view(50,10);
hold off;

end