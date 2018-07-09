disp('avg_interference_noise')
figure
grid on
hold on
title('worst case R factor')
xlabel('average noise power level')
ylabel(('average worst case R factor'))
plot((-60-(((1:z)-1)*10)),average_R_factor,'r-')

plot((-60-(((1:z)-1)*10)),average_R_factor_f_1)
%  plot((-60-(((1:z)-1)*10)),average_R_factor_f_11,'k-')
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('average R factor')
xlabel('average noise power level')
ylabel(('average R factor'))
plot((-60-(((1:z)-1)*10)),average_R_factor_f2_1)


plot((-60-(((1:z)-1)*10)),average_R_factor_f2_1_PSR,'r-')
% plot((-60-(((1:z)-1)*10)),average_R_factor_f2_11,'k-')
legend('Gale-Shapley','Boston')

figure
grid on
hold on
title('average percentage of cells utilized')
xlabel('average noise power level')
ylabel(('average percentage of cells utilized'))
plot((-60-(((1:z)-1)*10)),average_utilize_f2,'r-')


plot((-60-(((1:z)-1)*10)),average_utilize_f)
% plot((-60-(((1:z)-1)*10)),average_utilize_f1,'k-')
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('average excution time')
xlabel('average noise power level')
ylabel(('average excution time'))
plot((-60-(((1:z)-1)*10)),avg_timing2,'r-')

plot((-60-(((1:z)-1)*10)),avg_timing)
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('average capacity per user')
xlabel('average noise power level')
ylabel(('average capacity per user'))
plot((-60-(((1:z)-1)*10)),rate_PSR,'r-')

plot((-60-(((1:z)-1)*10)),rate)
% plot((-60-(((1:z)-1)*10)),avg_timing1,'k-')
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('maximum number of iterations')
xlabel('average noise power level')
ylabel(('maximum number of iterations'))
plot((-60-(((1:z)-1)*10)),max_num_iter_PSR-1,'r-')

plot((-60-(((1:z)-1)*10)),max_num_iter-1)%(-60-(((1:z)-1)*10))
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('average number of users in macro cells')
xlabel('average noise power level')
ylabel(('average number of users in macro cells'))
plot((-60-(((1:z)-1)*10)),avg_utilization_macro_PSR,'r-')

plot((-60-(((1:z)-1)*10)),avg_utilization_macro)%(-60-(((1:z)-1)*10))
legend('Boston','Gale-Shapley')

figure
grid on
hold on
title('average number of users in small cells')
xlabel('average noise power level')
ylabel(('average number of users in small cells'))
plot((-60-(((1:z)-1)*10)),avg_utilization_SBS_PSR,'r-')

plot((-60-(((1:z)-1)*10)),avg_utilization_SBS)%(-60-(((1:z)-1)*10))
legend('Boston','Gale-Shapley')
