disp('num_users_90noise')
%20+((1:z)-1)*10
figure
hold on
grid on
title('worst case R factor')
xlabel('number of users')
ylabel(('average worst case R factor'))
plot(20+((1:z)-1)*10,average_R_factor,'r-')

plot(20+((1:z)-1)*10,average_R_factor_f_1)
%  plot(20+((1:z)-1)*10,average_R_factor_f_11,'k-')
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average R factor')
xlabel('number of users')
ylabel(('average R factor'))
plot(20+((1:z)-1)*10,average_R_factor_f2_1)


plot(20+((1:z)-1)*10,average_R_factor_f2_1_PSR,'r-')
% plot(20+((1:z)-1)*10,average_R_factor_f2_11,'k-')
legend('Gale-Shapley','Boston')

figure
hold on
grid on
title('average percentage of cells utilized')
xlabel('number of users')
ylabel(('average percentage of cells utilized'))
plot(20+((1:z)-1)*10,average_utilize_f2,'r-')


plot(20+((1:z)-1)*10,average_utilize_f)
legend('Boston','Gale-Shapley')

% plot(20+((1:z)-1)*10,average_utilize_f1,'k-')

figure
hold on
grid on
title('average excution time')
xlabel('number of users')
ylabel(('average excution time'))
plot(20+((1:z)-1)*10,avg_timing2,'r-')

plot(20+((1:z)-1)*10,avg_timing)
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average capacity per user')
xlabel('number of users')
ylabel(('average capacity per user'))
plot(20+((1:z)-1)*10,rate_PSR,'r-')

plot(20+((1:z)-1)*10,rate)
legend('Boston','Gale-Shapley')

% plot(20+((1:z)-1)*10,avg_timing1,'k-')
figure
hold on
grid on
title('maximum number of iterations')
xlabel('number of users')
ylabel(('maximum number of iterations'))
plot(20+((1:z)-1)*10,max_num_iter_PSR-1,'r-')

plot(20+((1:z)-1)*10,max_num_iter-1)%20+((1:z)-1)*10
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average number of users in macro cells')
xlabel('number of users')
ylabel(('average number of users in macro cells'))
plot(20+((1:z)-1)*10,avg_utilization_macro_PSR,'r-')


plot(20+((1:z)-1)*10,avg_utilization_macro)%20+((1:z)-1)*10
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average number of users in small cells')
xlabel('number of users')
ylabel(('average number of users in small cells'))
plot(20+((1:z)-1)*10,avg_utilization_SBS_PSR,'r-')

plot(20+((1:z)-1)*10,avg_utilization_SBS)
legend('Boston','Gale-Shapley')
