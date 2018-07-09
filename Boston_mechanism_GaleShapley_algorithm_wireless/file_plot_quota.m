disp('qouta')
%(1:z)+1
figure
hold on
grid on
title('worst case R factor')
xlabel('quota')
ylabel(('average worst case R factor'))
plot((1:z)+1,average_R_factor,'r-')

plot((1:z)+1,average_R_factor_f_1)
%  plot((1:z)+1,average_R_factor_f_11,'k-')
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average R factor')
xlabel('quota')
ylabel(('average R fator'))

plot((1:z)+1,average_R_factor_f2_1)


plot((1:z)+1,average_R_factor_f2_1_PSR,'r-')
legend('Gale-Shapley','Boston')

% plot((1:z)+1,average_R_factor_f2_11,'k-')

figure
hold on
grid on
title('average percentage of cells utilized')
xlabel('quota')
ylabel(('average percentage of cells utilized'))
plot((1:z)+1,average_utilize_f2,'r-')


plot((1:z)+1,average_utilize_f)
legend('Boston','Gale-Shapley')

% plot((1:z)+1,average_utilize_f1,'k-')

figure
hold on
grid on
title('average excution time')
xlabel('quota')
ylabel(('average excution time'))
plot((1:z)+1,avg_timing2,'r-')


plot((1:z)+1,avg_timing)
legend('Boston','Gale-Shapley')

figure
hold on
grid on
title('average rate per user')
xlabel('quota')
ylabel(('average rate per user'))
plot((1:z)+1,rate_PSR,'r-')


plot((1:z)+1,rate)
legend('Boston','Gale-Shapley')

% plot((1:z)+1,avg_timing1,'k-')
figure
hold on
grid on
title('maximum number of iterations')
xlabel('quota')
ylabel(('maximum number of iterations'))
plot((1:z)+1,max_num_iter_PSR-1,'r-')


plot((1:z)+1,max_num_iter-1)%(1:z)+1
legend('Boston','Gale-Shapley')
figure
hold on
grid on
plot((1:z)+1,avg_utilization_macro_PSR,'r-')


plot((1:z)+1,avg_utilization_macro)%(1:z)+1
title('average number of users in macro cells')
legend('Boston','Gale-Shapley')

xlabel('quota')
ylabel('average number of users in macro cells')
figure
hold on
grid on
title('average number of users in small cells')
xlabel('quota')
ylabel('average number of users in small cells')
plot((1:z)+1,avg_utilization_SBS_PSR,'r-')

plot((1:z)+1,avg_utilization_SBS)
legend('Boston','Gale-Shapley')
