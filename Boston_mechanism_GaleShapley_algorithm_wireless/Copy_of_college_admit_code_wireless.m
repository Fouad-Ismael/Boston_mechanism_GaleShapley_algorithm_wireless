clear all
warning off
%% general parameters utilized
L = 12;%total number of users
L_M = 2;%number of macro cell
wsxc = 0;
t_M = 15;%scale factor for macro-cell backhaul delay(in ms)
t_S = 40;%scale factor for small-cell backhaul delay(in ms)
k_pareto = 1/1.16;% pareto index
sigma_M=t_M*k_pareto;
sigma_S=t_S*k_pareto;
delay = zeros(1,L);
delay_PSR = zeros(1,L);
number_packets = 256;%number of bits in packet
propagation_loss=3;
path_loss_exponent = 3;
average_intrference_noise = (10^(-90/10));% in milliwats
slot_time = 20;%(in ms)
epson3 = 177.3;
%% main code
tic
for z = 9:9
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    N = 20+(z-1)*10;%number of users(varies by step 10 )
   % qouta(:,1:L_M) = 0.25*N*ones(size(qouta(:,1:L_M)));
    rewq = [];
    average_R_factor_f = [];
    average_R_factor_f2 = [];   

    for zz = 1:0.2*10^3
        zz
        z
        qouta = 4*ones(1,L);%uplink and downlink qouta for each cell

        % the part of calculating the position%
        pos_ue_cell = [randi(1001,1,N+L);randi(2001,1,N+L)];% assume that the area of interest is 1000*2000
        % so here we put the coordinates of all users and cells.
        % i utilize the convention of having the users from 1 to N and
        % cells from N+1 to N+L.
        
        %modifying the position of the macro-cells to be at the center of
        %each cell:
        % X&Y coordinates for first station
        pos_ue_cell(1,N+1)=501;
        pos_ue_cell(2,N+1)=501;
        % X&Y coordinates for second station
        pos_ue_cell(1,N+2)=501;
        pos_ue_cell(2,N+2)=1501;
        % calculating distances
        distance1 = zeros(L+N,N+L);
        for i = 1:N+L
            distance1(i,:) = (((pos_ue_cell(1,:)-pos_ue_cell(1,i)).^2+(pos_ue_cell(2,:)-pos_ue_cell(2,i)).^2).^0.5);
            % here we apply euclidean distance between all points (users or
            % cells), so , at the end we will have distances between cells
            % and each other, users and each other or between users and cells.
        end
        near_cell = pos_ue_cell(2,1:N)>1001; % by this we determine which users are in the second cell
        % so that they can be served by its macro base station
        qouta(1:2)= [length(find(~near_cell)) length(find(near_cell))];% modification of the qouta of macro-cell
        % according to how many users in each cell .
        % the guaranteed delay for

        BackHual_Delay=[gprnd(k_pareto,sigma_M,sigma_M/k_pareto,1,L_M) gprnd(k_pareto,sigma_S,sigma_S/k_pareto,1,L-L_M)];
        Wireless_Delay= ((qouta-1)*slot_time);
        maximum_delay = Wireless_Delay + BackHual_Delay;
%         H = (maximum_delay>epson3);
        %epson3 = 177.3;./([qouta(:,1:L_M)-1  ones(1,L-L_M)]) .*([15*ones(1,L_M)  ones(1,L-L_M)])+[gprnd(k_pareto,sigma_M,0,1,L_M) gprnd(k_pareto,sigma_S,0,1,L-L_M)]
%         if N<=70     
%         x=find(maximum_delay>3.3*epson3);
%                 maximum_delay(x)=3.3*epson3;
%         else
%             x=find(maximum_delay>3.3*epson3);
%                 maximum_delay(x)=3.3*epson3;
%         end
%         
        is_NOT_associated = ones(2,N);
        
        %for now; we assume that in the channel gain matrix, coordinates from 1:N
        %are to users and coordinates from N+1 : N+L belong to cells. this goes for
        %both coordinates(x,y)
        P_U = 0.02; % power of  UE in milliwatts
        P_UE_matrix = P_U*ones(1,N);
        channel_gain = ((distance1+eps).^(-path_loss_exponent));
        % 1:N represent the users , N+1 to N+L as cells
        % so (channel_gain(1:N,N+1:N+L)) means the channel gain between
        % users and cells
        sinr_up_2 = repmat(P_UE_matrix,L,1)'.*channel_gain(1:N,N+1:N+L)/average_intrference_noise;
        
        % i'll use BPSK or QPSK modulation and this is its packet success rate
        PSR_up_2= (1-(0.5*erfc(sqrt((sinr_up_2))))).^number_packets;
         
        % construct cell preference:
        % the small cells get a separate benefit from the macro cell.
        small_cell_benfit_UP= PSR_up_2./ repmat(max(PSR_up_2(:,1:L_M)')',1,L);
        % we do this because the paper suggested that we divide by the PSR
        % of the associated macro base station, that's why I used the
        % minimum function
        cell_benfit_UP = [PSR_up_2(:,1:L_M), small_cell_benfit_UP(:,L_M+1:end)];
       
       % Id = (0.024*repmat(maximum_delay(1,:),N,1)) - (0.11*repmat((maximum_delay(1,:)-epson3).*H,N,1))
       frew = 177.3;
       mT = 150;
       sT = 0.4;
       for gg = 1: L
            if maximum_delay(1,gg)<=150
                Iddd(gg) = 0;
            else
            XXX = log10(maximum_delay(1,gg)/mT)/log10(2);
        Iddd(gg) = 25*((1+XXX^(6*sT))^(1/(6*sT))-3*(1+(XXX/3)^(6*sT))^(1/(6*sT))+2);
            end
           if maximum_delay(1,gg)<=frew
       Id1(gg) = -2.468*(10^-14)*(maximum_delay(1,gg).^6);
       Id2(gg) = 5.062*(10^-11)*(maximum_delay(1,gg).^5);
       Id3(gg) = -3.903*(10^-8)*(maximum_delay(1,gg).^4);
       Id4(gg) = 1.344*(10^-5)*(maximum_delay(1,gg).^3);
       Id5(gg) = -0.001802*(maximum_delay(1,gg).^2);
       Id6(gg) = 0.103*maximum_delay(1,gg);
       Id7(gg) = -0.1698;
       Iddddd(gg) = Id1(gg)+Id2(gg)+Id3(gg)+Id4(gg)+Id5(gg)+Id6(gg)+Id7(gg);
       Id(gg) = Iddddd(gg)+ Iddd(gg);
           else
       Id1(gg) = -2.468*(10^-14)*(frew^6);
       Id2(gg) = 5.062*(10^-11)*(frew^5);
       Id3(gg) = -3.903*(10^-8)*(frew^4);
       Id4(gg) = 1.344*(10^-5)*(frew^3);
       Id5(gg) = -0.001802*(frew^2);
       Id6(gg) = 0.103*frew;
       Id7(gg) = -0.1698;
        Iddddd(gg) = Id1(gg)+Id2(gg)+Id3(gg)+Id4(gg)+Id5(gg)+Id6(gg)+Id7(gg); 
        Id(gg) = Iddddd(gg)+ Iddd(gg);
           end
       end
         ID_Final = repmat(Id, N, 1);      
        R_factor=94.2-ID_Final-12-(15*log(1+60*(1-PSR_up_2)));
        %preference lists for both users and cells
%         [~,R_factor_index]=sort(PSR_up_2,2,'descend');
        [~,PSR_up_2_index]=sort(PSR_up_2,2,'descend');
        
        %%%%%%%college admission begins%%%%%%%%
        %preference lists for all users
        for i = 1:N
%             user_preference.([string('user'),num2str(i)]) =R_factor_index(i,:);
             user_preference_PSR.([string('user'),num2str(i)]) =PSR_up_2_index(i,:);
        end
        %dynamic list name of users per cell
        for i = 1:L
%             cell_association.([string('cell'),num2str(i)]) =[];
%             cell_association_rank.([string('cell'),num2str(i)]) =[];
                            cell_association_PSR.([string('cell'),num2str(i)]) =[];
            %                 cell_association_rank_PSR.([string('cell'),num2str(i)]) =[];
        end
        %this flag will represent whether the user still have preference list
%         has_preference = ones(1,N);
        %this flag will represent whether the user is associated or not
%         is_NOT_associated(1,:) = ones(1,N);
                     has_preference_PSR = ones(1,N);
                    is_NOT_associated_PSR = ones(1,N);
%         tic
%         while (sum(has_preference.* is_NOT_associated(1,:))>0)
%             % (no user that has a preference list and not associated
%             % exists)
%             for i = 1:N
%                 if has_preference(i)==0
%                     continue;
%                 else
%                     if is_NOT_associated(1,i) == 1
%                         application = user_preference.([string('user'),num2str(i)])(1);% submit application
%                         rank = R_factor(i,application);% application evaluated by cell
%                         if length(cell_association.([string('cell'),num2str(application)]))+1<= qouta(1,application)
%                             cell_association.([string('cell'),num2str(application)]) =[cell_association.([string('cell'),num2str(application)]) i];
%                             cell_association_rank.([string('cell'),num2str(application)]) =[cell_association_rank.([string('cell'),num2str(application)]) rank];
%                             is_NOT_associated(1,i) = 0;
%                         else
%                             [p,j] = min(cell_association_rank.([string('cell'),num2str(application)]));
%                             if rank > p
%                                 x = cell_association.([string('cell'),num2str(application)])(j);
%                                 cell_association.([string('cell'),num2str(application)])(j)=i;
%                                 cell_association_rank.([string('cell'),num2str(application)])(j)=rank;
%                                 is_NOT_associated(1,i) = 0;
%                                 user_preference.([string('user'),num2str(x)]) = user_preference.([string('user'),num2str(x)])(2:end);
%                                 has_preference(x)=~isempty(user_preference.([string('user'),num2str(x)]));
%                                 is_NOT_associated(1,x) = 1;
%                             else
%                                 user_preference.([string('user'),num2str(i)]) = user_preference.([string('user'),num2str(i)])(2:end);
%                                 has_preference(i)=~isempty(user_preference.([string('user'),num2str(i)]));
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         toc
        
                quota_actual = zeros(size(qouta));
                
                    while (sum(has_preference_PSR.* is_NOT_associated_PSR(1,:))>0)
                        for i = 1:L
                          cell_associate_boston.(([string('cell'),num2str(i)])) = [];
                          cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];            
                        end
                        for i_PSR = 1:N
                            if has_preference_PSR(i_PSR)==0
                                continue;
                            else
                                if is_NOT_associated_PSR(1,i_PSR) ==1
                                    application_PSR = user_preference_PSR.([string('user'),num2str(i_PSR)])(1);
                                    rank_PSR = R_factor(i_PSR,application_PSR);
                                    if length(cell_associate_boston.([string('cell'),num2str(application_PSR)]))+1<= (qouta(1,application_PSR)-quota_actual(1,application_PSR))
                                        cell_associate_boston.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston.([string('cell'),num2str(application_PSR)]) i_PSR];
                                        cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) rank_PSR];
                                        is_NOT_associated_PSR(1,i_PSR) = 0;
                                    elseif (qouta(1,application_PSR)-quota_actual(1,application_PSR))>0
                                                   [p,j] = min(cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]));
                                         if rank_PSR > p
                                            x = cell_associate_boston.([string('cell'),num2str(application_PSR)])(j);
                                            cell_associate_boston.([string('cell'),num2str(application_PSR)])(j)=i_PSR;
                                            cell_associate_boston_rank.([string('cell'),num2str(application_PSR)])(j)=rank_PSR;
                                            is_NOT_associated_PSR(1,i_PSR) = 0;
                                            user_preference_PSR.([string('user'),num2str(x)]) = user_preference_PSR.([string('user'),num2str(x)])(2:end);
                                            has_preference_PSR(x)=~isempty(user_preference_PSR.([string('user'),num2str(x)]));
                                            is_NOT_associated_PSR(1,x) =1;
                                        else
                                            user_preference_PSR.([string('user'),num2str(i_PSR)]) = user_preference_PSR.([string('user'),num2str(i_PSR)])(2:end);
                                            has_preference_PSR(i_PSR)=~isempty(user_preference_PSR.([string('user'),num2str(i_PSR)]));
                                        end
                                    else
%                                         [p,j] = min(cell_association_rank_PSR.([string('cell'),num2str(application_PSR)]));
%                                         if rank_PSR > p
%                                             x = cell_association_PSR.([string('cell'),num2str(application_PSR)])(j);
%                                             cell_association_PSR.([string('cell'),num2str(application_PSR)])(j)=i_PSR;
%                                             cell_association_rank_PSR.([string('cell'),num2str(application_PSR)])(j)=rank_PSR;
%                                             is_NOT_associated_PSR(1,i_PSR) = 0;
%                                             user_preference_PSR.([string('user'),num2str(x)]) = user_preference_PSR.([string('user'),num2str(x)])(2:end);
%                                             has_preference_PSR(x)=~isempty(user_preference_PSR.([string('user'),num2str(x)]));
%                                             is_NOT_associated_PSR(1,x) = 1;
%                                         else
                                            user_preference_PSR.([string('user'),num2str(i_PSR)]) = user_preference_PSR.([string('user'),num2str(i_PSR)])(2:end);
                                            has_preference_PSR(i_PSR)=~isempty(user_preference_PSR.([string('user'),num2str(i_PSR)]));
%                                         end
                                     end
                                end
                            end
                        end
                        for i = 1:L
                cell_association_PSR.([string('cell'),num2str(i)]) = [cell_association_PSR.([string('cell'),num2str(i)]), cell_associate_boston.(([string('cell'),num2str(i)]))];
                  quota_actual(i) =  length(cell_association_PSR.([string('cell'),num2str(i)]));
                      
                        end
                    end
                    
        rewq_1=[];
         rewq_1_PSR=[];
%         rewq_psr=[];
%         for i = 1:L
%             if isempty(cell_association.([string('cell'),num2str(i)]))
%                 continue
%             else
%                 %                     if i > L_M
%                 New_Wireless_Delay(i)= slot_time*(length(cell_association.([string('cell'),num2str(i)]))-1);
%                 delay(i) = New_Wireless_Delay(i) + BackHual_Delay(i);
%                 %                     else
%                 %                       delay(i) =maximum_delay(i) + slot_time*(length(cell_association.([string('cell'),num2str(i)]))-1) - 15*slot_time + (N^2/1000)*slot_time;
%                 %                     end
%                 h = delay(i)>epson3;
% %                 if N<=70
% %                                    if (delay(i)>3.3*epson3)
% %                                        delay(i)=3.3*epson3;
% %                                    end
% %                 else
% %                     if (delay(i)>3.3*epson3)
% %                         delay(i)=3.3*epson3;
% %                     end
% %                 end
% 
%             if delay(i) <= 150
%                 Idd(i)=0;
%             else
%             X = log10(delay(i)/mT)/log10(2);
%             Idd(i) = 25*((1+X^(6*sT))^(1/(6*sT))-3*(1+(X/3)^(6*sT))^(1/(6*sT))+2);
%             end
%       if delay(i)<= frew
%        Id1i = -2.468*(10^-14)*delay(i)^6;
%        Id2i = 5.062*(10^-11)*delay(i)^5;
%        Id3i = -3.903*(10^-8)*delay(i)^4;
%        Id4i = 1.344*(10^-5)*delay(i)^3;
%        Id5i = -0.001802*delay(i)^2;
%        Id6i = 0.103*delay(i)^1;
%        Id7i = -0.1698;
%         Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
%         Idi(i) =  Idii(i) + Idd(i);
%       else
%        Id1i = -2.468*(10^-14)*(frew)^6;
%        Id2i = 5.062*(10^-11)*(frew)^5;
%        Id3i = -3.903*(10^-8)*(frew)^4;
%        Id4i = 1.344*(10^-5)*(frew)^3;
%        Id5i = -0.001802*(frew)^2;
%        Id6i = 0.103*(frew)^1;
%        Id7i = -0.1698;
%         Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
%         Idi(i) =  Idii(i) + Idd(i);
%      end
%                 for jj = 1:length(cell_association.([string('cell'),num2str(i)]))
% %                     rewq = [rewq  R_factor(cell_association.([string('cell'),num2str(i)])(jj),i)];
%                     
%                     rewq_1= [rewq_1  94.2-Idi(i)-12-(15*(log(1+60*(1-PSR_up_2(cell_association.([string('cell'),num2str(i)])(jj),i)))))];
%                 end
%             end
%         end
                     for i = 1:L
                        if isempty(cell_association_PSR.([string('cell'),num2str(i)]))
                            continue
                        else
        %                    if i <= L_M
        %                    delay_PSR(i) =maximum_delay(i) + slot_time*(length(cell_association_PSR.([string('cell'),num2str(i)]))-1) - 20*slot_time; %+ (N/20)*slot_time;
        %                    else
                         New_Wireless_Delay(i)= slot_time*(length(cell_association_PSR.([string('cell'),num2str(i)]))-1);

                              delay_PSR(i) = New_Wireless_Delay(i) + BackHual_Delay(i);
        %                    end
        %                     if delay_PSR(i)>4*epson3
        %                        delay_PSR(i)=4*epson3;
        %                    end
         if delay_PSR(i) <= 150
                Idd(i)=0;
            else
            X = log10(delay_PSR(i)/mT)/log10(2);
            Idd(i) = 25*((1+X^(6*sT))^(1/(6*sT))-3*(1+(X/3)^(6*sT))^(1/(6*sT))+2);
            end
      if delay_PSR(i)<= frew
       Id1i = -2.468*(10^-14)*delay_PSR(i)^6;
       Id2i = 5.062*(10^-11)*delay_PSR(i)^5;
       Id3i = -3.903*(10^-8)*delay_PSR(i)^4;
       Id4i = 1.344*(10^-5)*delay_PSR(i)^3;
       Id5i = -0.001802*delay_PSR(i)^2;
       Id6i = 0.103*delay_PSR(i)^1;
       Id7i = -0.1698;
        Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
        Idi(i) =  Idii(i) + Idd(i);
      else
       Id1i = -2.468*(10^-14)*(frew)^6;
       Id2i = 5.062*(10^-11)*(frew)^5;
       Id3i = -3.903*(10^-8)*(frew)^4;
       Id4i = 1.344*(10^-5)*(frew)^3;
       Id5i = -0.001802*(frew)^2;
       Id6i = 0.103*(frew)^1;
       Id7i = -0.1698;
        Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
        Idi(i) =  Idii(i) + Idd(i);
     end
                            for jj = 1:length(cell_association_PSR.([string('cell'),num2str(i)]))
%                                 rewq_psr = [rewq_psr  R_factor(cell_association_PSR.([string('cell'),num2str(i)])(jj),i)];
                                rewq_1_PSR= [rewq_1_PSR  94.2-Idi(i)-12-15*(log(1+60*(1-PSR_up_2(cell_association_PSR.([string('cell'),num2str(i)])(jj),i))))];
                            end
                        end
                    end
        %
         % average_R_factor_f(zz) = mean((rewq_psr));
%     if sum(rewq_1<0)
%         disp('negative benifite')
%         wsxc = wsxc+1;
%     end
    average_R_factor_f2(zz) = mean((rewq_1_PSR));
    % average_R_factor_f2_PSR(zz) = mean((rewq_1_PSR));

    end
   
wsxc
% average_R_factor(z) = mean(average_R_factor_f);
% average_R_factor_f_1(z)=mean((rewq));
average_R_factor_f2_1(z) = mean(average_R_factor_f2);
% average_R_factor_f2_1_PSR(z) = mean(average_R_factor_f2_PSR);
end
toc
%% results plotted
% figure
%plot(20+((1:7)-1)*10,average_R_factor,'r-')
%hold on
% plot(20+((1:z)-1)*10,average_R_factor_f_1)
figure
plot(20+((1:z)-1)*10,average_R_factor_f2_1)
%hold on
%plot(20+((1:7)-1)*10,average_R_factor_f2_1_PSR,'r-')
