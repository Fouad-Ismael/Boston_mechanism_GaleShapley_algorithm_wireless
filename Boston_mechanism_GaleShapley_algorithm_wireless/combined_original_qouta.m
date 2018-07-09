clear all
warning off
% college admission based on paper of Dr/Walid Saad "A College Admissions Game for Uplink User Association
%in Wireless Small Cell Networks"
% boston mechanism is based on 
%% general parameters utilized

L = 12;%total number of cells
L_M = 2;%number of macro cell
t_M = 15;%scale factor for macro-cell backhaul delay(in ms)
t_S = 40;%scale factor for small-cell backhaul delay(in ms)
k_pareto = 1/1.16;% pareto index
sigma_M=t_M*k_pareto;
sigma_S=t_S*k_pareto;
delay = zeros(1,L); % delay in the cells in the college admission
delay_PSR = zeros(1,L);% delay in the cells in the boston mechanism
number_packets = 160;%number of bits in packet
propagation_loss=3;
path_loss_exponent = 3;
slot_time = 20;%(in ms)
epson3 = 177.3;
%       qouta = (4)*ones(1,L);%uplink and downlink qouta for each cell
         N = 60;
      average_intrference_noise = (10^((-90/10)));% in milliwats
num_runs = 500;
%% main code
tic
for z = 1:9

    %   Detailed explanation goes here
%        N = 20+(z-1)*10;%number of users(varies by step 10 )
%     rewq = [];
    qouta = (z+1)*ones(1,L);%uplink and downlink qouta for each cell
%  average_intrference_noise = (10^(((-60-(z-1)*10)/10)));% in milliwats

fdsa = [zeros(num_runs,45)];
fdsa2 = [zeros(num_runs,45)];
cell_utilize_tot = [zeros(num_runs,L)];
cell_utilize_tot_PSR = [zeros(num_runs,L)];
    min_R_factor_f = []; % minimum R factor in the system in boston 
    average_R_factor_f2 = [];% average R factor for college admission 
    utilize = []; % utilization of cells in system college admission 
          utilize1 = [];
    utilize2= [];% utilization of cells in system boston admission 
    timing = [];% for college admission
    timing2=[];% for boston mechanism 
          timing1 = [];
iter_tot = [];
iter_tot_PSR = [];
    for zz = 1:num_runs
        zz
        z
        fdsa_lk = [];
        fdsa2_lk = [];
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
      % since backhaul delay is pareto distribution
        Wireless_Delay= ((qouta-1)*slot_time); % equation (10) in paper 
        maximum_delay = Wireless_Delay + BackHual_Delay;
       
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
        PSR_up_2= (1-(0.5*(erfc(sqrt((sinr_up_2)))))).^number_packets; % equation (1) in paper walid saad
        
        % construct cell preference:
        % the small cells get a separate benefit from the macro cell.
        small_cell_benfit_UP= PSR_up_2./ repmat(max(PSR_up_2(:,1:L_M)')',1,L);
        % we do this because the paper suggested that we divide by the PSR
        % of the associated macro base station, that's why I used the
        % maximum function
        cell_benfit_UP = [PSR_up_2(:,1:L_M), small_cell_benfit_UP(:,L_M+1:end)];% equation (6,7)
        %%%%%%%% R factor guarantee calculations%%%%%%%%%%%
        frew = epson3;
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
%         small_cell_benfit_R= R_factor./ repmat(max(R_factor(:,1:L_M)')',1,L);
%         R_factor_benefit = [R_factor(:,1:L_M), small_cell_benfit_R(:,L_M+1:end)];
        %preference lists for both users and cells
        [~,R_factor_index]=sort(R_factor,2,'descend');
        [~,PSR_up_2_index]=sort(R_factor,2,'descend');
        
        %%%%%%%college admission begins%%%%%%%%
        %preference lists for all users
        for i = 1:N
            user_preference.([string('user'),num2str(i)]) =R_factor_index(i,:);
                         user_preference1.([string('user'),num2str(i)]) =R_factor_index(i,:);

            user_preference_PSR.([string('user'),num2str(i)]) =PSR_up_2_index(i,:);
        end
        %dynamic list name of users per cell
        for i = 1:L
            cell_association.([string('cell'),num2str(i)]) =[];
            cell_association_rank.([string('cell'),num2str(i)]) =[];
            cell_association_PSR.([string('cell'),num2str(i)]) =[];
                                          cell_association_rank_PSR.([string('cell'),num2str(i)]) =[];
                     cell_association1.([string('cell'),num2str(i)])  = [];   
%                     cell_association1_rank.([string('cell'),num2str(i)]) = [];
             

        end
        %this flag will represent whether the user still have preference list
        has_preference = ones(1,N);
                has_preference_PSR = ones(1,N);

        %this flag will represent whether the user is associated or not
        is_NOT_associated(1,:) = ones(1,N);
        is_NOT_associated_PSR = ones(1,N);
        iter = 0;
        tstart= tic;
        %%%%%%%%%%%%% college admission %%%%%%%%%%%%
                while (sum( is_NOT_associated(1,:))>0)
                                        iter =  iter+1;
 fdsa_lk = [ fdsa_lk ,sum( is_NOT_associated(1,:)) ];
            for i = 1:N
                if ( is_NOT_associated(1,i) == 0)
                    continue;
                else
                    
                        application = user_preference.([string('user'),num2str(i)])(1);% submit application
                        rank = cell_benfit_UP(i,application);% application evaluated by cell
                        if length(cell_association.([string('cell'),num2str(application)]))+1<= qouta(1,application) % if there is a place 
                            cell_association.([string('cell'),num2str(application)]) =[cell_association.([string('cell'),num2str(application)]) i]; 
                            cell_association_rank.([string('cell'),num2str(application)]) =[cell_association_rank.([string('cell'),num2str(application)]) rank];
                            is_NOT_associated(1,i) = 0;
                        else
                            [p,j] = min(cell_association_rank.([string('cell'),num2str(application)]));
                            if rank > p
                                x = cell_association.([string('cell'),num2str(application)])(j);
                                cell_association.([string('cell'),num2str(application)])(j)=i;
                                cell_association_rank.([string('cell'),num2str(application)])(j)=rank;
                                is_NOT_associated(1,i) = 0;
                                user_preference.([string('user'),num2str(x)]) = user_preference.([string('user'),num2str(x)])(2:end);
                                is_NOT_associated(1,x) = 1;
                            else
                                user_preference.([string('user'),num2str(i)]) = user_preference.([string('user'),num2str(i)])(2:end);
                            end
                        end
                    
                end
            end
%             end
                end

%%%%%%%%%%%%%%% end of college admission%%%%%%%%%%%%%%%%
        timing = [timing toc(tstart)];
                     iter_tot = [iter_tot ,iter];

        quota_actual = ((qouta));
%          for i = 1:L % define temporary list for applicants in all users 
%                 
%                     cell_associate_boston.(([string('cell'),num2str(i)])) = [];
%                     cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];
%                 end
%         tstart2 = tic;
%%%%%%%%%% boston school choice %%%%%%%%%%%%%%%%%%%%
%   
%         for e = 1:L
%             if (sum(is_NOT_associated_PSR(1,:))==0)
%                 % terminate algorithm if we no longer have users to
%                 % associate
%                 break;
%             else
%                 
%                
%                 for i_PSR = 1:N
%                     if  (is_NOT_associated_PSR(1,i_PSR)==0)
%                         continue;
%                     else
%                         application_PSR = user_preference_PSR.([string('user'),num2str(i_PSR)])(e);% application submitted
%                         
%                         if (quota_actual(1,application_PSR))>0 % cells must have places available for this round
%                             
%                         rank_PSR = PSR_up_2(i_PSR,application_PSR); % application evaluated by cell
%                                                     cell_associate_boston.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston.([string('cell'),num2str(application_PSR)]) i_PSR];
%                             cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) rank_PSR];
%                         end                        
%                     end
%                 end
%                 %add accepted users in this round to the list of accepted
%                 %users from previous rounds and cells calculate remaining
%                 %seats for next round
%                 for i = 1:L
%                     if isempty([cell_associate_boston.(([string('cell'),num2str(i)]))])
%                         continue;
%                     elseif (length([cell_associate_boston.(([string('cell'),num2str(i)]))])<= quota_actual(i))&&~isempty([cell_associate_boston.(([string('cell'),num2str(i)]))])
%                         cell_association_PSR.([string('cell'),num2str(i)]) = [cell_association_PSR.([string('cell'),num2str(i)]), cell_associate_boston.(([string('cell'),num2str(i)]))];
%                         is_NOT_associated_PSR(1,cell_associate_boston.(([string('cell'),num2str(i)]))) = 0;
%                          quota_actual(i) = qouta(i) - length(cell_association_PSR.([string('cell'),num2str(i)]));
%                         cell_associate_boston.(([string('cell'),num2str(i)])) = [];
%                     cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];
%                     else
%                           [~,index] = sort(cell_associate_boston_rank.(([string('cell'),num2str(i)])),2,'descend');
%                           cell_association_PSR.([string('cell'),num2str(i)]) = [cell_association_PSR.([string('cell'),num2str(i)]), cell_associate_boston.(([string('cell'),num2str(i)]))(1:index(quota_actual(i)))];
%                       is_NOT_associated_PSR(1,cell_associate_boston.(([string('cell'),num2str(i)]))(1:index(quota_actual(i)))) = 0;
%                      quota_actual(i) = qouta(i) - length(cell_association_PSR.([string('cell'),num2str(i)])); 
%                       cell_associate_boston.(([string('cell'),num2str(i)])) = [];
%                     cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];
%                     end
%                     end
%                     
%                 end
%             end
%         
%         timing2 = [timing2 toc(tstart2)];
        %%%%%%%%%% boston school choice end %%%%%%%%%%%%%%%%%%%
        for i = 1:L % define temporary list for applicants in all users 
                
                    cell_associate_boston.(([string('cell'),num2str(i)])) = [];
                    cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];
        end
                iter_PSR = 0;
        tstart2 = tic;
%%%%%%%%%% boston school choice %%%%%%%%%%%%%%%%%%%%
  
        for e = 1:L
             fdsa2_lk = [ fdsa2_lk ,sum( is_NOT_associated_PSR(1,:)) ];

            if (sum(is_NOT_associated_PSR(1,:))==0)
                % terminate algorithm if we no longer have users to
                % associate
                break;
            else
                iter_PSR = iter_PSR+1;
                for i = 1:L % define temporary list for applicants in all users 
                if ~isempty(cell_associate_boston.(([string('cell'),num2str(i)])) )
                    cell_associate_boston.(([string('cell'),num2str(i)])) = [];
                    cell_associate_boston_rank.(([string('cell'),num2str(i)])) = [];
                end
                end
                for i_PSR = 1:N
                    if  (is_NOT_associated_PSR(1,i_PSR)==0)
                        continue;
                    else
                        application_PSR = user_preference_PSR.([string('user'),num2str(i_PSR)])(e);% application submitted
                         if ((quota_actual(1,application_PSR))==0) % cells must have places available for this round
                            continue;
                         else
                                               
%                         else
                        rank_PSR = cell_benfit_UP(i_PSR,application_PSR); % application evaluated by cell
                        if length(cell_associate_boston.([string('cell'),num2str(application_PSR)]))+1<= (quota_actual(1,application_PSR)) 
                            % accept user if cell has a place for it
                            cell_associate_boston.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston.([string('cell'),num2str(application_PSR)]) i_PSR];
                            cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]) rank_PSR];
                            is_NOT_associated_PSR(1,i_PSR) =0;
                        else
                            % get user with minimum rank associated in this
                            % round 
                            [p,j] = min(cell_associate_boston_rank.([string('cell'),num2str(application_PSR)]));
                            %compare it to new applicant
                            if rank_PSR > p
                                % accept the new applicant , reject the
                                % previously accepted user in this round
                                x = cell_associate_boston.([string('cell'),num2str(application_PSR)])(j);
                                cell_associate_boston.([string('cell'),num2str(application_PSR)])(j)=i_PSR;
                                cell_associate_boston_rank.([string('cell'),num2str(application_PSR)])(j)=rank_PSR;
                                is_NOT_associated_PSR(1,i_PSR) = 0;
                                is_NOT_associated_PSR(1,x) =1;
                            end
                        end
                         end
                    end
                end
                %add accepted users in this round to the list of accepted
                %users from previous rounds and cells calculate remaining
                %seats for next round
                for i = 1:L
                    if isempty([cell_associate_boston.(([string('cell'),num2str(i)]))])
                        continue;
                    else
                        cell_association_PSR.([string('cell'),num2str(i)]) = [cell_association_PSR.([string('cell'),num2str(i)]), cell_associate_boston.(([string('cell'),num2str(i)]))];
                        quota_actual(i) = qouta(i) - length(cell_association_PSR.([string('cell'),num2str(i)]));
                    end
                end
            end
        end
        timing2 = [timing2 toc(tstart2)];
         iter_tot_PSR = [iter_tot_PSR ,iter_PSR];

        %%%%%%%%%%% boston school choice end %%%%%%%%%%%%%%%%%%%
%                    is_NOT_associated1 = ones(1,N);
%         quota_actual = ((qouta));
% 
%   for i = 1:L % define temporary list for applicants in all users 
%                 
%                     cell_associate_boston1.(([string('cell'),num2str(i)])) = [];
%                     cell_associate_boston1_rank.(([string('cell'),num2str(i)])) = [];
%                 end
% %%%%%%%%%% boston school choice %%%%%%%%%%%%%%%%%%%%
%           tstart1 = tic;           
% 
%         for e = 1:L
%             if (sum(is_NOT_associated1(1,:))==0)
%                 % terminate algorithm if we no longer have users to
%                 % associate
%                 break;
%             else
%                 
% 
%                 for i_PSR = 1:N
%                     if  (is_NOT_associated1(1,i_PSR)==0)
%                         continue;
%                     else
%                         for f = e:L
%                         application_PSR = user_preference1.([string('user'),num2str(i_PSR)])(f);% application submitted
%                          if ((quota_actual(1,application_PSR))==0) % cells must have places available for this round
%                             continue;
%                          else
%                              break;
%                          end
%                         end
%                            
%                         rank_PSR = cell_benfit_UP(i_PSR,application_PSR); % application evaluated by cell
%                         if length(cell_associate_boston1.([string('cell'),num2str(application_PSR)]))+1<= (quota_actual(1,application_PSR)) 
%                             % accept user if cell has a place for it
%                             cell_associate_boston1.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston1.([string('cell'),num2str(application_PSR)]) i_PSR];
%                             cell_associate_boston1_rank.([string('cell'),num2str(application_PSR)]) =[cell_associate_boston1_rank.([string('cell'),num2str(application_PSR)]) rank_PSR];
%                             is_NOT_associated1(1,i_PSR) =0;
%                         else
%                             % get user with minimum rank associated in this
%                             % round 
%                             [p,j] = min(cell_associate_boston1_rank.([string('cell'),num2str(application_PSR)]));
%                             %compare it to new applicant
%                             if rank_PSR > p
%                                 % accept the new applicant , reject the
%                                 % previously accepted user in this round
%                                 x = cell_associate_boston1.([string('cell'),num2str(application_PSR)])(j);
%                                 cell_associate_boston1.([string('cell'),num2str(application_PSR)])(j)=i_PSR;
%                                 cell_associate_boston1_rank.([string('cell'),num2str(application_PSR)])(j)=rank_PSR;
%                                 is_NOT_associated1(1,i_PSR) = 0;
%                                 is_NOT_associated1(1,x) =1;
%                             end
%                         end
% %                         end
%                     end
%                 end
%                 %add accepted users in this round to the list of accepted
%                 %users from previous rounds and cells calculate remaining
%                 %seats for next round
%                 for i = 1:L
%                     if isempty([cell_associate_boston1.(([string('cell'),num2str(i)]))])
%                         continue;
%                     else
%                         cell_association1.([string('cell'),num2str(i)]) = [cell_association1.([string('cell'),num2str(i)]), cell_associate_boston1.(([string('cell'),num2str(i)]))];
%                         quota_actual(i) = qouta(i) - length(cell_association1.([string('cell'),num2str(i)]));
%                          cell_associate_boston1.(([string('cell'),num2str(i)])) = [];
%                     cell_associate_boston1_rank.(([string('cell'),num2str(i)])) = [];
%                     end
%                 end
%             end
%         end
%          timing1 = [timing1 toc(tstart1)];

        rewq_1=[];
                rewq_11=[];
rewq_2 = [];
rewq_2_PSR = [];

        rewq_1_PSR=[];
                  rewq_psr=[];
        minrf = [];
        count = 0;
        cell_utilize = zeros(1,L);
        for i = 1:L
            if isempty(cell_association.([string('cell'),num2str(i)]))
                continue
            else
                count = count + (1/L);
                cell_utilize(i) = length(cell_association.([string('cell'),num2str(i)]));
                %                     if i > L_M
                New_Wireless_Delay(i)= slot_time*(length(cell_association.([string('cell'),num2str(i)]))-1);
                delay(i) = New_Wireless_Delay(i) + BackHual_Delay(i);
               
                if delay(i) <= 150
                    Idd(i)=0;
                else
                    X = log10(delay(i)/mT)/log10(2);
                    Idd(i) = 25*((1+X^(6*sT))^(1/(6*sT))-3*(1+(X/3)^(6*sT))^(1/(6*sT))+2);
                end
                if delay(i)<= frew
                    Id1i = -2.468*(10^-14)*delay(i)^6;
                    Id2i = 5.062*(10^-11)*delay(i)^5;
                    Id3i = -3.903*(10^-8)*delay(i)^4;
                    Id4i = 1.344*(10^-5)*delay(i)^3;
                    Id5i = -0.001802*delay(i)^2;
                    Id6i = 0.103*delay(i)^1;
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
                rewq = [];
                for jj = 1:length(cell_association.([string('cell'),num2str(i)]))
                    rewq = [rewq  94.2-Idi(i)-12-(15*(log(1+60*(1-PSR_up_2(cell_association.([string('cell'),num2str(i)])(jj),i)))))];
                    
                    rewq_1= [rewq_1  94.2-Idi(i)-12-(15*(log(1+60*(1-PSR_up_2(cell_association.([string('cell'),num2str(i)])(jj),i)))))];
               rewq_2 = [rewq_2, ((log2(1+(sinr_up_2(cell_association.([string('cell'),num2str(i)])(jj),i)))))];
                end
                minrf = [minrf,min(rewq)];
            end
            
        end
        utilize = [utilize count];
%         minrf1 = [];
%         count1 = 0;
%                 for i = 1:L
%             if isempty(cell_association1.([string('cell'),num2str(i)]))
%                 continue
%             else
%                 count1 = count1 + (1/L);
%                 %                     if i > L_M
%                 New_Wireless_Delay1(i)= slot_time*(length(cell_association1.([string('cell'),num2str(i)]))-1);
%                 delay1(i) = New_Wireless_Delay1(i) + BackHual_Delay(i);
%                
%                 
%                 if delay1(i) <= 150
%                     Idd(i)=0;
%                 else
%                     X = log10(delay1(i)/mT)/log10(2);
%                     Idd(i) = 25*((1+X^(6*sT))^(1/(6*sT))-3*(1+(X/3)^(6*sT))^(1/(6*sT))+2);
%                 end
%                 if delay1(i)<= frew
%                     Id1i = -2.468*(10^-14)*delay1(i)^6;
%                     Id2i = 5.062*(10^-11)*delay1(i)^5;
%                     Id3i = -3.903*(10^-8)*delay1(i)^4;
%                     Id4i = 1.344*(10^-5)*delay1(i)^3;
%                     Id5i = -0.001802*delay1(i)^2;
%                     Id6i = 0.103*delay1(i)^1;
%                     Id7i = -0.1698;
%                     Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
%                     Idi(i) =  Idii(i) + Idd(i);
%                 else
%                     Id1i = -2.468*(10^-14)*(frew)^6;
%                     Id2i = 5.062*(10^-11)*(frew)^5;
%                     Id3i = -3.903*(10^-8)*(frew)^4;
%                     Id4i = 1.344*(10^-5)*(frew)^3;
%                     Id5i = -0.001802*(frew)^2;
%                     Id6i = 0.103*(frew)^1;
%                     Id7i = -0.1698;
%                     Idii(i) = Id1i+Id2i+Id3i+Id4i+Id5i+Id6i+Id7i;
%                     Idi(i) =  Idii(i) + Idd(i);
%                 end
%                 rewq1 = [];
%                 for jj = 1:length(cell_association1.([string('cell'),num2str(i)]))
%                     rewq1 = [rewq1  94.2-Idi(i)-12-(15*(log(1+60*(1-PSR_up_2(cell_association1.([string('cell'),num2str(i)])(jj),i)))))];
%                     
%                     rewq_11= [rewq_11  94.2-Idi(i)-12-(15*(log(1+60*(1-PSR_up_2(cell_association1.([string('cell'),num2str(i)])(jj),i)))))];
%                 end
%                 minrf1 = [minrf1,min(rewq1)];
%             end
%             
%         end
%         utilize1 = [utilize1 count1];

       %%%%%%%%%% measurements for boston admission%%%%%%%%%%%%%%%%% 
        minrf_t = [];
        count2 = 0;
                cell_utilize_PSR = zeros(1,L);
        for i = 1:L
            if isempty(cell_association_PSR.([string('cell'),num2str(i)]))
                continue
            else
                                rewq_psr = [];

                count2 = count2 + (1/L);
               cell_utilize_PSR(i) = length(cell_association_PSR.([string('cell'),num2str(i)]));
                New_Wireless_Delay_psr(i)= slot_time*(length(cell_association_PSR.([string('cell'),num2str(i)]))-1);
                delay_PSR(i) = New_Wireless_Delay_psr(i) + BackHual_Delay(i);
                h_psr= delay_PSR(i)>epson3;
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
                    rewq_psr = [rewq_psr  94.2-Idi(i)-12-15*(log(1+60*(1-PSR_up_2(cell_association_PSR.([string('cell'),num2str(i)])(jj),i))))];
                    rewq_1_PSR= [rewq_1_PSR  94.2-Idi(i)-12-15*(log(1+60*(1-PSR_up_2(cell_association_PSR.([string('cell'),num2str(i)])(jj),i))))];
                             rewq_2_PSR = [rewq_2_PSR, (log2(1+(sinr_up_2(cell_association_PSR.([string('cell'),num2str(i)])(jj),i))))];

                end
                minrf_t = [minrf_t  min(rewq_psr)];
            end
        end
        utilize2 = [utilize2 count2];
        
cell_utilize_tot_PSR(zz,:)=        cell_utilize_PSR;
cell_utilize_tot(zz,:)=        cell_utilize;
        
%         sum(mean(cell_utilize_tot)(1:L_M)) mean(cell_utilize_tot_PSR)
        min_R_factor_f(zz) = min((minrf_t));
        
        average_R_factor_f2(zz) = mean((rewq_1));
                average_R_factor_f21(zz) = mean((rewq_11));

        average_R_factor_f2_PSR(zz) = mean((rewq_1_PSR));
        average_R_factor_f21rate(zz) = mean((rewq_2));

        average_R_factor_f2_PSRrate(zz) = mean((rewq_2_PSR));
        average_R_factor_f2_min(zz) = min((minrf));
%                 average_R_factor_f2_min1(zz) = min((minrf1));

        average_utilize(zz)=mean(utilize);
        average_utilize2(zz)=mean(utilize2);
%                 average_utilize1(zz)=mean(utilize1);
 fdsa(zz,1:length(fdsa_lk)) = [fdsa_lk];
 fdsa2(zz,1:length(fdsa2_lk)) = [fdsa2_lk];

    end
   max_num_iter_PSR(z) =find(mean(fdsa2)==0,1);
     max_num_iter(z) =  find(mean(fdsa)==0,1);
avg_utilization_macro(z) =sum((mean(cell_utilize_tot(:,1:L_M))) )/(N);
avg_utilization_macro_PSR(z) =sum((mean(cell_utilize_tot_PSR(:,1:L_M))) )/(N);
avg_utilization_SBS(z) =mean((mean(cell_utilize_tot(:,L_M+1:L))) )/qouta(L_M+1);
avg_utilization_SBS_PSR(z) =mean((mean(cell_utilize_tot_PSR(:,L_M+1:L))) )/qouta(L_M+1);
    average_R_factor(z) = mean(min_R_factor_f);
    average_R_factor_f_1(z)=mean((average_R_factor_f2_min));
%         average_R_factor_f_11(z)=mean((average_R_factor_f2_min1));
    average_R_factor_f2_1(z) = mean(average_R_factor_f2);
        average_R_factor_f2_11(z) = mean(average_R_factor_f21);
rate_PSR(z) = mean(average_R_factor_f2_PSRrate);
rate(z) = mean(average_R_factor_f21rate);
    average_R_factor_f2_1_PSR(z) = mean(average_R_factor_f2_PSR);
    average_utilize_f(z) = mean(average_utilize);
%         average_utilize_f1(z) = mean(average_utilize1);
    average_utilize_f2(z) = mean(average_utilize2);
    avg_timing2(z) = mean(timing2);
    avg_timing(z) = mean(timing);
%         avg_timing1(z) = mean(timing1);
avg_num_iter(z) = mean(iter_tot);
        avg_num_iter_PSR(z) = mean(iter_tot_PSR);
end
toc
%% results plotted
disp('qouta')
%(1:z)+1
figure
title('worst case R factor')
plot((1:z)+1,average_R_factor,'r-')
grid on
hold on
plot((1:z)+1,average_R_factor_f_1)
%  plot((1:z)+1,average_R_factor_f_11,'k-')

figure
title('average R factor')
plot((1:z)+1,average_R_factor_f2_1)
hold on
grid on

plot((1:z)+1,average_R_factor_f2_1_PSR,'r-')
% plot((1:z)+1,average_R_factor_f2_11,'k-')

figure
title('average percentage of cells utilized')

plot((1:z)+1,average_utilize_f2,'r-')
hold on
grid on

plot((1:z)+1,average_utilize_f)
% plot((1:z)+1,average_utilize_f1,'k-')

figure
title('average excution time')

plot((1:z)+1,avg_timing2,'r-')
hold on
grid on

plot((1:z)+1,avg_timing)
figure
title('sum rate')

plot((1:z)+1,rate_PSR,'r-')
hold on
grid on

plot((1:z)+1,rate)
% plot((1:z)+1,avg_timing1,'k-')
figure
title('maximum number of iterations')

plot((1:z)+1,max_num_iter_PSR-1,'r-')
hold on
grid on

plot((1:z)+1,max_num_iter-1)%(1:z)+1
figure
title('average number of users in macro cells')

plot((1:z)+1,avg_utilization_macro_PSR,'r-')
hold on
grid on

plot((1:z)+1,avg_utilization_macro)%(1:z)+1
figure
title('average number of users in small cells')

plot((1:z)+1,avg_utilization_SBS_PSR,'r-')
hold on
grid on
plot((1:z)+1,avg_utilization_SBS)
figure
hold on
grid on
title('average number of iterations')
xlabel('quota')
ylabel('average number of iterations')
plot((1:z)+1,avg_num_iter_PSR,'r-')

plot((1:z)+1,avg_num_iter)
legend('Boston','Gale-Shapley')