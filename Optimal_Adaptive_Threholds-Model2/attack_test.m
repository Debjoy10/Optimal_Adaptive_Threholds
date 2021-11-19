clc;
clear;


% System: Trajectory tracking
Ac = [-0.313 56.7 0;-0.0139 -0.426 0;0 56.7 0];
Bc = [0.232;0.0203;0];
Cc = [0 0 1];
Dc = zeros(size(Cc,1),size(Bc,2));

Ts = 0.2; %Probably in seconds

sys_c = ss(Ac,Bc,Cc,Dc);
sys = c2d(sys_c,Ts);

[A,B,C,D] = dssdata(sys);

safex = [0.3;0;1.5];
depth = 0.1;

size_x = [size(A,2) 1];
size_y = [size(C,1) 1]; 


%p = 500; 
%Q = p*C'*C;
%Q = zeros(size(A));
%R = zeros(size(C,1),size(C,1)); 
% Q = [ 0.5 0 0 ;
%       0 0.1 0;
%       0 0 0.12;
%     ];
id = 3;
Q = [ 500 0 0 ;
      0 10 0;
      0 0 100;
    ];

R = 100;
R = R*Ts*Ts;
Q = Q*Ts*Ts;
[K] = lqrd(Ac,Bc,Q,R,Ts);

% sys_cl = ss(A-B*K, B, C, D);
% step(0.2*sys_cl)
% 
% sys_cl = ss(A-B*K,B,C,D,Ts);
% QN = 500;
% RN = 0.01*eye(1);
% [kest,L,P] = kalman(sys_cl,QN,RN);
QN = 1;
RN = 10;
[kalmf,L,P] = kalman(sys,QN,RN);
disp("L = ");
disp(L);
th_all = [0.1, 0.2, 0.5, 1, 1.5, 2, 3]; % Threshold Values to run with



timeWindow = 80/Ts;
sensorAttack = 0.01;
actuatorAttack = 0.02;

cusum_true = true;
cusum_cost_mat = [1]; %In case Y is also a vector, then we would require to normalize it

 th_arr = zeros(timeWindow,1); %Update it by running the optimal threholds function or read it from a saved file
 if isfile('files/optimal_thresholds.csv')
    th_arr = readmatrix("files/optimal_thresholds.csv");
 else
    [TCP_opt, optimal_delays, th_arr] = optimal_adaptive_thresholds();
     fout = sprintf('files/optimal_thresholds.csv');
    writematrix(th_arr, fout);    
end

x_plot =zeros(size_x(1),timeWindow);
close all;

%Attack-less plot
sensorAttack_zero = 0;
actuatorAttack_zero = 0;
x0 = [0.0; 0; 0.3];

%     x_a = depth*safex;
%     xhat_a = x_a;'
    x_a = x0;
    xhat_a = x0;
    u_a = -K*xhat_a;

%     x = depth*safex;
%     xhat = x;
    x = x0;
    xhat = x0;
    u = -K*xhat;
                
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            
            for i=1:timeWindow
                % Non-attack scenario
                u_a = u + actuatorAttack_zero;
                u_a_plot(:,i) = u_a;
                x = A*x + B*u_a; % state updattion in plant side
                x_plot(:,i) = x; % This will be the actual state of the plant, which is not observable to the controller though
                y = C*x + sensorAttack_zero; %sensor output % sensor measurement in plant side + attacker values sensed by the controller
                r = y - C*xhat; % residue computation in controller side
                %The estimator gets the non-actuator attack value of u but
                %sensor attack value is added
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                xhat_plot(:,i) = xhat;
                u = - K*xhat; % control signal computation in controller side for next actuation
                u_plot(:,i) = u;
%                 th = th_arr(i); %We have an optimal value of threshold for all time steps
%                 
%                 if cusum_true 
%                     %S_p = max([0, S_p + r_a]); %If residue is +ve, it is added to S_p, if it was -ve, then S_p drops to a min 0
%                     %S_n = min([0, S_n + r_a]); %If residue is -ve, it is added to S_n, making it more -ve. If it is +ve, S_n increases, moving towards 0 
%                       for j = 1:size_y(1)
%                           S_p(j) = max(0,S_p(j) + r(j));
%                           S_n(j) = min(0,S_n(j) + r(j));
%                       end
%                       S_p_single = cusum_cost_mat*abs(S_p);
%                       S_n_single = cusum_cost_mat*abs(S_n);
%                       if(max(S_p_single,S_n_single)<th)                  
%                           %We also need to count the damage in this case as
%                           %it will be plotted.
%                             d = d + 1;
%                       else
%                           %Attack is detected, need to break the loop.
%                         break;
%                       end
%                 else
%                     if norm(r,inf)<th
%                        % if i >= k_a% you can use any norm here. But CUSUM is suggested as it is given in the paper
%                            d = d + 1;
%                        % end
%                     else
%                         break;
%                     end
%                 
%                 end

            end
      figure();
      plot((1:timeWindow)*Ts,x_plot(id,:),(1:timeWindow)*Ts,xhat_plot(id,:));
      title("No Error State");
      
      figure();
      plot((1:timeWindow)*Ts,u_plot(1,:),(1:timeWindow)*Ts,u_a_plot(1,:));
      title("No Error Actuation");
      
sensorAttack_zero = 0; %0.05 is a decently large value
actuatorAttack_zero = 0.01; %0.01 is a very very large value of attack

    x_a = depth*safex;
    xhat_a = x_a;
    u_a = -K*xhat_a;

    x = depth*safex;
    xhat = x;
    u = -K*xhat;
    
    k_a = 50*1/Ts;
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            
            for i=1:timeWindow
                % Non-attack scenario
                if i>k_a
                    u_a = u + actuatorAttack_zero;
                else
                    u_a = u;
                end
                
                x = A*x + B*u_a; % state updattion in plant side
                x_plot(:,i) = x; % This will be the actual state of the plant, which is not observable to the controller though
                
                if i >k_a
                    y = C*x + sensorAttack_zero; %sensor output % sensor measurement in plant side + attacker values sensed by the controller
                else
                    y = C*x;
                end
                %Call Kalman to get value of L
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                xhat_plot(:,i) = xhat;
                
                r = y - C*xhat; % residue computation in controller side
                r_plot(:,i) = r; 
                %The estimator gets the non-actuator attack value of u but
                %sensor attack value is added
                %Call LQR to get value of K for controller gain
                u = - K*xhat; % control signal computation in controller side for next actuation
                %th = th_arr(i); %We have an optimal value of threshold for all time steps   
                
                for j = 1:size_y(1)
                    S_p(j) = max(0,S_p(j) + r(j));
                    S_n(j) = min(0,S_n(j) + r(j));
                end
                S_p_single = cusum_cost_mat*abs(S_p);
                S_n_single = cusum_cost_mat*abs(S_n);

                S_p_plot(:,i) = S_p_single;
                S_n_plot(:,i) = S_n_single;

            end
            
      figure();
      plot(1:timeWindow,x_plot(1,:),1:timeWindow,xhat_plot(1,:));
      title("No Detector state");
      legend("x plant", "x estimate");
      figure();
      plot(1:timeWindow,r_plot(1,:));
      title("No Detector Residues");
      figure();
      plot(1:timeWindow,S_p_plot(1,:),1:timeWindow,S_n_plot(1,:));
      title("No Detector CuSum");
      legend("S_p","S_n");
      
for k_a = (1:10)*1/Ts
    
    x_a = depth*safex;
    %xhat_a = zeros(size(x_a)); %This is probably wrong, because if the
    %plant was running from -infty, and no attack took place, the initial
    %estimate in the controller and the plant should match.
    xhat_a = x_a;
    u_a = -K*xhat_a;

    x = depth*safex;
    %xhat = zeros(size(x)); %Same argument as above
    xhat = x;
    u = -K*xhat;
            
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            
            r_plot = zeros(size_y(1),timeWindow);
            x_plot = zeros(size_x(1),timeWindow);
            
            S_n_plot(size_x(1),timeWindow);
            S_p_plot(size_x(1),timeWindow);
            
            for i=1:timeWindow
                
                if i>k_a
                    u_a = u + actuatorAttack;
                else
                    u_a = u;
                end
                
                x = A*x + B*u_a; % state updattion in plant side
                x_plot(:,i) = x; % This will be the actual state of the plant, which is not observable to the controller though
                
                if i >k_a
                    y = C*x + sensorAttack; %sensor output % sensor measurement in plant side + attacker values sensed by the controller
                else
                    y = C*x;
                end
                

                
                r = y - C*xhat; % residue computation in controller side
                r_plot(:,i) = r; 
                %The estimator gets the non-actuator attack value of u but
                %sensor attack value is added
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                xhat_plot(:,i) = xhat;
                
                u = - K*xhat; % control signal computation in controller side for next actuation
                th = th_arr(i); %We have an optimal value of threshold for all time steps   
                
                if cusum_true 
                    %S_p = max([0, S_p + r_a]); %If residue is +ve, it is added to S_p, if it was -ve, then S_p drops to a min 0
                    %S_n = min([0, S_n + r_a]); %If residue is -ve, it is added to S_n, making it more -ve. If it is +ve, S_n increases, moving towards 0 
                      for j = 1:size_y(1)
                          S_p(j) = max(0,S_p(j) + r(j));
                          S_n(j) = min(0,S_n(j) + r(j));
                      end
                      S_p_single = cusum_cost_mat*abs(S_p);
                      S_n_single = cusum_cost_mat*abs(S_n);
                      
                      S_p_plot(:,i) = S_p_single;
                      S_n_plot(:,i) = S_n_single;
                      
                      if(max(S_p_single,S_n_single)<th)                  
                          %We also need to count the damage in this case as
                          %it will be plotted.
                            d = d + 1;
                      else
                          %Attack is detected, need to break the loop.
                        break;
                      end
                else
                    if norm(r_a,inf)<th
                        if i >= k_a% you can use any norm here. But CUSUM is suggested as it is given in the paper
                           d = d + 1;
                        end
                    else
                        break;
                    end
                
                end

                
            end
            figure();        
            plot(x_plot(1,:));
            title("Attack start at t = " + k_a + " detected at t = " + i);
            xlabel("Time Step");
            ylabel("Trajectory of system");
            figure();
            plot(1:timeWindow,S_p_plot(1,:),1:timeWindow,S_n_plot(1,:));
            title("Cusum Errors");
            legend("Sp","Sn");
end

