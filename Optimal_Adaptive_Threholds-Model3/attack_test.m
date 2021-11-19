clc;
clear;

A= [0.4450 -0.0458;1.2939 0.4402];
B= [0.0550;4.5607];
C= [0 1];
D= [0];
Ts = 0.04;

Q = [ 5 0;
	 0 1/0.1*2 
	 ];


R = 1;
R = R*Ts*Ts;
Q = Q*Ts*Ts;
N = zeros(size(A,2),size(B,2));
[K] = dlqr(A,B,Q,R,N)

sys_gain = ss(A- B*K,B,C,D,Ts);
[kalmf,L,P] = kalman(sys_gain,QN,RN);


size_x = [2 1];
size_y = [1 1]; 

timeWindow = 15;
sensorAttack = 2;
actuatorAttack = 0;
safex = [-1 ; 2];
depth = 0.1;
cusum_true = true;
cusum_cost_mat = [1]; %In case Y is also a vector, then we would require to normalize it

th_arr = zeros(timeWindow,1); %Update it by running the optimal threholds function or read it from a saved file
if isfile('files/optimal_thresholds.csv')
   th_arr = readmatrix("files/optimal_thresholds.csv");
else
    [TCP_opt, optimal_delays, th_arr] = optimal_adaptive_thresholds;
    fout = sprintf('files/optimal_thresholds.csv');
   writematrix(th_arr, fout);    
end

x_plot =zeros(size_x(1),timeWindow);
close all;

%Attack-less plot
sensorAttack_zero = 0;
actuatorAttack_zero = 0;

    x_a = depth*safex;
    xhat_a = x_a;
    u_a = -K*xhat_a;

    x = depth*safex;
    xhat = x;
    u = -K*xhat;
                
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            
            for i=1:timeWindow
                % Non-attack scenario
                u_a = u + actuatorAttack_zero;
                x = A*x + B*u_a; % state updattion in plant side
                x_plot(:,i) = x; % This will be the actual state of the plant, which is not observable to the controller though
                y = C*x + sensorAttack_zero; %sensor output % sensor measurement in plant side + attacker values sensed by the controller
                r = y - C*xhat; % residue computation in controller side
                %The estimator gets the non-actuator attack value of u but
                %sensor attack value is added
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                u = - K*xhat; % control signal computation in controller side for next actuation
                th = th_arr(i); %We have an optimal value of threshold for all time steps
                
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
      plot(x_plot(1,:));
      title("No Error");
      
sensorAttack_zero = 0;
actuatorAttack_zero = 5;

    x_a = depth*safex;
    xhat_a = x_a;
    u_a = -K*xhat_a;

    x = depth*safex;
    xhat = x;
    u = -K*xhat;
    
    k_a = 4;
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
                th = th_arr(i); %We have an optimal value of threshold for all time steps   
                
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
      
for k_a = 1:timeWindow
    
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