clc;
clear;

A= [1.0000    0.1000;0    1.0000];
B= [0.0050;0.1000];
C= [1 0];
D= [0];

K = [16.0302    5.6622]; % LQR gain
L = [0.9902;0.9892]; % Kalman gain
th_all = [0.1, 0.2, 0.5, 1, 1.5, 2, 3]; % Threshold Values to run with

size_x = [2 1];
size_y = [1 1]; 

timeWindow = 15;
sensorAttack = 2;
actuatorAttack = 2;
safex = [25;30];
depth = 0.1;
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
for k_a = 1:timeWindow
    
    x_a = depth*safex;
    xhat_a = zeros(size(x_a));
    u_a = -K*xhat_a;

    x = depth*safex;
    xhat = zeros(size(x));
    u = -K*xhat;
            
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            
            for i=1:timeWindow
                % Non-attack scenario
                x = A*x + B*u; % state updattion in plant side
                x_plot(:,i) = x; % This will be the actual state of the plant, which is not observable to the controller though
                y = C*x; %sensor output % sensor measurement in plant side
                r = y - C*xhat; % residue computation in controller side
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                u = - K*xhat; % control signal computation in controller side
                th = th_arr(i); %We have an optimal value of threshold for all time steps
                
                % attack scenario
                if i>=k_a % attack window is start till the time it is detected, to maximize the damage
                    x_a = A*x_a + B*(u_a + actuatorAttack);
                else
                    x_a = A*x_a + B*u_a;
                end
                y_a = C*x_a; %sensor output
                if i>=k_a % attack window
                    y_a = y_a + sensorAttack;
                end
                r_a = y_a - C*xhat_a; % residue
                xhat_a = A*xhat_a + B*u_a + L*r_a;
                u_a = - K*xhat_a;
                if cusum_true 
                    %S_p = max([0, S_p + r_a]); %If residue is +ve, it is added to S_p, if it was -ve, then S_p drops to a min 0
                    %S_n = min([0, S_n + r_a]); %If residue is -ve, it is added to S_n, making it more -ve. If it is +ve, S_n increases, moving towards 0 
                      for j = 1:size_y(1)
                          S_p(j) = max(0,S_p(j) + r_a(j));
                          S_n(j) = min(0,S_n(j) + r_a(j));
                      end
                      S_p_single = cusum_cost_mat*abs(S_p);
                      S_n_single = cusum_cost_mat*abs(S_n);
                      if(max(S_p_single,S_n_single)<th)                  
                          %We also need to count the damage in this case as
                          %it will be plotted.
                            d = d + 1;
                            p = p + norm(x-x_a,inf);
                      else
                          %Attack is detected, need to break the loop.
                        break;
                      end
                else
                    if norm(r_a,inf)<th
                        if i >= k_a% you can use any norm here. But CUSUM is suggested as it is given in the paper
                           d = d + 1;
                          p = p + norm(x-x_a,inf);
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
end

