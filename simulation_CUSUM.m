clc;
clear;



% System: Trajectory tracking
A= [1.0000    0.1000;0    1.0000];
B= [0.0050;0.1000];
C= [1 0];
D= [0];

K = [16.0302    5.6622]; % LQR gain
L = [0.9902;0.9892]; % Kalman gain

safex = [25;30];
depth = 0.1;
sensorRange = 30;
actuatorRange = 36;
th = 2; %threshold
timeWindow = 15;
sensorAttack = 2;
actuatorAttack = 2;

delay = zeros(timeWindow,timeWindow);
damage = zeros(timeWindow,timeWindow);
for k_a=2:timeWindow %attack start
    for k_e=1:k_a-1 %attack end
        delay(k_e,k_a) = 100; %a large value indicating attack not possible in this sub-window
    end
end

% initialize
x_a = depth*safex;
xhat_a = zeros(size(x_a));
u_a = -K*xhat_a;

x = depth*safex;
xhat = zeros(size(x));
u = -K*xhat;

size_x = size(x_a);
size_y = [1; 1]; 

cusum_true = true;
cusum_cost_mat = [1]; %In case Y is also a vector, then we would require to normalize it

for k_a=1:timeWindow % attack start
    for k_e=k_a:timeWindow % attack end
        d = 0;
        p = 0;
        S_p = zeros(size_x);
        S_n = zeros(size_x);
        
        for i=1:timeWindow
            % Non-attack scenario
            x = A*x + B*u; % state updattion in plant side
            y = C*x; %sensor output % sensor measurement in plant side
            r = y - C*xhat; % residue computation in controller side
            xhat = A*xhat + B*u + L*r; % state estimation in controller side
            u = - K*xhat; % control signal computation in controller side
            
            % attack scenario
            if i>=k_a && i<=k_e % attack window
                x_a = A*x_a + B*(u_a + actuatorAttack);
            else
                x_a = A*x_a + B*u_a;
            end
            y_a = C*x_a; %sensor output
            if i>=k_a && i<=k_e % attack window
                y_a = y_a + sensorAttack;
            end
            r_a = y_a - C*xhat_a; % residue can be positive or negative, noise added will be zero mean
            xhat_a = A*xhat_a + B*u_a + L*r_a;
            u_a = - K*xhat_a;
            
            %%CUSUM detector
            %CUSUM detects a change in the mean value of readings. 
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
                     d = d + 1;
                     p = p + norm(x-x_a,inf);
                 else
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
        if d>0
            delay(k_e,k_a) = d;
            damage(k_e,k_a) = p;
        else
            delay(k_e,k_a) = timeWindow-k_a;
        end
    end
end