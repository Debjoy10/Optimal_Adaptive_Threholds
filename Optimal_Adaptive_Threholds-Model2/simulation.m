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

for k_a=1:timeWindow % attack start
    for k_e=k_a:timeWindow % attack end
        d = 0;
        p = 0;
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
            r_a = y_a - C*xhat_a; % residue
            xhat_a = A*xhat_a + B*u_a + L*r_a;
            u_a = - K*xhat_a;
            
            if norm(r_a,inf)<th
                if i >= k_a% you can use any norm here. But CUSUM is suggested as it is given in the paper
                    d = d + 1;
                    p = p + norm(x-x_a,inf);
                end
            else
                break;
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