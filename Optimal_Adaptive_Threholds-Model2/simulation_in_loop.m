% CPS Group-2: Adaptive attack detector for dynamical systems
% Runs simulation and stores delay/damage values

clc;
clear;

clc;
clear;


% System: Trajectory tracking
A = [-0.313 56.7 0;-0.0139 -0.426 0;0 56.7 0];
B = [0.232;0.0203;0];
C = [0 0 1];
D = zeros(size(C,1),size(B,2));


p = 500;
Q = p*C'*C
R = 0.1;
[K] = lqr(A,B,Q,R)
sys_cl = ss(A-B*K, B, C, D);
step(0.2*sys_cl)

Ts = 0.1;
sys_cl = ss(A-B*K,B,C,D,Ts);
QN = 500;
RN = 0.01*eye(1);
[kest,L,P] = kalman(sys_cl,QN,RN);
disp("L = ");
disp(L);
th_all = [0.1,0.2, 0.5, 1, 1.5, 2, 3]; % Threshold Values to run with

for th = th_all
    safex = [25;30;35];
    depth = 0.1;
    sensorRange = 30;
    actuatorRange = 36;
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
    
    % Save Threshold Values in CSVs
    fout = sprintf('files/delay_th=%.2f.csv', th);
    writematrix(delay, fout);
    fout = sprintf('files/damage_th=%.2f.csv', th);
    writematrix(damage, fout);
end

% Extract Delay to Threshold Mapping
d_to_th = containers.Map('KeyType','double','ValueType','double');

% Store delay - th = max threshold that can lead to delay
for th = sort(th_all)
    fcsv = sprintf('files/delay_th=%.2f.csv', th);
    delay = readmatrix(fcsv);
    unique_delays = unique(delay);
    for d = unique_delays.'
        d_int = typecast(d, 'uint32');
        if d ~= 100
            d_to_th(d) = th;
        end
    end
end
all_delays = cell2mat(keys(d_to_th));

% Store Damage Values
P_all = -100;
th_to_damage  = containers.Map('KeyType', 'double', 'ValueType', 'any');
for th = sort(th_all)
    fcsv = sprintf('files/damage_th=%.2f.csv', th);
    P = readmatrix(fcsv);
    th_to_damage(th) = P;
    if P_all == -100
        P_all = P;
    else
        P_all = cat(1, P_all, P);
    end 
end
all_damages = unique(P_all).';

% Saving as Files
% 1. Delay to max-threshold that achieves that delay
% 2. Threshold to Damage Payoff Matrix
% 3. All Delay Values
% 4. All Damage Values
save('files/delay_to_thresh', 'd_to_th'); % load('files/delay_to_thresh', 'd_to_th')
save('files/thresh_to_damage', 'th_to_damage'); % load('files/thresh_to_damage', 'th_to_damage')
save('files/all_delays', 'all_delays');
save('files/all_damages', 'all_damages');