A = [0.4450 -0.0458;
      1.2939 0.4402];
B = [0.0550;
      4.5607];
C = [0 1];
D = 0;
acc = ss(A,B,C,D);
sys_ss = acc;
x0 = [10;1]; % initial state
Ts = 0.04;
ref = 0; % reference to be tracked
QN = 500;
RN = 0.01*eye(1);
[kalmf,L,P,M] = kalmd(sys_ss,QN,RN,Ts);
% [kest,L,P] = kalman(sys_cl,QN,RN);
p = 50000;
Q = 1*p*(C'*C);
R = 0.01;
[K,S,CLP] = dlqr(A,B,Q,R)
L
time = 20;
x = [0.1;1];
z = [0;0];
u = -K*z;

plot_vectorx1hat=zeros(1,time);
plot_vectorx2hat=zeros(1,time);

for i=1:time
    r = C*x-C*z;
    z = A*z + B*u + L*r;
    x = A*x + B*u;    
    u = -K*z;
    plot_vectorx1(i) = x(1);
    plot_vectorx1hat(i) = z(1);
    plot_vectorx2(i) = x(2);
    plot_vectorx2hat(i) = z(2);
end
disp('Closed loop plant response using Observer');
clf
hold on
plot(plot_vectorx1)
plot(plot_vectorx1hat)
plot(plot_vectorx2)
plot(plot_vectorx2hat)
legend({'x1','z1','x2','z2'});
grid on
hold off