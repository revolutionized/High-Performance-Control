%% Example LQR infinte horizon problem
% Consider a cart of mass m and linear viscous coefficient d, which is
% acted by a control force F. Let us focus on the dynamics for the velocity
% v. From Newton's 2nd law, it follows that:
%		m*(dv/dt) = F - dv
% If we define x = v and u = F, then the above model can be repesented as
%		(dx/dt) = Ax + Bu
% where A = -d/m and B = 1/m. For the parameter values m = 1kg and b = 2
% Ns/m, we seek the optimal control law that minimises the following cost:
%		J(x_0, u) = integral[0, inf] (5x^2 + u^2) dt
global A B T h t
A=-2;
B=1;
%Control parameters
Q=5;
R=1;
%Solve the LQR OCP
[K,S,E] = lqr(A,B,Q,R);

% Simulation - Able to compare with MCM

%Simulation - NOT ABLE TO COMPARE THIS FORM WITH MCM
SYS=ss(A-B*K,0,1,h); %the Bcl, Ccl, Dcl are irrelevant
x0=2;
[Y,T,X] = initial(SYS,x0,t);
figure(1)
plot(T',X)
xlabel('Time [s]')
ylabel('State X(t)')
grid on
figure(2)
plot(T',-K*X')
xlabel('Time [s]')
ylabel('Control u(t)')
grid on