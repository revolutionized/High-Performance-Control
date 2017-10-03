close all
clear all
%% Comparison of LQR and MCM (for 1D problem) via Euelers rule formed simulation
% Eulers rule, simply spoken:
%		dx/dt = f(x,alpha)
% Let xdot = dx/dt.
%		alpha__(n+1) = alpha__(n) + \Delta t
%		x_(n+1) = x_(n) + ( \Delta t * f(x_(n),alpha_(n)) )

%% LQR Problem scenario
% Manual solution:
% We have A = -2, B = 1, Q = 5, R = 1
global A B h t x_discretised
A=-2;
B=1;
L = 3; % arbitrary length
bound_loc_west = 0;
bound_loc_east = L;
h = 10^-4; % delta_t (not delta x)
t = linspace(bound_loc_west, bound_loc_east, 10^4 + 1); % discretised 1-D state space 
% Ricatti equation is
%		A'S + SA - (SB)R^(-1)(B'S) + Q = 0
% Thus we get S = 1, and the controller is
%		u = -Kx
% where K = R^(-1)(B'S) = 1
% LQR_1D_Script
% So closed-loop system i dx/dt = (A-BK)x -> dx/dt = -3x
% Solution of ODE is x(t) = x_0 * exp(-t/3)
% and optimal cost is 
%		J(x_0) = (x_0)'(S)(x_0) = (x_0)^2

% Now we look at Matlab solution with x_(0) = 2
x_0 = 2;
k = 1;

%% Use Eulers method to solve LQR
x_dot = zeros(size(t)); % derivative approximation xdot = f(x,y) = Ax + Bu
x_approx = x_0 * ones(size(x_dot)); 
u = -x_approx;
x_dot(1) = (A*x_approx(1)) + (B*u(1));

for i = 2:length(t)
	x = x_approx(i-1);
	u(i) = -k * x;
	x_dot(i) = A*x + B*u(i);
	x_approx(i) = x + (t(i) - t(i-1))*x_dot(i);
end

x_exact = x_0 * exp(-3*t);
u_exact = -x_exact;

%% Use Markov Chain Approximation to solve LQR
MCM_1D_Solver
x_mc_approx = x_0 * ones(size(x_dot));
alpha_s = size(alpha_min);

for i = 2:length(t)	
	x = x_mc_approx(i-1);
    distindex=abs(x_discretised-x);
    [dummy,ci]=min(distindex);  % find index closest to current x value
	x_dot(i) = A*x + B*alpha_min(ci);
	x_mc_approx(i) = x + (t(i) - t(i-1))*x_dot(i);
    alpha_s(i)=alpha_min(ci);
end
alpha_s(1) = -2;

%% Plot results for comparison

figure(1)
hold on
plot(t, x_exact)
plot(t, x_approx, '--o', 'MarkerIndices', 1:(1/(10*h)):length(t))
plot(t, x_mc_approx, 'k--o', 'MarkerIndices', 1:(1/(10*h)):length(t))
legend('Exact','Eulers Approximate','Markov Approximate')
xlabel('Time [s]')
ylabel('State X(t)')
grid on
hold off
figure(2)
hold on
plot(t,u_exact)
plot(t, u, '--o', 'MarkerIndices', 1:(1/(10*h)):length(t))
plot(t, alpha_s, 'k--o', 'MarkerIndices', 1:(1/(10*h)):length(t))
legend('Exact','Eulers Approximate','Markov Approximate')
xlabel('Time [s]')
ylabel('Control u(t)')
grid on
hold off

%% MCM Problem
% MCM_1D_Solver