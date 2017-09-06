% Now we look at Matlab solution with x_(0) = 2
x_0 = 2;
% h = 10^-3; % step size
LQR_1D_Script

%% Now use Eulers method to solve LQR the two
Delta_t_euler = 10^-3;
L = 3; % arbitrary length
bound_loc_west = 0;
bound_loc_east = L;
h = 10^-6;
t = linspace(bound_loc_west, bound_loc_east, (1/h)+1); % discretised 1-D state space 

x_dot = zeros(size(t)); % derivative approximation xdot = f(x,y) = Ax + Bu
x_approx = x_0 * ones(size(x_dot)); 
u = -x_approx;

for i = 2:length(t)	
	x = x_approx(i-1);
	u(i) = -k * x;
	x_dot(i) = A*x + B*u(i);
	x_approx(i) = x + h*x_dot(i);
end