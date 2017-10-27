%% Example LQR infinte horizon problem
% Consider a cart of mass m and linear viscous coefficient d, which is
% acted by a control force F. Let us focus on the dynamics for the velocity
% v. From Newton's 2nd law, it follows that:
%		m*(dv/dt) = F - dv
% If we define x = v and u = F, the the above model can be repesented as
%		(dx/dt) = Ax + Bu
% where A = -d/m and B = 1/m. For the parameter values m = 1kg and b = 2
% Ns/m, we seek the optimal control law that minimises the following cost:
%		J(x_0, u) = integral[0, inf] (5x^2 + u^2) dt

%% Algorithm 1 - Initial Setup [Markov Chain Approximation (MCA) 1-D]
% (InitScenario)
alpha_discretised = -2.5:0.01:0.5; % control parameter
global h x_discretised
% L = 3; % arbitrary length
% bound_loc_west = 0;
% bound_loc_east = L;
xdiv = 0.05;
x_discretised = 0:xdiv:3;
% x_discretised = linspace(bound_loc_west, bound_loc_east, (1/h)+1); % discretised 1-D state space
a_diags = sigma(0)*ones(length(x_discretised),1);

k_func = zeros(length(x_discretised),length(alpha_discretised)); % both alpha_iteration and $\Delta t$ rely on same variables (x and \alpha_discretised)
deltat = zeros(size(k_func)); 

y_discretised = zeros(3,1); % in 1-D can have 3 different state space jumps -> west, east, or stay same

%p_hat = zeros(length(x_discretised), length(y_discretised), length(alpha_discretised));
w_hat = zeros(length(x_discretised), length(alpha_discretised));
% the next two are to split the code up to be easier to read and follow
v_possiblities = zeros(length(y_discretised), length(alpha_discretised));
v_summed = zeros(length(alpha_discretised),1);
% BIG = 1e5;
BIG = 2;

V = x_0 * ones(2, length(x_discretised)); % dynamic programming equation
alpha_min = zeros(length(x_discretised),1); % the control values which produce optimal control

V(1,1) = BIG;
V(1,end) = BIG;
% V(2,1) = BIG;
% V(2,end) = BIG;

epsilon = 10^-6;

n = 1; % For indexing the first and previous iterations
count = 0;

while (true)

	% For this example, the BC's are simple and can be calculated before looping through x	
% 	V(n+1,1) = 0;
	alpha_min(1) = 0; % TODO: Unsure of this
% 	V(n+1,end) = 0;
	alpha_min(end) = 0; % TODO: Also unsure of this
	V(n+1,:) = V(n,:);
	
	
	for x_iteration = 2:length(x_discretised)-1 % loop over x (excluding boundaries)

        x = x_discretised(x_iteration);
        y_discretised(1) = x - h;
        y_discretised(2) = x + h;
        y_discretised(3) = x;
        
        p_hat = zeros(length(y_discretised),1);
        
		for alpha_iteration = 1:length(alpha_discretised) % loop over \alpha_discretised
            
            alpha = alpha_discretised(alpha_iteration);
			den = sigma(x)^2 + h*B_func(x, alpha_discretised); %h*abs(b(x,alpha));
			deltat(x_iteration,alpha_iteration) = Delta_t_hat(den);
			k_func(x_iteration,alpha_iteration) = k_solver(x,alpha);

			for y_iteration = 1:length(y_discretised) % loop over y_discretised
                y = y_discretised(y_iteration);
				
				% Since p_hat might not sum over y to unity we set the case
				% where y == x to be the remainder of the two probabilities
				% summed together
                if y == x
                    p_hat(y_iteration) = 1 - sum(p_hat(1:2));
                else
                    p_hat(y_iteration) = p_transition_solver(x,y,alpha,den);
				end
                
				% Need to match y_iteration to the x_iteration of same
				% value
				if y_iteration == 1
					V_y_iteration = x_iteration - 1;
				elseif y_iteration == 2
					V_y_iteration = x_iteration + 1;
				else
					V_y_iteration = x_iteration;
				end
				
				v_possiblities(y_iteration,alpha_iteration) = p_hat(y_iteration) * V(n,V_y_iteration);
			end
           
			% Find Dynamic Programming Equation
			v_summed(alpha_iteration) = sum(v_possiblities(:,alpha_iteration))	...
										+ (deltat(x_iteration,alpha_iteration)	...
										* k_func(x_iteration,alpha_iteration));
									
%             figure(4)
% 			plot(v_summed)
% 			pause(0.04)
		end
		
	
		[V(n+1,x_iteration), min_index]  = min(v_summed);
		alpha_min(x_iteration) = alpha_discretised(min_index);

	end
    
	count = count + 1;
	rel_error = norm(V(n+1,:) - V(n,:),inf);
	if rel_error <= epsilon
		break;
	end

	if count > 100
		break;
	end

	n = n+1;
	
end

%% Function definitions

function result = sigma(~)
% Given in differential equation definition
result = 0; % sigma(x) = 0 for this particular example
end

function result = B_func(x, alpha_discretised) % currently unused
% Defined as max_(\alpha) { |b(x,\alpha) | }, but here we pass it the 
% whole discretised portion of alpha
b_result = b(x, alpha_discretised);
result = max(abs(b_result));
end

function result = b(x, alpha)
% Defined in differential equation definiton
global A B
result = A*x + B*alpha;
end

function result = Delta_t_hat(den)
% Defined as (h^2)/[sigma(x)^2 + h*B_func(x,alpha)]
global h
result = (h^2)/den;
end

function result = k_solver(x, alpha)
% Defined in cost function definition
result = (5*x^2) + (alpha^2);
end

function result = p_transition_solver(x, y, alpha, den)
global h
num = (sigma(x)^2) / 2;

if y > x
    bpart = max(b(x,alpha),0);
else 
    bpart = max(-b(x,alpha),0);
end
num = num + (h*bpart);

result = num/den;
if result > 1
	disp('p_hat is > 1')
end
if result < 0
	disp('p_hat is < 0')
end
end