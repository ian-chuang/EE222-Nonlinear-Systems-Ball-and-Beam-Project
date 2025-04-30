classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        dt = 0.01;

        rg_val = 0.0254;
        L_val = 0.4255;
        g_val = 9.81;
        K_val = 1.5;
        tau_val = 0.025;

        x_hat = [0.00; 0.00; 0.00; 0.00]; % Initial state estimate for observer
        u = 0;

        % LQR gains
%         Q = [1,0,0,0;...
%             0,0,0,0;...
%             0,0,0,0;...
%             0,0,0,0];
%         R = 0.001;
        
        u_fn;
        h_fn;
        Lfh_fn;
        Lf2h_fn;
        Lf3h_fn;
        Lf4h_fn;
        K = [31.6228   34.8466   19.1995    6.1967];

        % N_MHE = 5;
        history = zeros(4, 5); % [dt, u, x, th]
        %Q_est = diag([5,1,5,1]);
        Q_est = diag([1, 1e-2, 1e-2]);
        R_est = diag([1e3, 1e3]);
        controller = 'TV-LQR';
        observer = 'EKF';
        x_eq = [0, 0, 0, 0, 0];
        xh = [0; 0; -57*pi/180; 0]
        P = diag([0.5 0.01 0.5 0.01]);
        xhat_prev = [0.0; 0; -57*pi/180; 0];
        Q = diag([1e-2,1e-2,1e-2,1e-2]);
        R = diag([3e-3,3e-3]);

    end
    properties
        % A_fn = @(x1,x2,x3,x4,x5) eye(4);
        % B_fn = @(x1,x2,x3,x4,x5) zeros(4);
        % f_fn = @(x1,x2,x3,x4,x5) [x1;x2;x3;x4];
        Q_tvlqr = diag([100 0 0 0]);
        R_tvlqr = diag([10]);
        % opti;
        % X_opt;
        % U_opt;
        % Y_opt;
        % DT_opt;
        X_prior = [-0.19; 0; 0; 0];
        % X_prior_num;
        % W_opt;
        % P_est;
        % P_est_num;
        V_servo = 0.0;
        % ekf;
        initialState = [-0.19; 0.00; 0; 0];
    end
    methods(Access = protected)
        function A = computeJacobianF(obj, x, u, dt, g, L, rg, tau, K)
            A = eye(4);
            A(1,2) = dt;
            A(2,1) = (5/7) * dt * ((rg/L)^2 * x(4)^2 * cos(x(3))^2);
            A(2,3) = dt * ((5 * g * rg)/(7 * L)) * cos(x(3)) + ...
                     (10/7) * dt * ((L/2 - x(1)) * (rg/L)^2 * x(4)^2 * cos(x(3)) * sin(x(3)));
            A(2,4) = -(10/7) * dt * ((L/2 - x(1)) * (rg/L)^2 * x(4) * cos(x(3))^2);
            A(3,4) = dt;
            A(4,4) = 1 - dt / tau;
        end
        function [V_servo, x_hat] = stepImpl(obj, t, p_ball, theta)
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.   
            t_prev = obj.t_prev;
            obj.t_prev = t;
            obj.dt = t-t_prev;
            u_prev = obj.u;
            %% Sample Controller: Simple Proportional Controller
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            V_servo = 0;
            x_hat = zeros(4);
            % Compute state estimate
            if obj.observer == "ELO"
                obj.x_hat = obj.extendedLuenbergerObserver(obj.x_hat, obj.u, [p_ball;theta;]);
            elseif obj.observer == "EKF"
                obj.x_hat = obj.EKFpredict(obj.x_hat, obj.u, [p_ball;theta;]);
            elseif obj.observer == "pleb"
                obj.x_hat = obj.plebObserver(obj.x_hat, obj.u, [p_ball; theta]);
            else
                % error("invalid observer")
                disp("invalid observer")
            end

            x_hat(1) = obj.x_hat(1);
            x_hat(2) = obj.x_hat(2);
            x_hat(3) = obj.x_hat(3);
            x_hat(4) = obj.x_hat(4);
            p_ball_obs = obj.x_hat(1);
            v_ball_obs = obj.x_hat(2);
            theta_obs = obj.x_hat(3);
            dtheta_obs = obj.x_hat(4);

            % Compute control
            % V_servo = obj.feedbackLinearizationController(p_ball, v_ball, theta, dtheta, ...
            %     p_ball_ref, v_ball_ref, a_ball_ref);
            if obj.controller == "FBL"
                V_servo = obj.feedbackLinearizationController(p_ball_obs, v_ball_obs, theta_obs, dtheta_obs, ...
                    p_ball_ref, v_ball_ref, a_ball_ref) + 1+0.7*cos(theta);
            elseif obj.controller == "TV-LQR"
                x_eq = obj.getEqPoint(t);
                [V_servo, theta_d] = obj.LQRController(x_hat, x_eq);
            else
                % error("invalid controller")
                disp("invalid controller")
            end
            % theta_saturation = 56 * pi / 180;
            % V_servo = clip(V_servo, -10, 10);
            V_servo = min(max(V_servo, -2), 3);

            
            % disp(V_servo);
            % disp([p_ball, v_ball, theta, dtheta]);
            % disp([p_ball_obs, v_ball_obs, theta_obs, dtheta_obs]);
            % disp([V_servo]);
            % disp("-----");
            
            % % (OPTIONAL) Define safe limits for theta
            % theta_min = -3*pi/8;  % Example lower bound
            % theta_max = 3*pi/8;   % Example upper bound
            % gain = 10;  % Scaling factor for control
            % if theta < theta_min
            %     V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction
            % elseif theta > theta_max
            %     V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction
            % end

            obj.u = V_servo;
            

            % % (DEFAULT) Decide desired servo angle based on simple proportional feedback.
            % k_p = 3;
            % theta_d = - k_p * (p_ball - p_ball_ref);
            % 
            % % Make sure that the desired servo angle does not exceed the physical
            % % limit. This part of code is not necessary but highly recommended
            % % because it addresses the actual physical limit of the servo motor.
            % theta_saturation = 56 * pi / 180;    
            % theta_d = min(theta_d, theta_saturation);
            % theta_d = max(theta_d, -theta_saturation);
            % 
            % % Simple position control to control servo angle to the desired
            % % position.
            % k_servo = 10;
            % V_servo = k_servo * (theta_d - theta);
            
            % % Update class properties if necessary.
            % obj.t_prev = t;
            % obj.theta_d = theta_d;
        end
    end
    
    methods(Access = public)
        function x_hat_next = plebObserver(obj, x_hat_cur, u_cur, y_next)
            % just takes derivatives with 1st order finite difference and
            % then uses an exponential filter to smooth everything
            
            a1 = exp(-0.2/obj.dt); % set time constant to 0.1s
            a2 = exp(-0.3/obj.dt); % time constant for derivatives is 0.2s
            x_hat_next = x_hat_cur;
            x_hat_next(1) = x_hat_cur(1)*a1 + y_next(1)*(1-a1);
            x_hat_next(2) = x_hat_cur(2)*a2 + ((y_next(1)-x_hat_cur(1))/obj.dt)*(1-a2);
            x_hat_next(3) = x_hat_cur(3)*a1 + y_next(2)*(1-a1);
            x_hat_next(4) = x_hat_cur(4)*a2 + ((y_next(2)-x_hat_cur(3))/obj.dt)*(1-a2);
        end
        function x_hat_next = extendedLuenbergerObserver(obj, x_hat_curr, u_curr, y_next)
            hat_x = x_hat_curr;
            hat_u = u_curr;
        
            A_eval = A_func(hat_x, hat_u);
            C_eval = C_func(hat_x);
        
            % Co = ctrb(A_eval', C_eval');
            % rank_Co = rank(Co);
            % disp(['Rank of Controllability Matrix: ', num2str(rank_Co)]);
        
            poles = [-4, -4.5, -5, -35.5];
            % if rank_Co < 4
            %     poles = [-4, -4.5, -5, -5.5];
            % else
            %     poles = [-1, -1.5, -2, -2.5];
            % end
        
%             L = place_fn(A_eval', C_eval', poles)';
            if hat_x(4)==0.0
                hat_x(4) = eps();
            end
            L = luenberger_gains_func(hat_x, hat_u);
            % disp(L)
            disp(hat_x)
            hat_y = h_fn(hat_x, 0);
            y = y_next;
            % disp(hat_y)
            
            dot_hat_x = f_fn(hat_x, hat_u) + L * (y - hat_y);
            x_hat_next = hat_x + obj.dt * dot_hat_x;
        end

        function u_next = feedbackLinearizationController(obj, p_ball, v_ball, theta, dtheta, p_ball_ref, v_ball_ref, a_ball_ref)
            z1 = h_fn([p_ball, v_ball, theta, dtheta]', 0) - p_ball_ref;
            z2 = Lfh_fn(p_ball, v_ball, theta, dtheta) - v_ball_ref;
            z3 = Lf2h_fn(p_ball, v_ball, theta, dtheta) - a_ball_ref;
            z4 = Lf3h_fn(p_ball, v_ball, theta, dtheta) - 0;
            z = [z1; z2; z3; z4];
            v = - obj.K * z;
            u_next = obj.u_fn(p_ball, v_ball, theta, dtheta, v);
        end

        function obj = studentControllerInterface()  
            % obj.controller = 'TV-LQR';
            % obj.observer = 'ELO';
            %%%%%% EKF %%%%%%%%
            % discrete_f = @(x, u) x + obj.dt * [ x(2);...
            % ((5 * obj.g_val * obj.rg_val)/(7 * obj.L_val)) * sin(x(3)) - (5/7) * ((obj.L_val/2) - x(1)) * ((obj.rg_val/obj.L_val) * x(4))^2 * cos(x(3))^2;...
            % x(4); -1 * x(4)/obj.tau_val + (obj.K_val/obj.tau_val) * u ];
            % 
            % syms x1 x2 x3 x4 u_sym dt_sym real
            % 
            % x_sym = [x1; x2; x3; x4];
            % 
            % f1_sym = x1 + dt_sym * x2;
            % f2_sym = x2 + dt_sym * (((5 * obj.g_val * obj.rg_val)/(7 * obj.L_val)) * sin(x3) - (5/7) * ((obj.L_val/2) - x1) * ((obj.rg_val/obj.L_val) * x4)^2 * cos(x3)^2);
            % f3_sym = x3 + dt_sym * x4;
            % f4_sym = x4 + dt_sym * (-1 * x4/obj.tau_val + (obj.K_val/obj.tau_val) * u_sym);
            % f_sym = [f1_sym; f2_sym; f3_sym; f4_sym];
            % 
            % J_sym = jacobian(f_sym, x_sym);
            % 
            % J_handle = matlabFunction(J_sym, 'Vars', {[x1; x2; x3; x4], u_sym, dt_sym});
            % 
            % f_Jacobian = @(x, u) J_handle(x, u, obj.dt);
            % 
            % measurement_z = @(x) [x(1); x(3)];
            % 
            % z_Jacobian = @(x) [1 0 0 0; 0 0 1 0];
            % 
            % obj.ekf = extendedKalmanFilter(discrete_f, measurement_z, obj.initialState, ...
            %     'StateTransitionJacobianFcn', f_Jacobian, ...
            %     'MeasurementJacobianFcn', z_Jacobian);
            % 
            % obj.ekf.ProcessNoise = diag([1e-3, 1e-2, 1e-3, 1e-3]);
            % obj.ekf.MeasurementNoise = diag([1e-6, 1e-6]);

            %%%%%% MHE %%%%%%%%
            % obj.setupDynamics();
            % obj.setupMHE();

        end

        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d, x_hat] = stepController(obj, t, p_ball, theta)        
            [V_servo, x_hat] = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
        function setupDynamics(obj)
            [obj.A_fn, obj.B_fn, obj.f_fn] = symbolic_dynamics();
        end
        function [V_servo, theta_d] = LQRController(obj, xhat, x_eq)
            A = A_func([x_eq(1), x_eq(2), x_eq(3), x_eq(4)]', x_eq(5));
            B = B_func([x_eq(1), x_eq(2), x_eq(3), x_eq(4)]', x_eq(5));
            
            K = lqr_fn(A, B, obj.Q_tvlqr, obj.R_tvlqr);

            V_servo = -K*[xhat(1)-x_eq(1); xhat(2)-x_eq(2); xhat(3)-x_eq(3); xhat(4)-x_eq(4)] + x_eq(5);
            theta_d = xhat(3);
        end
        function x_eq = getEqPoint(obj, t)
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            opts = optimoptions("fsolve", "Algorithm", "levenberg-marquardt", "Display", "none");
            x_eq = fsolve(@(x) [1, 1, 0, 0]*(f_fn([p_ball_ref, v_ball_ref, x(1), x(2)]', x(3)) ...
                - [v_ball_ref; a_ball_ref; 0.0; 0.0]), ...
                [obj.x_eq(3), obj.x_eq(4), obj.x_eq(5)], ...
                opts);
            x_eq = [p_ball_ref, v_ball_ref, x_eq(1), x_eq(2), x_eq(3)];
            obj.x_eq = x_eq;
        end
        % function setupMHE(obj)
        %     import casadi.*
        % 
        %     g = 9.81;
        %     r_arm = 0.0254;
        %     L = 0.4255;
        % 
        %     a = 5 * g * r_arm / (7 * L);
        %     b = (5 * L / 14) * (r_arm / L)^2;
        %     c = (5 / 7) * (r_arm / L)^2;
        % 
        %     K = 1.5;
        %     tau = 0.025;
        % 
        %     vars = MX.sym('x', 6);
        %     x = MX.sym('x', 4);
        %     u = MX.sym('u');
        %     dt = MX.sym('dt');
        %     f = Function('f', {x, u}, {[
        %         x(2)
        %         u * sin(x(3)) - b * x(4)^2 * cos(x(3))^2 + c * x(1) * x(4)^2 * cos(x(3))^2
        %         x(4)
        %         (- x(4) + K * u) / tau
        %     ]}, {'x', 'u'}, {'xdot'});
        % 
        %     % rk4 integration for discretization
        %     k1 = dt*f(x, u);
        %     k2 = dt*f(x+dt*k1/2, u);
        %     k3 = dt*f(x+dt*k2/2, u);
        %     k4 = dt*f(x+dt*k3, u);
        %     xf = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        % 
        %     xf = x + dt*f(x+(dt/2)*f(x, u), u); % also works if we need it
        %     %a bit faster. doesn't really hurt performance tbh
        %     F = Function('F', {x, u, dt}, {xf}).map(obj.N_MHE-1);
        % 
        % 
        %     opti = Opti();
        %     % X0 = opti.variable(4);
        %     % U0 = opti.variable(1);
        %     vars = opti.variable(7, obj.N_MHE);
        %     X = vars(1:4, :);
        %     W = vars(5:7, 1:(end-1));
        %     %X = opti.variable(4, obj.N_MHE);
        %     %W = opti.variable(obj.N_MHE-1);
        %     U = opti.parameter(1, obj.N_MHE-1);
        %     Y = opti.parameter(2, obj.N_MHE);
        %     DT = opti.parameter(obj.N_MHE-1);
        %     X_prior = opti.parameter(4);
        %     P_est = opti.parameter(4, 4);
        %     disp(size(W(2, :)));
        %     disp(obj.N_MHE);
        %     perturbation = [zeros(1, obj.N_MHE-1); W(2, :); zeros(1, obj.N_MHE-1); W(3, :)];
        %     disp(perturbation);
        %     dynamics_gap = X(:, 2:end) - F(X(:,1:end-1), U + W(1, :), DT) + perturbation;
        %     observation_gap = X([1,3], :) - Y;
        %     %dummy_param = opti.parameter(1);
        %     %opti.set_value(dummy_param, 0);
        % 
        %     % disp(size(dynamics_gap));
        %     % disp(size(observation_gap));
        %     theta_max = deg2rad(60); % physical limit
        %     x_max = 0.20; % physical limit
        % 
        %     cost = 0.0;
        %     %opti.subject_to( X0==[0;0;0;0])
        %     %opti.subject_to(X(:, 1) - X0*dummy_param + U0*dummy_param == 0)
        %     for i=1:(obj.N_MHE)
        %         if(i<obj.N_MHE)
        %             cost = cost + bilin(obj.Q_est, W(:, i))*DT(i) + bilin(obj.R_est, observation_gap(:, i))*DT(i);
        %             opti.subject_to(dynamics_gap(:, i) == [0; 0; 0; 0]);
        %         end
        %         opti.subject_to(-x_max <= X(1, i));
        %         opti.subject_to(          X(1, i) <= x_max);
        %         opti.subject_to(-theta_max <= X(3, i));
        %         opti.subject_to(              X(3, i) <= theta_max);
        %     end
        %     cost = cost + bilin(obj.R_est, observation_gap(:, obj.N_MHE))*0.01;
        %     cost = cost + bilin(P_est, X(:, 1) - X_prior)*DT(1);
        %     cost = cost + sum(vars(5:7, end) .^ 2);
        %     opti.minimize(cost);
        %     opts = struct;
        %     opts.expand = true;
        %     opts.ipopt.linear_solver = 'mumps'; % default; comes preinstalled. Small problem so mumps is good
        %     opts.ipopt.print_level = 0;
        %     opts.print_time = 0;
        %     opts.ipopt.max_wall_time = 0.010; % 10ms is super safe - typ. is 2-3ms
        %     %opts = struct;
        %     %opts.structure_detection = 'auto';
        %     %opts.debug = true;
        %     opti.solver('ipopt', opts);
        % 
        %     obj.opti = opti;
        %     obj.X_opt = X;
        %     obj.U_opt = U;
        %     obj.Y_opt = Y;
        %     obj.DT_opt = DT;
        %     obj.W_opt = W;
        %     obj.X_prior = X_prior;
        %     obj.P_est = P_est;
        %     obj.X_prior_num = [0,0,0,0];
        %     obj.P_est_num = zeros(4, 4);
        % 
        % 
        % end
        % function xhat = MovingWindowEstimator(obj, dt, y)
        %     obj.history = [obj.history(:, 2:end), [0; 0; y(2); y(3)]]; % zero is unused.
        %     obj.history(1, end-1) = dt; % technically dt is dt_prev
        %     obj.history(2, end-1) = y(1); % y(1) is u_prev
        %     %disp(obj.history);
        %     obj.opti.set_value(obj.X_prior, obj.X_prior_num);
        %     obj.opti.set_value(obj.P_est, obj.P_est_num);
        % 
        %     obj.opti.set_value(obj.U_opt, clip(obj.history(2, 1:end-1), -10, 10));
        %     obj.opti.set_value(obj.Y_opt, obj.history(3:4, :));
        %     obj.opti.set_value(obj.DT_opt, obj.history(1, 1:end-1));
        %     sol = obj.opti.solve_limited(); % solve_limited makes it not error if it hits time or iter limits
        %     Xhat = sol.value(obj.X_opt);
        %     What = sol.value(obj.W_opt);
        %     obj.opti.set_initial(obj.X_opt, Xhat); % set up warmstarting
        %     obj.opti.set_initial(obj.W_opt, What);
        %     obj.X_prior_num = Xhat(:, 2);
        %     xprior_cell = num2cell(obj.X_prior_num);
        %     Uhat = What(1, 2) + clip(obj.history(2, 2), -10, 10);
        %     G = [obj.B_fn(xprior_cell{:}, Uhat), [0;1;0;0], [0;0;0;1]]*obj.history(1, 1);
        %     A = obj.A_fn(xprior_cell{:}, Uhat)*obj.history(1, 1) + eye(4);
        %     P = obj.P_est_num;
        %     C = [1 0 0 0; 0 0 1 0];
        %     obj.P_est_num = G*obj.Q_est*G' + A*P*A' - A*P*C'*inv(obj.R_est + C*P*C')*C*P*A';
        %     % disp(obj.P_est_num)
        %     xhat = Xhat(:, end);
        %     % disp(Xhat);
        % end
        function x_hat = EKFpredict(obj, x_cur, u_cur, y_next)
            p_ball = y_next(1);
            theta = y_next(2);
            u_prev = obj.u;
            x_prev = obj.xhat_prev;
    
            f = @(x,u) x + obj.dt*[
                x(2);
                (5*obj.g_val*obj.rg_val/(7*obj.L_val))*sin(x(3)) ...
                  - (5/7)*((obj.L_val/2)-x(1))*((obj.rg_val/obj.L_val)*x(4))^2*cos(x(3))^2;
                x(4);
                -x(4)/obj.tau_val + (obj.K_val/obj.tau_val)*u
            ];
    
            x_pred = f(x_prev, u_prev);
            A = obj.computeJacobianF(x_prev, u_prev, obj.dt, ...
                                     obj.g_val, obj.L_val, obj.rg_val, obj.tau_val, obj.K_val);
    
    
            obj.P = A * obj.P * A' + obj.Q;
    
            y = [p_ball; theta];
            h = @(x) [x(1); x(3)];
            C = [1 0 0 0;
                0 0 1 0];
            y_pred = h(x_pred);
            S = C*obj.P*C' + obj.R ;
            K_gain = obj.P*C'/S;
            obj.xh = x_pred + K_gain*(y - y_pred);
            obj.P = (eye(4) - K_gain*C)*obj.P;
            obj.xhat_prev = obj.xh;

            x_hat = obj.xh;

        end
        % function x_predictEKF = EKFpredict(obj, dummy, u_curr, y_next)
        %     zTrue = y_next;
        %     z = zTrue + chol(obj.ekf.MeasurementNoise)' * randn(2,1);
        %     predict(obj.ekf, u_curr);
        %     ekfState = correct(obj.ekf, z);
        %     x_predictEKF = ekfState;
        % end

    end
    
end
