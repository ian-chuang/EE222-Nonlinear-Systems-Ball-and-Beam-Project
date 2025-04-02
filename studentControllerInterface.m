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

        x_hat = [-0.19; 0.00; 0; 0]; % Initial state estimate for observer
        u = 0;

        % LQR gains
        Q = [1,0,0,0;...
            0,0,0,0;...
            0,0,0,0;...
            0,0,0,0];
        R = 0.001;
        
        u_fn;
        h_fn;
        Lfh_fn;
        Lf2h_fn;
        Lf3h_fn;
        Lf4h_fn;
        K;

        f;
        h;
        A_func;
        C_func;

        N_MHE = 5;
        history = zeros(4, 5); % [dt, u, x, th]
        Q_est = diag([5,1,5,1]);
        R_est = diag([1, 1])*50;
    end
    properties
        A_fn = @(x1,x2,x3,x4,x5) eye(4);
        B_fn = @(x1,x2,x3,x4,x5) zeros(4);
        f_fn = @(x1,x2,x3,x4,x5) [x1;x2;x3;x4];
        Q_tvlqr = diag([100 0 0 0]);
        R_tvlqr = diag([0.3]);
        x_eq = [0, 0, 0, 0, 0];
        opti;
        X_opt;
        U_opt;
        Y_opt;
        DT_opt;
        V_servo = 0.0;
        controller;
        observer;
    end
    methods(Access = protected)
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
            u_prev = obj.u;
            %% Sample Controller: Simple Proportional Controller
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            % Compute state estimate
            if obj.observer == "ELO"
                obj.x_hat = obj.extendedLuenbergerObserver(obj.x_hat, obj.u, [p_ball;theta;]);
            elseif obj.observer == "MHE"
                obj.x_hat = obj.MovingWindowEstimator(t-obj.t_prev, [u_prev, p_ball, theta]);
            else
                error("invalid observer")
            end

            x_hat = obj.x_hat;
            p_ball_obs = obj.x_hat(1);
            v_ball_obs = obj.x_hat(2);
            theta_obs = obj.x_hat(3);
            dtheta_obs = obj.x_hat(4);

            % Compute control
            % V_servo = obj.feedbackLinearizationController(p_ball, v_ball, theta, dtheta, ...
            %     p_ball_ref, v_ball_ref, a_ball_ref);
            if obj.controller == "FBL"
                V_servo = obj.feedbackLinearizationController(p_ball_obs, v_ball_obs, theta_obs, dtheta_obs, ...
                    p_ball_ref, v_ball_ref, a_ball_ref);
            elseif obj.controller == "TV-LQR"
                x_eq = obj.getEqPoint(t);
                [V_servo, theta_d] = obj.LQRController(x_hat, x_eq);
            else
                error("invalid controller")
            end
            theta_saturation = 56 * pi / 180;
            V_servo = clip(V_servo, -theta_saturation, theta_saturation);

            
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
            obj.t_prev = t;

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

        function x_hat_next = extendedLuenbergerObserver(obj, x_hat_curr, u_curr, y_next)
            hat_x = x_hat_curr;
            hat_u = u_curr;
        
            A_eval = obj.A_func(hat_x, hat_u);
            C_eval = obj.C_func(hat_x);
        
            % Co = ctrb(A_eval', C_eval');
            % rank_Co = rank(Co);
            % disp(['Rank of Controllability Matrix: ', num2str(rank_Co)]);
        
            poles = [-4, -4.5, -5, -35.5];
            % if rank_Co < 4
            %     poles = [-4, -4.5, -5, -5.5];
            % else
            %     poles = [-1, -1.5, -2, -2.5];
            % end
        
            L = place(A_eval', C_eval', poles)';
            hat_y = obj.h(hat_x);
            
            y = y_next;
            
            dot_hat_x = obj.f(hat_x, hat_u) + L * (y - hat_y);
            x_hat_next = hat_x + obj.dt * dot_hat_x;
        end

        function u_next = feedbackLinearizationController(obj, p_ball, v_ball, theta, dtheta, p_ball_ref, v_ball_ref, a_ball_ref)
            z1 = obj.h_fn(p_ball, v_ball, theta, dtheta) - p_ball_ref;
            z2 = obj.Lfh_fn(p_ball, v_ball, theta, dtheta) - v_ball_ref;
            z3 = obj.Lf2h_fn(p_ball, v_ball, theta, dtheta) - a_ball_ref;
            z4 = obj.Lf3h_fn(p_ball, v_ball, theta, dtheta) - 0;
            z = [z1; z2; z3; z4];
            v = - obj.K * z;
            u_next = obj.u_fn(p_ball, v_ball, theta, dtheta, v);
        end

        function obj = studentControllerInterface(controller, observer)
            obj.controller = controller;
            obj.observer = observer;


            syms x1 x2 x3 x4 u g rg L K tau v real
            f = [ 
                x2;
                (5*g*rg/(7*L)) * sin(x3) - (5/7) * ((L/2) - x1) * ((rg/L)^2) * x4^2 * cos(x3)^2;
                x4;
                -x4/tau
            ];
            g_vec = [0; 0; 0; K/tau];

            
            h = x1; 
            Lfh = simplify(jacobian(h, [x1; x2; x3; x4]) * (f+g_vec*u));
            Lf2h = simplify(jacobian(Lfh, [x1; x2; x3; x4]) * (f+g_vec*u));
            Lf3h = simplify(jacobian(Lf2h, [x1; x2; x3; x4]) * (f+g_vec*u));
            Lf3h_discard = simplify(expand(Lf3h) - (-(5*K*rg^2*u*x4*cos(x3)^2)/(7*L*tau) + (10*K*rg^2*u*x1*x4*cos(x3)^2)/(7*L^2*tau)));
            Lf4h = jacobian(Lf3h_discard, [x1; x2; x3; x4]) * (f+g_vec*u);
            u_sol = simplify(solve(Lf4h == v, u));

            symbols = {g, rg, L, K, tau};
            consts = {obj.g_val, obj.rg_val, obj.L_val, obj.K_val, obj.tau_val};
            u_sol = subs(u_sol, symbols, consts);
            h = subs(h, symbols, consts);
            Lfh = subs(Lfh, symbols, consts);
            Lf2h = subs(Lf2h, symbols, consts);
            Lf3h_discard = subs(Lf3h_discard, symbols, consts);
                        
            obj.u_fn = matlabFunction(u_sol, 'Vars', [x1, x2, x3, x4, v]);
            obj.h_fn = matlabFunction(h, 'Vars', [x1, x2, x3, x4]);
            obj.Lfh_fn = matlabFunction(Lfh, 'Vars', [x1, x2, x3, x4]);
            obj.Lf2h_fn = matlabFunction(Lf2h, 'Vars', [x1, x2, x3, x4]);
            obj.Lf3h_fn = matlabFunction(Lf3h_discard, 'Vars', [x1, x2, x3, x4]);

            A = zeros(4,4);  % Initialize a 4x4 zero matrix
            A(1:3,2:4) = eye(3);  % Set the top-right 3x3 block to identity
            B = [0; 0; 0; 1];
            [K,S,P] = lqr(A,B,obj.Q,obj.R);
            obj.K = K;

            %%%%%%%%%% Extended Luenberger Observer

            rg = obj.rg_val;
            L = obj.L_val;
            g = obj.g_val;
            K = obj.K_val;
            tau = obj.tau_val;

            obj.f = @(x, u) [ 
                x(2); 
                ((5 * g * rg)/(7 * L)) * sin(x(3)) - (5/7) * ((L/2) - x(1)) * ((rg/L) * x(4))^2 * cos(x(3))^2;
                x(4); 
                -1 * x(4)/tau + (K/tau) * u
                ];

            obj.h = @(x) [
                x(1); 
                x(3)
                ]; % ball position and servo angle

            syms x1 x2 x3 x4 u_sym real

            x_sym = [x1; x2; x3; x4];
            f_sym_continuous = subs(obj.f(x_sym, u_sym));

            A_sym = jacobian(f_sym_continuous, x_sym);
            obj.A_func = matlabFunction(A_sym, 'Vars', {x_sym, u_sym});

            obj.C_func = @(x) [1 0 0 0; 0 0 1 0];


            obj.setupDynamics();
            obj.setupMHE();

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
            A = obj.A_fn(x_eq(1), x_eq(2), x_eq(3), x_eq(4), x_eq(5));
            B = obj.B_fn(x_eq(1), x_eq(2), x_eq(3), x_eq(4), x_eq(5));
            
            K = lqr(A, B, obj.Q_tvlqr, obj.R_tvlqr);

            V_servo = -K*[xhat(1)-x_eq(1); xhat(2)-x_eq(2); xhat(3)-x_eq(3); xhat(4)-x_eq(4)] + x_eq(5);
            theta_d = xhat(3);
        end
        function x_eq = getEqPoint(obj, t)
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            opts = optimoptions("fsolve", "Algorithm", "levenberg-marquardt", "OptimalityTolerance", 1e-4, "Display", "none");
            x_eq = fsolve(@(x) [1, 1, 0, 0]*(obj.f_fn(p_ball_ref, v_ball_ref, x(1), x(2), x(3)) ...
                - [v_ball_ref; a_ball_ref; 0.0; 0.0]), ...
                [obj.x_eq(3), obj.x_eq(4), obj.x_eq(5)], ...
                opts);
            x_eq = [p_ball_ref, v_ball_ref, x_eq(1), x_eq(2), x_eq(3)];
            obj.x_eq = x_eq;
        end
        function setupMHE(obj)
            import casadi.*

            g = 9.81;
            r_arm = 0.0254;
            L = 0.4255;
            
            a = 5 * g * r_arm / (7 * L);
            b = (5 * L / 14) * (r_arm / L)^2;
            c = (5 / 7) * (r_arm / L)^2;
            
            K = 1.5;
            tau = 0.025;

            vars = MX.sym('x', 6);
            x = MX.sym('x', 4);
            u = MX.sym('u');
            dt = MX.sym('dt');
            f = Function('f', {x, u}, {[
                x(2)
                u * sin(x(3)) - b * x(4)^2 * cos(x(3))^2 + c * x(1) * x(4)^2 * cos(x(3))^2
                x(4)
                (- x(4) + K * u) / tau
            ]}, {'x', 'u'}, {'xdot'});

            % rk4 integration for discretization
            k1 = dt*f(x, u);
            k2 = dt*f(x+dt*k1/2, u);
            k3 = dt*f(x+dt*k2/2, u);
            k4 = dt*f(x+dt*k3, u);
            xf = x + (k1 + 2*k2 + 2*k3 + k4)/6;

            %xf = x + dt*f(x+(dt/2)*f(x, u), u); % also works if we need it
            %a bit faster. doesn't really hurt performance tbh
            F = Function('F', {x, u, dt}, {xf}).map(obj.N_MHE-1);
            

            opti = Opti();
            X = opti.variable(4, obj.N_MHE);
            U = opti.parameter(obj.N_MHE-1);
            Y = opti.parameter(2, obj.N_MHE);
            DT = opti.parameter(obj.N_MHE-1);
            dynamics_gap = F(X(:,1:end-1), U, DT) - X(:, 2:end);
            observation_gap = X([1,3], :) - Y;

            % disp(size(dynamics_gap));
            % disp(size(observation_gap));

            cost = bilin(obj.R_est, observation_gap(:, obj.N_MHE))*0.01; % tune number here to be good guess for runtime frequency
            for i=1:(obj.N_MHE-1)
                cost = cost + bilin(obj.Q_est, dynamics_gap(:, i))*DT(i) + bilin(obj.R_est, observation_gap(:, i))*DT(i);
            end
            opti.minimize(cost);
            opts = struct;
            opts.ipopt.linear_solver = 'mumps'; % default; comes preinstalled. Small problem so mumps is good
            opts.ipopt.print_level = 0;
            opts.print_time = 0;
            opts.ipopt.max_wall_time = 0.010; % 10ms is super safe - typ. is 2-3ms
            opti.solver('ipopt', opts);

            obj.opti = opti;
            obj.X_opt = X;
            obj.U_opt = U;
            obj.Y_opt = Y;
            obj.DT_opt = DT;
            

        end
        function xhat = MovingWindowEstimator(obj, dt, y)
            obj.history = [obj.history(:, 2:end), [0; 0; y(2); y(3)]]; % zero is unused.
            obj.history(1, end-1) = dt; % technically dt is dt_prev
            obj.history(2, end-1) = y(1); % y(1) is u_prev
            %disp(obj.history);
            
            obj.opti.set_value(obj.U_opt, clip(obj.history(2, 1:end-1), -10, 10));
            obj.opti.set_value(obj.Y_opt, obj.history(3:4, :));
            obj.opti.set_value(obj.DT_opt, obj.history(1, 1:end-1));
            sol = obj.opti.solve_limited(); % solve_limited makes it not error if it hits time or iter limits
            Xhat = sol.value(obj.X_opt);
            obj.opti.set_initial(obj.X_opt, Xhat); % set up warmstarting

            xhat = Xhat(:, end);
        end
    end
    
end
