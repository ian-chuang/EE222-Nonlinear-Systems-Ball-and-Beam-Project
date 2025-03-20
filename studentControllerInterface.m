classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        extra_dummy1 = 0;
        extra_dummy2 = 0;
        N_MHE = 5;
        history = zeros(4, 5); % [dt, u, x, th]
        Q_est = diag([5,1,5,1]);
        R_est = diag([1, 1])*50;
    end
    properties
        A_fn = @(x1,x2,x3,x4,x5) eye(4);
        B_fn = @(x1,x2,x3,x4,x5) zeros(4);
        f_fn = @(x1,x2,x3,x4,x5) [x1;x2;x3;x4];
        Q = diag([100 0 0 0]);
        R = diag([0.3]);
        x_eq = [0, 0, 0, 0, 0];
        opti;
        X_opt;
        U_opt;
        Y_opt;
        DT_opt;
        V_servo = 0.0;
    end
    methods(Access = protected)


        function V_servo = stepImpl(obj, t, p_ball, theta)
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
            %% Sample Controller: Simple Proportional Controller
            t_prev = obj.t_prev;
            u_prev = obj.V_servo;
            % Extract reference trajectory at the current timestep.
            xhat = obj.MovingWindowEstimator(t-obj.t_prev, [u_prev, p_ball, theta]);
            %disp(xhat);
            x_eq = obj.getEqPoint(t);
            [V_servo, theta_d] = obj.LQRController(xhat, x_eq);
            obj.V_servo = V_servo;
            obj.t_prev = t;
            return

            theta_d = theta; % bad hack

            % Decide desired servo angle based on simple proportional feedback.
            k_p = 3;
            theta_d = - k_p * (p_ball - p_ball_ref);

            % Make sure that the desired servo angle does not exceed the physical
            % limit. This part of code is not necessary but highly recommended
            % because it addresses the actual physical limit of the servo motor.
            theta_saturation = 56 * pi / 180;    
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            % Simple position control to control servo angle to the desired
            % position.
            k_servo = 10;
            V_servo = k_servo * (theta_d - theta);
            
            % Update class properties if necessary.
            obj.t_prev = t;
            obj.theta_d = theta_d;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)  
            % fprintf("t=%f s\n", t);
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
        function setupDynamics(obj)
            [obj.A_fn, obj.B_fn, obj.f_fn] = symbolic_dynamics();
        end
        function [V_servo, theta_d] = LQRController(obj, xhat, x_eq)
            A = obj.A_fn(x_eq(1), x_eq(2), x_eq(3), x_eq(4), x_eq(5));
            B = obj.B_fn(x_eq(1), x_eq(2), x_eq(3), x_eq(4), x_eq(5));
            
            K = lqr(A, B, obj.Q, obj.R);

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
