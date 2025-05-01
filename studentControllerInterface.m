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

        Q_est = diag([1, 1e-2, 1e-2]);
        R_est = diag([1e3, 1e3]);
        x_eq = [0, 0, 0, 0, 0];
        xh = [0; 0; -57*pi/180; 0]
        P = diag([0.5 0.01 0.5 0.01]);
        xhat_prev = [0.0; 0; -57*pi/180; 0];
        Q = diag([1e-2,1e-2,1e-2,1e-2]);
        R = diag([3e-3,3e-3]);

    end
    properties
        Q_tvlqr = diag([100 0 0 0]);
        R_tvlqr = diag([10]);
        X_prior = [-0.19; 0; 0; 0];
        V_servo = 0.0;
        % ekf;
        initialState = [0.0; 0.00; 0; 0];
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
            obj.t_prev = t;
            obj.dt = t-t_prev;
            u_prev = obj.u;
            %% Sample Controller: Simple Proportional Controller
            % Extract reference trajectory at the current timestep.
            x_hat = zeros(4);
            obj.x_hat = obj.EKFpredict(obj.x_hat, obj.u, [p_ball;theta;]);

            x_hat(1) = obj.x_hat(1);
            x_hat(2) = obj.x_hat(2);
            x_hat(3) = obj.x_hat(3);
            x_hat(4) = obj.x_hat(4);
            p_ball_obs = obj.x_hat(1);
            v_ball_obs = obj.x_hat(2);
            theta_obs = obj.x_hat(3);
            dtheta_obs = obj.x_hat(4);

            x_eq = obj.getEqPoint(t);
            [V_servo, theta_d] = obj.LQRController(x_hat, x_eq);
           
            V_servo = min(max(V_servo, -2), 3);
            
            % (OPTIONAL) Define safe limits for theta
            theta_min = -3*pi/8;  % Example lower bound
            theta_max = 3*pi/8;   % Example upper bound
            gain = 10;  % Scaling factor for control
            if theta < theta_min
                V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction
            elseif theta > theta_max
                V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction
            end

            obj.u = V_servo;
            
        end
    end
    
    methods(Access = public)

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

        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d, x_hat] = stepController(obj, t, p_ball, theta)        
            [V_servo, x_hat] = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end

    end
    
end
