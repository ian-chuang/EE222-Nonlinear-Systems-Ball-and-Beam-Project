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

        x_hat = [-0.19; 0.00; 0; 0];
        u = 0;

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

    end
    methods(Access = protected)
        function [V_servo, x_hat] = stepImpl(obj, t, p_ball, v_ball, theta, dtheta)
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
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            obj.x_hat = obj.extendedLuenbergerObserver(obj.x_hat, obj.u, [p_ball;theta;]);
            x_hat = obj.x_hat;

            % disp(obj.x_hat);
            p_ball_obs = obj.x_hat(1);
            v_ball_obs = obj.x_hat(2);
            theta_obs = obj.x_hat(3);
            dtheta_obs = obj.x_hat(4);

            V_servo = obj.feedbackLinearizationController(p_ball, v_ball, theta, dtheta, ...
                p_ball_ref, v_ball_ref, a_ball_ref);
            % V_servo = obj.feedbackLinearizationController(p_ball_obs, v_ball_obs, theta_obs, dtheta_obs, ...
            %     p_ball_ref, v_ball_ref, a_ball_ref);
            % disp(V_servo);
            disp([p_ball, v_ball, theta, dtheta]);
            disp([p_ball_obs, v_ball_obs, theta_obs, dtheta_obs]);
            disp([V_servo]);
            disp("-----");
            

            % Define safe limits for theta
            theta_min = -3*pi/8;  % Example lower bound
            theta_max = 3*pi/8;   % Example upper bound
            gain = 10;  % Scaling factor for control
            if theta < theta_min
                V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction
            elseif theta > theta_max
                V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction
            end

            obj.u = V_servo;

            % % Decide desired servo angle based on simple proportional feedback.
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
            
            % Update class properties if necessary.
            obj.t_prev = t;
            % obj.theta_d = theta_d;
        end
    end
    
    methods(Access = public)

        function x_hat_next = extendedLuenbergerObserver(obj, x_hat_curr, u_curr, y_next)
            rg = obj.rg_val;
            L = obj.L_val;
            g = obj.g_val;
            K = obj.K_val;
            tau = obj.tau_val;

            f = @(x, u) [ 
                x(2); 
                ((5 * g * rg)/(7 * L)) * sin(x(3)) - (5/7) * ((L/2) - x(1)) * ((rg/L) * x(4))^2 * cos(x(3))^2;
                x(4); 
                -1 * x(4)/tau + (K/tau) * u
                ];

            h = @(x) [
                x(1); 
                x(3)
                ]; % ball position and servo angle

            syms x1 x2 x3 x4 u_sym real
            
            x_sym = [x1; x2; x3; x4];
            f_sym_continuous = subs(f(x_sym, u_sym));
            
            A_sym = jacobian(f_sym_continuous, x_sym);
            A_func = matlabFunction(A_sym, 'Vars', {x_sym, u_sym});
            
            C_func = @(x) [1 0 0 0; 0 0 1 0];

            hat_x = x_hat_curr;
            hat_u = u_curr;
        
            A_eval = A_func(hat_x, hat_u);
            C_eval = C_func(hat_x);
        
            Co = ctrb(A_eval', C_eval');
            rank_Co = rank(Co);
            % disp(['Rank of Controllability Matrix: ', num2str(rank_Co)]);
        
            if rank_Co < 4
                poles = [-4, -4.5, -5, -5.5];
            else
                poles = [-1, -1.5, -2, -2.5];
            end
        
            L = place(A_eval', C_eval', poles)';
            hat_y = h(hat_x);
            
            y = y_next;
            
            dot_hat_x = f(hat_x, hat_u) + L * (y - hat_y);
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

        function obj = studentControllerInterface
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

        end

        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d, x_hat] = stepController(obj, t, p_ball, v_ball, theta, dtheta)        
            [V_servo, x_hat] = stepImpl(obj, t, p_ball, v_ball, theta, dtheta);
            theta_d = obj.theta_d;
        end
    end
    
end
