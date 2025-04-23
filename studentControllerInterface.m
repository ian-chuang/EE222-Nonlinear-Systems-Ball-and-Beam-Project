
classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;

        K = [44 45 22 6];
        % x_hat = [-0.19; 0.00; 0; 0]; % Initial state estimate for observer
        x_hat = [0; 0.00; -57 * pi / 180; 0];
        u = 0; % Initial state estimate for observer
        elo_poles = [-4, -4.5, -5, -35.5];
        
    end
    methods(Access = protected)
        % function setupImpl(obj)
        %    disp("You can use this function for initializaition.");
        % end

        function [V_servo, p_ball_obs, theta_obs]  = stepImpl(obj, t, p_ball, theta)
            t_prev = obj.t_prev;
            %% Sample Controller: Simple Proportional Controller
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            obj.x_hat = obj.extendedLuenbergerObserver(obj.x_hat, obj.u, [p_ball;theta;]);

            p_ball_obs = obj.x_hat(1);
            v_ball_obs = obj.x_hat(2);
            theta_obs = obj.x_hat(3);
            dtheta_obs = obj.x_hat(4);

            V_servo = obj.feedbackLinearizationController(p_ball_obs, v_ball_obs, theta_obs, dtheta_obs, ...
                p_ball_ref, v_ball_ref, a_ball_ref);

            % (OPTIONAL) Define safe limits for theta
            theta_min = -3*pi/8;  % Example lower bound
            theta_max = 3*pi/8;   % Example upper bound
            gain = 10;  % Scaling factor for control
            if theta < theta_min
                V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction
            elseif theta > theta_max
                V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction
            end

            V_servo = max(min(V_servo, 2.5), -2.5);

            obj.u = V_servo;
            obj.t_prev = t;

            p_ball_obs = 100 * p_ball_obs;
            theta_obs = theta_obs * 180 / pi;
        end
    end

    



    
    methods(Access = public)
        function u_next = feedbackLinearizationController(obj, p_ball, v_ball, theta, dtheta, p_ball_ref, v_ball_ref, a_ball_ref)
            z1 = obj.h_fn(p_ball, v_ball, theta, dtheta) - p_ball_ref;
            z2 = obj.Lfh_fn(p_ball, v_ball, theta, dtheta) - v_ball_ref;
            z3 = obj.Lf2h_fn(p_ball, v_ball, theta, dtheta) - a_ball_ref;
            z4 = obj.Lf3h_fn(p_ball, v_ball, theta, dtheta) - 0;
            z = [z1; z2; z3; z4];
            v = - obj.K * z;
            u_next = obj.u_fn(p_ball, v_ball, theta, dtheta, v);
        end
        function u = u_fn(obj, x1, x2, x3, x4, v)
            u = (4000000000*((25698887331649*v)/12800000000000000 + (537476460632559*x4^2*sin(x3))/640000000000000000 - (1710053628273*x4^2*(sin(x3) - sin(x3)^3))/800000000000000000 + (69581560143053*x4^2*cos(x3)^2)/10000000000000000 - (69581560143053*x4^4*cos(x3)^2)/16000000000000000000 + (221383089491*x4^4*cos(x3)^4)/80000000000000000000 + (69581560143053*x4^3*sin(2*x3))/320000000000000000 - (81764465503*x1*x4^4)/8000000000000000 + (537476460632559*x4*cos(x3))/16000000000000000 + (69581560143053*x4^4)/32000000000000000000 - (81764465503*x1*x4^2*cos(x3)^2)/2500000000000 + (81764465503*x2*x4^2*cos(x3)^2)/100000000000000 + (81764465503*x1*x4^4*cos(x3)^2)/4000000000000000 - (260144641*x1*x4^4*cos(x3)^4)/20000000000000000 - (81764465503*x1*x4^3*sin(2*x3))/80000000000000 + (81764465503*x2*x4^3*sin(2*x3))/8000000000000000))/(1931444067*cos(x3)*((324231*x4^2*sin(x3))/400000000 + (108077*x4*cos(x3))/5000000 - (127*x1*x4*cos(x3))/1250 + (127*x2*x4*cos(x3))/100000 - (381*x1*x4^2*sin(x3))/100000 + 834831/8000000));
        end
        function h = h_fn(obj, x1, x2, x3, x4)
            h = x1;
        end
        function Lfh = Lfh_fn(obj, x1, x2, x3, x4)
            Lfh = x2;
        end
        function Lf2h = Lf2h_fn(obj, x1, x2, x3, x4)
            Lf2h = (124587*sin(x3))/297850 + (64516*x4^2*cos(x3)^2*((5*x1)/7 - 851/5600))/18105025;
        end
        function Lf3h = Lf3h_fn(obj, x1, x2, x3, x4)
            Lf3h = (20320000*x4*cos(x3)*((108077*x4^2*sin(x3))/400000000 + (108077*x4*cos(x3))/10000000 - (127*x1*x4*cos(x3))/2500 + (127*x2*x4*cos(x3))/200000 - (127*x1*x4^2*sin(x3))/100000 + 834831/8000000))/5069407;
        end


        function x_hat_next = extendedLuenbergerObserver(obj, x_hat_curr, u_curr, y_next)
            coder.extrinsic('place')
            dt = 0.01;
            hat_x = x_hat_curr;
            hat_u = u_curr;
        
            A_eval = obj.A_func(hat_x(1), hat_x(2), hat_x(3), hat_x(4));
            C_eval = obj.C_func(hat_x);
%             L = zeros(2, 4);
%             At = A_eval';
%             Ct = C_eval';
%             L(:) = place(At, Ct, obj.elo_poles);
%             L = L';
            % display(L);

            L = [ ...
                     9.0061,    0.9367; 
                    20.0456,    4.6994; 
                     0.0352,   -0.0061; 
                    -0.6457,  159.7919 ...
                ];

            hat_y = obj.h_func(hat_x);
            
            y = y_next;
            
            dot_hat_x = obj.f_func(hat_x, hat_u) + L * (y - hat_y);
            x_hat_next = hat_x + dt * dot_hat_x;
        end
        function A = A_func(obj, x1, x2, x3, x4)
            A = [                              0, 1,                                                                                                           0,                                                    0;...
(64516*x4^2*cos(x3)^2)/25347035, 0, (7535201836833413*cos(x3))/18014398509481984 - (129032*x4^2*cos(x3)*sin(x3)*((5*x1)/7 - 851/5600))/18105025, (129032*x4*cos(x3)^2*((5*x1)/7 - 851/5600))/18105025;...
                             0, 0,                                                                                                           0,                                                    1;...
                             0, 0,                                                                                                           0,                                                  -40];
        end
        function C = C_func(obj, hat_x)
            C = [1 0 0 0; 0 0 1 0];
        end
        function f = f_func(obj,x,u)
            rg = 0.0254;
            L = 0.4255;
            g = 9.81;
            K = 1.5;
            tau = 0.025;
            f = [ 
                x(2); 
                ((5 * g * rg)/(7 * L)) * sin(x(3)) - (5/7) * ((L/2) - x(1)) * ((rg/L) * x(4))^2 * cos(x(3))^2;
                x(4); 
                -1 * x(4)/tau + (K/tau) * u
                ];
        end
        function h = h_func(obj, x)
            h = [
                x(1); 
                x(3)
                ];
        end









        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d, obs_pos] = stepController(obj, t, p_ball, theta)        
            [V_servo, obs_pos] = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end




        
    end
    
end