classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        
        u_fn;
        kp = 100;
        kd = 1;

    end
    methods(Access = protected)
        

        function V_servo = stepImpl(obj, t, p_ball, v_ball, theta, dtheta)
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

            e_pos = p_ball - p_ball_ref;
            e_vel = v_ball - v_ball_ref;

            v = - e_pos * obj.kp - e_vel * obj.kd;

            V_servo = obj.u_fn(p_ball, v_ball, theta, dtheta, v);

            % Define safe limits for theta
            theta_min = -3*pi/8;  % Example lower bound
            theta_max = 3*pi/8;   % Example upper bound
            % Define a proportional gain (tune this value)
            K = 10;  % Scaling factor for control
            % Apply scaled control correction
            if theta < theta_min
                V_servo = max(V_servo, K * (theta_min - theta));  % Proportional correction
            elseif theta > theta_max
                V_servo = min(V_servo, K * (theta_max - theta));  % Proportional correction
            end

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
            
            Lfh = simplify(jacobian(h, [x1; x2; x3; x4]) * (f+g_vec*u))
            Lf2h = simplify(jacobian(Lfh, [x1; x2; x3; x4]) * (f+g_vec*u))
            Lf3h = simplify(jacobian(Lf2h, [x1; x2; x3; x4]) * (f+g_vec*u))
            
            Lf3h = expand(Lf3h)
            
            Lf3h_discard = Lf3h - (- (5*K*rg^2*u*x4*cos(x3)^2)/(7*L*tau) + (10*K*rg^2*u*x1*x4*cos(x3)^2)/(7*L^2*tau))
            Lf3h_discard = simplify(Lf3h_discard)
            
            Lf4h = jacobian(Lf3h_discard, [x1; x2; x3; x4]) * (f+g_vec*u)
            
            expand(Lf4h)
            
            u_sol = simplify(solve(Lf4h == v, u))

            rg_val = 0.0254;
            L_val = 0.4255;
            g_val = 9.81;
            K_val = 1.5;
            tau_val = 0.025;
            u_sol = subs(u_sol, {g, rg, L, K, tau}, {g_val, rg_val, L_val, K_val, tau_val})
            
            % Convert to a function handle for numerical evaluation
            obj.u_fn = matlabFunction(u_sol, 'Vars', [x1, x2, x3, x4, v]);

        end

        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, v_ball, theta, dtheta)        
            V_servo = stepImpl(obj, t, p_ball, v_ball, theta, dtheta);
            theta_d = obj.theta_d;
        end
    end
    
end
