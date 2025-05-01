classdef studentControllerInterface < matlab.System
    
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        rg_val = 0.0254;
        L_val = 0.4255;
        g_val = 9.81;
        K_val = 1.5;
        tau_val = 0.025;
        dt = 0.01

        max_V = 3;
        K = [100 143 52.1 10.2]; % feed back linearization gains
        Kx = [32.6074; 31.6228]; % observer gains
        Kz = [31.6283; 0.1745]; % observer gains

        
        xh = [0.0; 0; 0.0; 0];

        
        t_prev = 0;
        u = 0;  
        theta_d = 0;

        P = diag([0.5, 0.01, 0.5, 0.01]);
        xhat_prev = [0.0; 0; -57*pi/180; 0];
        Q = diag([1e-2,1e-2,1e-2,1e-2]);
        R = diag([3e-3,3e-3]);
    end
    
    methods(Access = protected)
        % function setupImpl(obj)
            % 
        % end

        function [V_servo, x1h, x2h, x3h, x4h] = stepImpl(obj, t, p_ball, theta)
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
        if obj.t_prev ~= 0
            obj.dt = t - obj.t_prev;
        end
        
        % EKF
        obj.EKF(p_ball, theta);
        
        % Feedback lin to follow ref traj
        [p_ref, v_ref, a_ref] = get_ref_traj(t);
        xr = [p_ref; v_ref; a_ref; 0];
        V_servo = obj.feedbackLinearizationController(obj.xh, xr); % + cos(theta)*0.25;%E(0.1 + 0.5*(-p_ball+0.2)/0.4);

        % Safety 
        theta_min = -3*pi/8;
        theta_max =  3*pi/8;
        gain = 10;
        if theta < theta_min
            V_servo = max(V_servo, gain*(theta_min - theta));
        elseif theta > theta_max
            V_servo = min(V_servo, gain*(theta_max - theta));
        end
        V_servo = max(min(V_servo, obj.max_V), -obj.max_V);
        
        % logging
        x1h = obj.xh(1);
        x2h = obj.xh(2);
        x3h = obj.xh(3);
        x4h = obj.xh(4);

        % record
        obj.u = V_servo;
        obj.theta_d = obj.xh(3);
        obj.t_prev = t;

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

        function EKF(obj, p_ball, theta)
            u_prev = obj.u;
            x_prev = obj.xh;
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
        end
    
        function u = feedbackLinearizationController(obj, x, xr)
    
            z_fn = @(x1,x2,x3,x4,x1r,x2r,x3r,x4r) [
                x1-x1r;
                x2-x2r;
                (124587*sin(x3))/297850 + (64516*x4^2*cos(x3)^2*((5*x1)/7 - 851/5600))/18105025 - x3r;
                (20320000*x4*cos(x3)*((108077*x4^2*sin(x3))/400000000 + (108077*x4*cos(x3))/10000000 - (127*x1*x4*cos(x3))/2500 + (127*x2*x4*cos(x3))/200000 - (127*x1*x4^2*sin(x3))/100000 + 834831/8000000))/5069407 - x4r;
            ];
            u_fn = @(x1,x2,x3,x4,v) (4000000000*((25698887331649*v)/12800000000000000 + (537476460632559*x4^2*sin(x3))/640000000000000000 - (1710053628273*x4^2*(sin(x3) - sin(x3)^3))/800000000000000000 + (69581560143053*x4^2*cos(x3)^2)/10000000000000000 - (69581560143053*x4^4*cos(x3)^2)/16000000000000000000 + (221383089491*x4^4*cos(x3)^4)/80000000000000000000 + (69581560143053*x4^3*sin(2*x3))/320000000000000000 - (81764465503*x1*x4^4)/8000000000000000 + (537476460632559*x4*cos(x3))/16000000000000000 + (69581560143053*x4^4)/32000000000000000000 - (81764465503*x1*x4^2*cos(x3)^2)/2500000000000 + (81764465503*x2*x4^2*cos(x3)^2)/100000000000000 + (81764465503*x1*x4^4*cos(x3)^2)/4000000000000000 - (260144641*x1*x4^4*cos(x3)^4)/20000000000000000 - (81764465503*x1*x4^3*sin(2*x3))/80000000000000 + (81764465503*x2*x4^3*sin(2*x3))/8000000000000000))/(1931444067*cos(x3)*((324231*x4^2*sin(x3))/400000000 + (108077*x4*cos(x3))/5000000 - (127*x1*x4*cos(x3))/1250 + (127*x2*x4*cos(x3))/100000 - (381*x1*x4^2*sin(x3))/100000 + 834831/8000000));
    
            z = z_fn(x(1),x(2),x(3),x(4),xr(1),xr(2),xr(3),xr(4));
            v = - obj.K * z;
            u = u_fn(x(1), x(2), x(3), x(4), v);
    
        end




        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            [V_servo,x1h,x2h,x3h,x4h] = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end
