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

        K = [ 10.0000   14.6947   10.7967    4.6469]; % feed back linearization gains
        Kx = [32.6074; 31.6228]; % observer gains
        Kz = [31.6283; 0.1745]; % observer gains

        t_prev = 0;
        theta_d = 0;
        xh = [-0.19; 0.0;0;0];
        u = 0; % Initial state estimate for observer
    end
    methods(Access = protected)
        % function setupImpl(obj)
        %    disp("You can use this function for initializaition.");
        % end

        function [V_servo,x1h,x2h,x3h,x4h] = stepImpl(obj, t, p_ball, theta)
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
            xh = obj.xh;
            u = obj.u;
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            
            y1 = p_ball;
            y2 = theta;
            dt = t-t_prev;
            xh = obj.homeworkObserver(xh, u, y1, y2, dt);

            xr = [p_ball_ref; v_ball_ref; a_ball_ref; 0];
            V_servo = obj.feedbackLinearizationController(xh, xr);
            
            % SAFETY
            theta_min = -3*pi/8;  % Example lower bound
            theta_max = 3*pi/8;   % Example upper bound
            gain = 10;  % Scaling factor for control
            if theta < theta_min
                V_servo = max(V_servo, gain * (theta_min - theta));  % Proportional correction
            elseif theta > theta_max
                V_servo = min(V_servo, gain * (theta_max - theta));  % Proportional correction
            end
            V_servo = max(min(V_servo, 2), -2);
            x1h = xh(1);
            x2h = xh(2);
            x3h = xh(3);
            x4h = xh(4);

            % update attributes
            obj.t_prev = t;
            obj.u = V_servo;
            obj.xh = xh;
        end
    end

    
    
    methods(Access = public)


        function xh = homeworkObserver(obj, xh, u, y1, y2, dt)

            phi = @(x1, x3, x4) [ 
                0;
                (5*obj.g_val*obj.rg_val/(7*obj.L_val)) * sin(x3) - (5/7) * ((obj.L_val/2) - x1) * ((obj.rg_val/obj.L_val)^2) * x4^2 * cos(x3)^2
            ];
            
            A = [0.0 1; 0 0];        
            C = [1.0 0];             
            
            F = [0 1; 0 -1/obj.tau_val];   
            G = [0; obj.K_val/obj.tau_val];        
            H = [1.0 0];         
    
            xh_dot = [
                A * xh(1:2) + phi(y1,y2,xh(4)) + obj.Kx * (y1 - C * xh(1:2));
                F * xh(3:4) + G * u + obj.Kz * (y2 - H * xh(3:4))
            ];
    
            xh = xh + xh_dot*dt;
    
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
