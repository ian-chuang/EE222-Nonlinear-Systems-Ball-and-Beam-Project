function [A, B, f] = symbolic_dynamics()
    syms p_ball v_ball theta d_theta u;
    
    g = 9.81;
    r_arm = 0.0254;
    L = 0.4255;
    
    a = 5 * g * r_arm / (7 * L);
    b = (5 * L / 14) * (r_arm / L)^2;
    c = (5 / 7) * (r_arm / L)^2;
    
    K = 1.5;
    tau = 0.025;
    
    % dynamics
    dx =[
        v_ball
        a * sin(theta) - b * d_theta^2 * cos(theta)^2 + c * p_ball * d_theta^2 * cos(theta)^2
        d_theta
        (- d_theta + K * u) / tau
    ];
    
    As = jacobian(dx, [p_ball; v_ball; theta; d_theta]);
    Bs = jacobian(dx, [u]);
    
    A = matlabFunction(As, 'vars', [p_ball; v_ball; theta; d_theta; u]);
    B = matlabFunction(Bs, 'vars', [p_ball; v_ball; theta; d_theta; u]);
    f = matlabFunction(dx, 'vars', [p_ball; v_ball; theta; d_theta; u]);
end