function [K_fbl, res] = exportDynamicsFuncs()
    g = 9.81;
    rg = 0.0254;
    L = 0.4255;
    
    K = 1.5;
    tau = 0.025;
    motor_torque_const = 0/tau;
    f = @(x, u) [ 
        x(2); 
        ((5 * g * rg)/(7 * L)) * sin(x(3)) - (5/7) * ((L/2) - x(1)) * ((rg/L) * x(4))^2 * cos(x(3))^2;
        x(4); 
        -1 * x(4)/tau + (K/tau) * u
        ];
    % f = @(x, u) [
    %     x(2);
    %     (1/(m_b - (J_b/(r_b^2))))*(m_b*g*(rg/L)*sin(x(3)) - m_b*(L/2 - x(1))*((rg/L)*cos(x(3))*x(4)))^2;
    %     x(4);
    %     -x(4)/tau + (K/tau)*u + motor_torque_const*m_b*x(2)/cos(x(3));
    % ];
    h = @(x) [
        x(1); 
        x(3)
    ]; % ball position and servo angle

    syms x1 x2 x3 x4 u_sym u v lambda p1 p2 p3 p4 k11 k12 k13 k14 k21 k22 k23 k24 real
    p = [p1; p2; p3; p4];
    K_sym = [k11 k21; k12 k22;k13 k23; k14 k24];
    x_sym = [x1; x2; x3; x4];
    f_sym_continuous = subs(f(x_sym, u_sym));
    h_sym_continuous = subs(h(x_sym));
    f_fn = matlabFunction(f_sym_continuous, 'vars', {x_sym, u_sym}, 'File', 'f_fn.m')
    h_fn = matlabFunction(h_sym_continuous, 'vars', {x_sym, u_sym}, 'File', 'h_fn.m')
    A_sym = jacobian(f_sym_continuous, x_sym);
    A_func = matlabFunction(A_sym, 'Vars', {x_sym, u_sym}, 'File', 'A_func.m');
    B_sym = jacobian(f_sym_continuous, u_sym);
    B_func = matlabFunction(B_sym, 'Vars', {x_sym, u_sym}, 'File', 'B_func.m');
    C_sym = jacobian(h_sym_continuous, x_sym);
    C_func = matlabFunction(C_sym, 'Vars', {x_sym, u_sym}, 'File', 'C_func.m');

    disp(A_sym)
    disp(B_sym)
    disp(C_sym)
    % res = solve( ...
    %     det(A_sym-K_sym*C_sym - p(1) .* eye(4, 4)), ...
    %     det(A_sym-K_sym*C_sym - p(2) .* eye(4, 4)), ...
    %     det(A_sym-K_sym*C_sym - p(3) .* eye(4, 4)), ...
    %     det(A_sym-K_sym*C_sym - p(4) .* eye(4, 4)), ...
    % K_sym);
    
    % disp(res);
    % luenberger_gains = matlabFunction([res.k11(2) res.k21(2);res.k12(2) res.k22(2);res.k13(2) res.k23(2);res.k14(2) res.k24(2)], 'Vars', {x_sym, u_sym, p}, 'File', 'luenberger_gains_func.m');
    
    g_vec = f([0,0,0,0], [1])
    f = subs(f(x_sym, 0))

    disp(f)


    % f = [ 
    %     x2;
    %     (5*g*rg/(7*L)) * sin(x3) - (5/7) * ((L/2) - x1) * ((rg/L)^2) * x4^2 * cos(x3)^2;
    %     x4;
    %     -x4/tau
    % ]
    % g_vec = [0; 0; 0; K/tau]


    
    h = x1; 
    Lfh = simplify(jacobian(h, [x1; x2; x3; x4]) * (f+g_vec*u));
    Lf2h = simplify(jacobian(Lfh, [x1; x2; x3; x4]) * (f+g_vec*u));
    Lf3h = simplify(jacobian(Lf2h, [x1; x2; x3; x4]) * (f+g_vec*u));
    Lf3h_discard = simplify(expand(Lf3h) - (-(5*K*rg^2*u*x4*cos(x3)^2)/(7*L*tau) + (10*K*rg^2*u*x1*x4*cos(x3)^2)/(7*L^2*tau)));
    Lf4h = jacobian(Lf3h_discard, [x1; x2; x3; x4]) * (f+g_vec*u);
    u_sol = simplify(solve(Lf4h == v, u));

    u_fn = matlabFunction(u_sol, 'Vars', [x1, x2, x3, x4, v], 'File', 'u_fn.m');
    % h_fn = matlabFunction(h, 'Vars', [x1, x2, x3, x4], 'File', 'h_fn.m');
    Lfh_fn = matlabFunction(Lfh, 'Vars', [x1, x2, x3, x4], 'File', 'Lfh_fn.m');
    Lf2h_fn = matlabFunction(Lf2h, 'Vars', [x1, x2, x3, x4], 'File', 'Lf2h_fn.m');
    disp(jacobian(simplify(Lf2h), u))
    Lf3h_fn = matlabFunction(Lf3h_discard, 'Vars', [x1, x2, x3, x4], 'File', 'Lf3h_fn.m');

    A = zeros(4,4);  % Initialize a 4x4 zero matrix
    A(1:3,2:4) = eye(3);  % Set the top-right 3x3 block to identity
    B = [0; 0; 0; 1];
    Q = [1,0,0,0;...
    0,0,0,0;...
    0,0,0,0;...
    0,0,0,0];
    R = 0.001;
    [K_fbl,S,P] = lqr(A,B,Q,R);
end