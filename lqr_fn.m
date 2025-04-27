function K = lqr_fn(A, B, Q, R)
%LQR_VIA_MATRIX_SIGN_FUNCTION Computes the LQR gain using matrix sign function method
%
%   K = LQR_VIA_MATRIX_SIGN_FUNCTION(A, B, Q, R) returns the optimal gain matrix K
%   for the continuous-time LQR problem:
%
%       minimize J = ∫ (x'Qx + u'Ru) dt
%       subject to dx/dt = Ax + Bu
%
%   Inputs:
%       A - System dynamics matrix (n x n)
%       B - Input matrix (n x m)
%       Q - State cost matrix (n x n), symmetric positive semi-definite
%       R - Input cost matrix (m x m), symmetric positive definite
%
%   Output:
%       K - Optimal state feedback gain matrix (m x n)

    % Invert R (assumes R is positive definite)
    R_inv = inv(R);
    G = B * R_inv * B';

    % Construct the Hamiltonian matrix Z
    Z = [ A      -G;
         -Q   -A' ];

    % Initialize matrix W for iteration
    W = Z;

    % Newton iteration to compute the matrix sign function
    for i = 1:1000
        W = W - 0.5 * (W - inv(W));
    end

    % Determine the size of the system
    n = size(A, 1);

    % Partition W into 4 submatrices
    W11 = W(1:n, 1:n);
    W12 = W(1:n, n+1:end);
    W21 = W(n+1:end, 1:n);
    W22 = W(n+1:end, n+1:end);

    % Solve for the unique positive semidefinite solution P to the Riccati equation
    M = [W12; W22 + eye(n)];
    N = [W11 + eye(n); W21];
    P = M \ (-N);

    % Compute the optimal LQR gain
    K = R_inv * B' * P;
end