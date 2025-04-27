%% SETUP_BALL_BEAM
% Sets the necessary parameters to run the Ball and Beam experiment.
%
clear;

%% System Parameters
% Load Ball and Beam model parameters.
[ L, r_arm, r_b, m_b, J_b, g, THETA_OFF, THETA_MIN, THETA_MAX, K_BS] = config_bb( );

fprintf("L=%f\n",L)
fprintf("r_arm=%f\n", r_arm)
fprintf("r_b=%f\n", r_b)
fprintf("m_b=%f\n", m_b)
fprintf("J_b=%f\n", J_b)
fprintf("g=%f\n", g)
fprintf("THETA_OFF=%f\n", THETA_OFF)
fprintf("THETA_MIN=%f\n", THETA_MIN)
fprintf("THETA_MAX=%f\n", THETA_MAX)
fprintf("K_BS=%f\n", K_BS)

%% Calculate Control Parameters
% Load model parameters based on servo configuration.    
K = 10; % enter correct value
tau = 0.1; % enter correct value
K_bb = 0.82;

fprintf("K=%f\n", K)
fprintf("tau=%f\n", tau)
fprintf("K_bb=%f\n", K_bb)

