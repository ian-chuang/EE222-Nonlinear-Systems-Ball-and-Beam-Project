function plot_states_report(ts, xs, x_hats, ref_ps, ref_vs, theta_ds, name)

if nargin < 5
    theta_ds = []; 
end

fig = figure();
set(fig, 'DefaultAxesFontSize', 16);
set(fig, 'DefaultTextFontName', 'Times New Roman');
set(fig, 'DefaultTextInterpreter', 'latex');

% Define colors
c1 = [0, 0.4470, 0.7410];   % blue
c2 = [0.8500, 0.3250, 0.0980]; % red-orange
c3 = [0.4660, 0.6740, 0.1880]; % green
c4 = [0.4940, 0.1840, 0.5560]; % purple

subplot(4, 1, 1);
plot(ts, 100 * xs(1, :), 'Color', c1, 'LineWidth', 1.5);
hold on;
plot(ts, 100 * x_hats(1, :), 'Color', c2, 'LineWidth', 1.5);
plot(ts, 100 * ref_ps, '-.', 'Color', c3, 'LineWidth', 1.5);
ylabel("z [cm]", 'Interpreter', 'latex');
grid on;
title(['State History for ', name]);
legend({'$z$', '$\hat{z}$', '$z_\mathrm{ref}$'}, 'Interpreter', 'latex');

subplot(4, 1, 2);
plot(ts, 100 * xs(2, :), 'Color', c1, 'LineWidth', 1.5);
hold on;
plot(ts, 100 * x_hats(2, :), 'Color', c2, 'LineWidth', 1.5);
plot(ts, 100 * ref_vs, '-.', 'Color', c3, 'LineWidth', 1.5);
grid on;
ylabel('$\dot{z}$ [cm / s]', 'Interpreter', 'latex');
legend({'$\dot{z}$', '$\dot{\hat{z}}$', '$\dot{z}_\mathrm{ref}$'}, 'Interpreter', 'latex');

subplot(4, 1, 3);
plot(ts, 180 * xs(3, :) / pi, 'Color', c1, 'LineWidth', 1.5);
hold on;
plot(ts, 180 * x_hats(3, :) / pi, 'Color', c2, 'LineWidth', 1.5);
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
if ~isempty(theta_ds)
    plot(ts, 180 * theta_ds / pi, 'r:', 'LineWidth', 1.5);
    legend({'$\theta$', '$\hat{\theta}$', '$\theta_d$'}, 'Interpreter', 'latex');
else
    legend({'$\theta$', '$\hat{\theta}$'}, 'Interpreter', 'latex');
end
grid on;

subplot(4, 1, 4);
plot(ts, 180 * xs(4, :) / pi, 'Color', c1, 'LineWidth', 1.5);
hold on;
plot(ts, 180 * x_hats(4, :) / pi, 'Color', c2, 'LineWidth', 1.5);
ylabel('$\dot{\theta}$ [deg/s]', 'Interpreter', 'latex');
xlabel('$t$ [sec]', 'Interpreter', 'latex');
legend({'$\dot{\theta}$', '$\dot{\hat{\theta}}$'}, 'Interpreter', 'latex');
grid on;

end