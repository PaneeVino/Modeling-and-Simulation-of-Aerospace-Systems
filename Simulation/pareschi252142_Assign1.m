% Modeling and Simulation of Aerospace systems (2024/2025)
% Assignment # 1
% Author: Marcello Pareschi 
%% Ex 1
% Script for Exercise 1: Solving the Freudenstein Equation
% This script solves the Freudenstein equation for alpha given a range of beta values using
% Newton's method (with analytical and finite difference derivatives) and compares results
% against MATLAB's fzero. It includes analytical validation, error analysis, and visualization,
% as specified in Exercise 1 of the assignment (Page 4).

clearvars; close all; clc;  % Clear variables, close figures, and clear command window
plotStyle;                  % Apply custom plotting style (assumed to be a user-defined function)

% Define rod lengths for the Freudenstein equation (Page 4)
a1 = 10;  % Length of rod a1 [cm]
a2 = 13;  % Length of rod a2 [cm]
a3 = 8;   % Length of rod a3 [cm]
a4 = 10;  % Length of rod a4 [cm]

% Define the Freudenstein function and its derivatives
g = @(beta, alpha) a1/a2*cos(beta) - a1/a4*cos(alpha) - cos(beta-alpha) + ...
    (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4);  % Freudenstein equation f(alpha, beta) = 0
dg = @(beta, alpha) a1/a4*sin(alpha) - sin(beta - alpha);  % First derivative df/dalpha

% Define beta range and initialize solution arrays
beta_vec = linspace(0, 2*pi/3, 1e3);  % Beta values from 0 to 2π/3 [rad]
alpha_vec1 = nan(size(beta_vec));     % Solutions with initial guess α0 = -0.1 (Newton)
Alpha_vec1 = nan(size(beta_vec));     % Solutions with initial guess α0 = -0.1 (fzero)
alpha_vec2 = nan(size(beta_vec));     % Solutions with initial guess α0 = 2π/3 (Newton)
Alpha_vec2 = nan(size(beta_vec));     % Solutions with initial guess α0 = 2π/3 (fzero)

% Set solver parameters
tol = 1e-5;   % Convergence tolerance (Page 4)
Nmax = 1e3;   % Maximum iterations

% Solve for alpha using Newton's method and fzero over beta range
for i = 1 : length(beta_vec)
    beta = beta_vec(i);  % Current beta value [rad]
    f = @(alpha) g(beta, alpha);   % Fix beta in f(alpha)
    df = @(alpha) dg(beta, alpha); % Fix beta in df/dalpha

    % Solve with Newton's method (analytical derivative) and fzero
    [alpha_vec1(i), ~] = newton(f, df, -0.1, tol, Nmax);  % α0 = -0.1
    Alpha_vec1(i) = fzero(f, -0.1);                       % Reference solution
    [alpha_vec2(i), ~] = newton(f, df, 2*pi/3, tol, Nmax); % α0 = 2π/3
    Alpha_vec2(i) = fzero(f, 2*pi/3);                     % Reference solution
end

% Compute absolute errors between Newton and fzero
err_alpha1 = abs(Alpha_vec1 - alpha_vec1);  % Error for α0 = -0.1
err_alpha2 = abs(Alpha_vec2 - alpha_vec2);  % Error for α0 = 2π/3

% Compute maximum errors
err_alpha1_max = norm(err_alpha1, Inf);  % Max error for α0 = -0.1
err_alpha2_max = norm(err_alpha2, Inf);  % Max error for α0 = 2π/3

figure
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact')
% Plot solutions
nexttile
plot(beta_vec, alpha_vec1, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Solution $(\alpha_0=-0.1)$')
hold on
plot(beta_vec, alpha_vec2, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5, 'DisplayName', 'Solution $(\alpha_0=2\pi/3)$')
xlabel('$\beta$ [rad]'); ylabel('$\alpha$ [rad]')
title('NS solution', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

% Plot errors 
nexttile
semilogy(beta_vec, err_alpha1, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=-0.1)$')
hold on
semilogy(beta_vec, err_alpha2, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=2\pi/3)$')
xlabel('$\beta$ [rad]'); ylabel('Absolute Error [rad]')
title('Absolute error between NS and FZERO', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

% Solve using Newton's method with finite differences (centered scheme)
alpha_vec1_fd = nan(size(beta_vec));  % FD solutions with α0 = -0.1
alpha_vec2_fd = nan(size(beta_vec));  % FD solutions with α0 = 2π/3

for i = 1 : length(beta_vec)
    beta = beta_vec(i);  % Current beta value [rad]
    f = @(alpha) g(beta, alpha);  % Fix beta in f(alpha)

    % Solve with Newton's method (finite difference) and update fzero
    [alpha_vec1_fd(i), ~] = newton_fd(f, -0.1, tol, Nmax, 'cen');  % α0 = -0.1
    Alpha_vec1(i) = fzero(f, -0.1);                                % Reference solution
    [alpha_vec2_fd(i), ~] = newton_fd(f, 2*pi/3, tol, Nmax, 'cen'); % α0 = 2π/3
    Alpha_vec2(i) = fzero(f, 2*pi/3);                             % Reference solution
end

% Compute errors and discrepancies for finite difference solutions
err_alpha1_fd = abs(Alpha_vec1 - alpha_vec1_fd);  % Error for α0 = -0.1 (FD)
err_alpha2_fd = abs(Alpha_vec2 - alpha_vec2_fd);  % Error for α0 = 2π/3 (FD)
err_alpha1_fd_max = norm(err_alpha1_fd, Inf);     % Max error for α0 = -0.1 (FD)
err_alpha2_fd_max = norm(err_alpha2_fd, Inf);     % Max error for α0 = 2π/3 (FD)
disc_1 = abs(alpha_vec1 - alpha_vec1_fd);         % Discrepancy between NS and FD (α0 = -0.1)
disc_2 = abs(alpha_vec2 - alpha_vec2_fd);         % Discrepancy between NS and FD (α0 = 2π/3)

% Plot FD errors and discrepancies
figure
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(beta_vec, err_alpha1_fd, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=-0.1)$')
hold on
semilogy(beta_vec, err_alpha2_fd, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=2\pi/3)$')
xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('Absolute Error [rad]', 'FontSize', 20)
title('Absolute error between FD and FZERO', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

nexttile;
semilogy(beta_vec, disc_1, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=-0.1)$')
hold on
semilogy(beta_vec, disc_2, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5, 'DisplayName', 'Error $(\alpha_0=2\pi/3)$')
xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('Absolute Error [rad]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)
title('Absolute Error between FD and NS', 'FontSize', 20)

% Display maximum errors
fprintf('Err_alpha1_NS = %e \n', err_alpha1_max)  % Max error NS (α0 = -0.1)
fprintf('Err_alpha2_NS = %e \n', err_alpha2_max)  % Max error NS (α0 = 2π/3)
fprintf('Err_alpha1_FD = %e \n', err_alpha1_fd_max)  % Max error FD (α0 = -0.1)
fprintf('Err_alpha2_FD = %e \n', err_alpha2_fd_max)  % Max error FD (α0 = 2π/3)
fprintf('\n')

% Analytical solution for beta = 0 (Point 2, Page 5)
beta = 0;
f = @(alpha) g(beta, alpha);
df = @(alpha) dg(beta, alpha);
% Num = (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4) + a1 / a2;  % Numerator of analytical solution
% Den = a1 / a4 + 1;                                       % Denominator of analytical solution
Num = a1^2 + a2^2 - a3^2 + a4^2 + 2*a1*a4;  % Numerator of analytical solution
Den = 2*a2*a4 * (1 + a1/a4);                                       % Denominator of analytical solution
alpha_anal = acos(Num / Den);                           % Analytical alpha [rad]

% Numerical solutions for beta = 0
[alpha_NS_1, it_NS_1] = newton(f, df, -0.1, tol, Nmax);        % NS with α0 = -0.1
[alpha_FD_1, it_FD_1] = newton_fd(f, -0.1, tol, Nmax, 'cen');  % FD with α0 = -0.1
[alpha_NS_2, it_NS_2] = newton(f, df, 2*pi/3, tol, Nmax);      % NS with α0 = 2π/3
[alpha_FD_2, it_FD_2] = newton_fd(f, 2*pi/3, tol, Nmax, 'cen'); % FD with α0 = 2π/3

% Compute errors relative to analytical solution
err_NS_1 = abs(-alpha_anal - alpha_NS_1); % Error NS (α0 = -0.1, negative root)
err_FD_1 = abs(-alpha_anal - alpha_FD_1); % Error FD (α0 = -0.1, negative root)
err_NS_2 = abs(alpha_anal - alpha_NS_2);  % Error NS (α0 = 2π/3)
err_FD_2 = abs(alpha_anal - alpha_FD_2);  % Error FD (α0 = 2π/3)

% Display errors and iterations for beta = 0
fprintf('alpha_anal = %e \n', alpha_anal)
fprintf('err_NS_1 = %e \n', err_NS_1)
fprintf('err_FD_1 = %e \n', err_FD_1)
fprintf('err_NS_2 = %e \n', err_NS_2)
fprintf('err_FD_2 = %e \n', err_FD_2)
fprintf('it_NS_1 = %d \n', it_NS_1)
fprintf('it_FD_1 = %d \n', it_FD_1)
fprintf('it_NS_2 = %d \n', it_NS_2)
fprintf('it_FD_2 = %d \n', it_FD_2)

% Extended beta range for further analysis
beta_vec = linspace(0, pi, 100);  % Beta from 0 to π [rad]
alpha_vec1 = nan(size(beta_vec)); % Solutions with α0 = -0.1 (NS)
Alpha_vec1 = nan(size(beta_vec)); % Solutions with α0 = -0.1 (fzero)
alpha_vec2 = nan(size(beta_vec)); % Solutions with α0 = 2π/3 (NS)
Alpha_vec2 = nan(size(beta_vec)); % Solutions with α0 = 2π/3 (fzero)

% Solve over extended beta range
for i = 1 : length(beta_vec)
    beta = beta_vec(i);  % Current beta value [rad]
    f = @(alpha) g(beta, alpha);   % Fix beta in f(alpha)
    df = @(alpha) dg(beta, alpha); % Fix beta in df/dalpha

    alpha_vec1(i) = newton(f, df, -0.1, tol, Nmax);    % NS with α0 = -0.1
    Alpha_vec1(i) = fzero(f, -0.1);                    % Reference solution
    alpha_vec2(i) = newton(f, df, 2*pi/3, tol, Nmax);  % NS with α0 = 2π/3
    Alpha_vec2(i) = fzero(f, 2*pi/3);                 % Reference solution
end

% Compute errors for extended range
err_alpha_1 = abs(alpha_vec1 - Alpha_vec1);  % Error for α0 = -0.1
err_alpha_2 = abs(alpha_vec2 - Alpha_vec2);  % Error for α0 = 2π/3

% Plot solutions and errors for extended range
figure
tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile; plot(beta_vec, alpha_vec1, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5)
xlim([0, 3.5]); xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('$\alpha$ [rad]', 'FontSize', 20)
title('NS solution ($\beta \in \left[0, \pi\right]$)', 'FontSize', 20)
nexttile; plot(beta_vec, alpha_vec2, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5)
xlim([0, 3.5]); xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('$\alpha$ [rad]', 'FontSize', 20)
title('NS solution ($\beta \in \left[0, \pi\right]$)', 'FontSize', 20)
nexttile; plot(beta_vec, err_alpha_1, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3.5)
xlim([0, 3.5]); xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('Absolute Error [rad]', 'FontSize', 20)
title('Absolute error between NS and FZERO ($\beta \in \left[0, \pi\right]$)', 'FontSize', 20)
nexttile; plot(beta_vec, err_alpha_2, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 3.5)
xlim([0, 3.5]); xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('Absolute Error [rad]', 'FontSize', 20)
title('Absolute error between NS and FZERO ($\beta \in \left[0, \pi\right]$)', 'FontSize', 20)

% Contour plot for solution locus and maximum beta
alpha_vector = linspace(-0.5, pi, 1000);  % Alpha range [rad]
beta_vector = linspace(0, pi, 1000);      % Beta range [rad]
[B, A] = meshgrid(beta_vector, alpha_vector);  % Grid for contour
C = g(B, A);  % Evaluate Freudenstein function over grid
%%
figure
hold on; box on; grid on
M = contour(B, A, C, [0, 0], 'ShowText', 'on', 'LineWidth', 2.5, 'DisplayName', 'Locus of Solutions');
B_vec = M(1, 2:end);  % Extract beta values from contour
beta_max = max(B_vec);  % Maximum feasible beta [rad]
xline(beta_max, '--', 'LineWidth', 2.5, 'DisplayName', '$\beta_{max} \approx 2.6339$')  % Plot beta_max (~2.6339)
legend('Orientation', 'vertical', 'Location', 'best', 'FontSize', 20)
xlabel('$\beta$ [rad]', 'FontSize', 20); ylabel('$\alpha$ [rad]', 'FontSize', 20)
ylim([-0.5, 2])
%% Ex 2
% Script for Exercise 2: Solving ODE for Moving Object with Decreasing Mass
% This script solves the ODE dv/dt = (F(v) + f(t) - alpha*v)/m(t) using RK2 and RK4 methods,
% compares results against an analytical solution, and analyzes accuracy and execution time for
% different step sizes. It also tests stability with high fluid density, as per Exercise 2 (Page 6).

clearvars; close all; clc;  % Clear variables, close figures, and clear command window
plotStyle;                  % Apply custom plotting style (assumed user-defined function)

% Define color palette for plotting
Colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% Define physical parameters (Page 6)
m0 = 20;      % Initial mass [kg]
cm = 0.1;     % Mass loss rate [kg/s]
m = @(t) m0 - cm * t;  % Mass as a function of time [kg]

rho = 0;      % Fluid density (initially 0 for Point 1) [kg/m^3]
Cd = 2.05;    % Drag coefficient [dimensionless]
Am = 1;       % Cross-sectional area [m^2]
F = @(v) -0.5 * rho * Cd * Am * v.^2;  % Drag force (zero when rho = 0) [N]
f = @(t) 1 * t.^0;  % External force (constant 1 N)
alpha = 0.01;  % Damping coefficient [s^-1]

% Define ODE and initial conditions
dvdt = @(t, v) (F(v) + f(t) - alpha * v) ./ m(t);  % ODE: dv/dt [m/s^2]
v0 = 0;       % Initial velocity [m/s]
tspan = [0 160];  % Time interval [s]

% Analytical solution (derived on Page 7)
v_anal = @(t) f(t)/alpha - (f(t)/alpha - v0) .* (1 - cm * t / m0) .^ (alpha/cm);  % [m/s]

% Define step sizes for testing (Point 3, Page 7)
h_vec = [50, 20, 10, 1];  % Time steps [s]

% Initialize structure array for RK2 solutions
Sols_RK2 = repmat(struct('tt', [], 'vv', [], 'delta', [], 'err_rel', [], 'max_err', [], 'vv_anal', [], 'h', []), length(h_vec), 1);

% Solve ODE using RK2 for each step size
for i = 1 : length(h_vec)
    h = h_vec(i);  % Current step size [s]
    [tt, vv] = heun_method(dvdt, tspan, v0, h);  % Solve with RK2
    vv_anal = v_anal(tt);  % Analytical solution at time points
    err_abs = abs(vv_anal - vv);  % Absolute error [m/s]
    err_rel = err_abs ./ vv_anal; % Relative error [dimensionless]
    err_rel(1) = 0;  % Set initial error to 0 (exact initial condition)

    % Store results in structure
    Sols_RK2(i).tt = tt;         % Time vector [s]
    Sols_RK2(i).vv = vv;         % Numerical velocity [m/s]
    Sols_RK2(i).err_abs = err_abs;  % Absolute error [m/s]
    Sols_RK2(i).err_rel = err_rel;  % Relative error
    Sols_RK2(i).vv_anal = vv_anal;  % Analytical velocity [m/s]
    Sols_RK2(i).max_err = err_rel(end);  % Final relative error
    Sols_RK2(i).h = h;           % Step size [s]
end

% Plot RK2 solutions vs analytical
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
for i = 1 : length(Sols_RK2)
    tt = Sols_RK2(i).tt;  % Time vector [s]
    vv = Sols_RK2(i).vv;  % Numerical velocity [m/s]
    vv_anal = Sols_RK2(i).vv_anal;  % Analytical velocity [m/s]
    h = Sols_RK2(i).h;    % Step size [s]
    title_str = sprintf('RK2 vs Analytical Solution (h=%d)', h);

    nexttile
    plot(tt, vv_anal, 'o-', 'Color', Colors(1, :), 'LineWidth', 2.5, ...
      'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
      'MarkerSize', 5, 'DisplayName', 'RK2')
    hold on
    plot(tt, vv, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'Sol. Anal.')
    xlabel('time [s]', 'FontSize', 20); ylabel('velocity [m/s]', 'FontSize', 20)
    legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)
    title(title_str)
end

% Plot RK2 absolute errors
figure
tt = Sols_RK2(1).tt; err_abs = Sols_RK2(1).err_abs; Color = Colors(1, :); h = Sols_RK2(1).h;
DisplayName = sprintf('Abs. Err. (h=%d)', h);
semilogy(tt, err_abs, 'Color', Color, 'LineWidth', 2.5, 'DisplayName', DisplayName)
hold on
for i = 2 : length(Sols_RK2)
    tt = Sols_RK2(i).tt; err_abs = Sols_RK2(i).err_abs; h = Sols_RK2(i).h; Color = Colors(i, :);
    DisplayName = sprintf('Abs. Err. (h=%d)', h);
    semilogy(tt, err_abs, 'Color', Color, 'LineWidth', 2.5, 'DisplayName', DisplayName)
end
legend('Orientation', 'vertical', 'Location', 'southeast')
xlabel('time [s]', 'FontSize', 20); ylabel('Absolute Error [m/s]', 'FontSize', 20)

% Initialize structure array for RK4 solutions
Sols_RK4 = repmat(struct('tt', [], 'vv', [], 'delta', [], 'err_rel', [], 'max_err', [], 'vv_anal', [], 'h', []), length(h_vec), 1);

% Solve ODE using RK4 for each step size
for i = 1 : length(h_vec)
    h = h_vec(i);  % Current step size [s]
    [tt, vv] = RK4(dvdt, tspan, v0, h);  % Solve with RK4
    vv_anal = v_anal(tt);  % Analytical solution at time points
    err_abs = abs(vv_anal - vv);  % Absolute error [m/s]
    err_rel = err_abs ./ vv_anal; % Relative error [dimensionless]
    err_rel(1) = 0;  % Set initial error to 0 (exact initial condition)

    % Store results in structure
    Sols_RK4(i).tt = tt;         % Time vector [s]
    Sols_RK4(i).vv = vv;         % Numerical velocity [m/s]
    Sols_RK4(i).err_abs = err_abs;  % Absolute error [m/s]
    Sols_RK4(i).err_rel = err_rel;  % Relative error
    Sols_RK4(i).vv_anal = vv_anal;  % Analytical velocity [m/s]
    Sols_RK4(i).max_err = err_rel(end);  % Final relative error
    Sols_RK4(i).h = h;           % Step size [s]
end

% Plot RK4 solutions vs analytical
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
for i = 1 : length(Sols_RK4)
    tt = Sols_RK4(i).tt;  % Time vector [s]
    vv = Sols_RK4(i).vv;  % Numerical velocity [m/s]
    vv_anal = Sols_RK4(i).vv_anal;  % Analytical velocity [m/s]
    h = Sols_RK4(i).h;    % Step size [s]
    title_str = sprintf('RK4 vs Analytical Solution (h=%d)', h);

    nexttile
    plot(tt, vv, 'o-', 'Color', Colors(1, :), 'LineWidth', 2.5, ...
      'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
      'MarkerSize', 5, 'DisplayName', 'RK4')
    hold on
    plot(tt, vv_anal, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'Sol. Anal.')
    xlabel('time [s]', 'FontSize', 20); ylabel('velocity [m/s]', 'FontSize', 20)
    legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)
    title(title_str)
end
%%
% Plot RK4 absolute errors
figure
tt = Sols_RK4(1).tt; err_abs = Sols_RK4(1).err_abs; Color = Colors(1, :); h = Sols_RK4(1).h;
DisplayName = sprintf('Abs. Err. (h=%d)', h);
semilogy(tt, err_abs, 'Color', Color, 'LineWidth', 2.5, 'DisplayName', DisplayName)
hold on
for i = 2 : length(Sols_RK4)
    tt = Sols_RK4(i).tt; err_abs = Sols_RK4(i).err_abs; h = Sols_RK4(i).h; Color = Colors(i, :);
    DisplayName = sprintf('Abs. Err. (h=%d)', h);
    semilogy(tt, err_abs, 'Color', Color, 'LineWidth', 2.5, 'DisplayName', DisplayName)
end
legend('Orientation', 'vertical', 'Location', 'southeast')
xlabel('time [s]', 'FontSize', 20); ylabel('Absolute Error [m/s]', 'FontSize', 20)

% Performance analysis: execution time and error (Point 4, Page 7)
N_run = 1e5; N_h = length(h_vec);  % Number of runs and step sizes
t_heun = nan(N_run, N_h); err_heun = nan(N_run, N_h);  % RK2 time and error arrays
t_RK4 = nan(N_run, N_h); err_RK4 = nan(N_run, N_h);    % RK4 time and error arrays

for ii = 1 : N_run
    for jj = 1 : N_h
        h = h_vec(jj);  % Current step size [s]
        % RK2 timing and error
        tic
        [tt, vv] = heun_method(dvdt, tspan, v0, h);
        Toc = toc;
        vv_anal = v_anal(tt); err_abs = abs(vv_anal - vv);
        t_heun(ii, jj) = Toc * 1e3;  % Execution time [ms]
        err_heun(ii, jj) = err_abs(end);  % Final absolute error [m/s]

        % RK4 timing and error
        tic
        [tt, vv] = RK4(dvdt, tspan, v0, h);
        Toc = toc;
        vv_anal = v_anal(tt); err_abs = abs(vv_anal - vv);
        t_RK4(ii, jj) = Toc * 1e3;  % Execution time [ms]
        err_RK4(ii, jj) = err_abs(end);  % Final absolute error [m/s]
    end
end

% Compute statistics for performance
t_heun_mean = nan(size(h_vec)); t_heun_std = nan(size(h_vec));
t_RK4_mean = nan(size(h_vec)); t_RK4_std = nan(size(h_vec));
err_heun_mean = nan(size(h_vec)); err_heun_std = nan(size(h_vec));
err_RK4_mean = nan(size(h_vec)); err_RK4_std = nan(size(h_vec));

for i = 1 : length(h_vec)
    t_heun_vec = t_heun(:, i); t_RK4_vec = t_RK4(:, i);
    err_heun_vec = err_heun(:, i); err_RK4_vec = err_RK4(:, i);

    t_heun_mean(i) = mean(t_heun_vec); t_heun_std(i) = std(t_heun_vec);
    t_RK4_mean(i) = mean(t_RK4_vec); t_RK4_std(i) = std(t_RK4_vec);
    err_heun_mean(i) = mean(err_heun_vec); err_heun_std(i) = std(err_heun_vec);
    err_RK4_mean(i) = mean(err_RK4_vec); err_RK4_std(i) = std(err_RK4_vec);
end

% Display performance statistics
fprintf('t_mean_heun: %e\n', t_heun_mean)  % Mean execution time RK2 [ms]
fprintf('t_std_heun: %e\n', t_heun_std)   % Std dev execution time RK2 [ms]
fprintf('\n')
fprintf('t_mean_RK4: %e\n', t_RK4_mean)   % Mean execution time RK4 [ms]
fprintf('t_std_RK4: %e\n', t_RK4_std)     % Std dev execution time RK4 [ms]
fprintf('\n')
fprintf('err_mean_heun: %e\n', err_heun_mean)  % Mean final error RK2 [m/s]
fprintf('err_std_heun: %e\n', err_heun_std)    % Std dev final error RK2 [m/s]
fprintf('\n')
fprintf('err_mean_RK4: %e\n', err_RK4_mean)    % Mean final error RK4 [m/s]
fprintf('err_std_RK4: %e\n', err_RK4_std)      % Std dev final error RK4 [m/s]

% Plot execution time vs step size
figure
hold on; grid on; box on
plot(h_vec, t_heun_mean, '--o', 'DisplayName', 'RK2', 'Color', Colors(1, :), ...
    'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), 'LineWidth', 2.5)
plot(h_vec, t_RK4_mean, '--o', 'DisplayName', 'RK4', 'Color', Colors(2, :), ...
    'MarkerEdgeColor', Colors(2, :), 'MarkerFaceColor', Colors(2, :), 'LineWidth', 2.5)
legend('Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 20)
xlabel('time step [s]', 'FontSize', 20); ylabel('Execution time [ms]', 'FontSize', 20)
title('Execution time vs time step', 'FontSize', 20)

% Plot error vs step size with convergence rates
figure
loglog(h_vec, err_heun_mean, 'o--', 'DisplayName', 'RK2', 'Color', Colors(1, :), ...
    'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), 'LineWidth', 2.5)
hold on; grid on
loglog(h_vec, err_RK4_mean, 'o--', 'DisplayName', 'RK4', 'Color', Colors(2, :), ...
    'MarkerEdgeColor', Colors(2, :), 'MarkerFaceColor', Colors(2, :), 'LineWidth', 2.5)
loglog(h_vec, h_vec.^2, '--', 'DisplayName', '$h^2$', 'Color', Colors(3, :), 'LineWidth', 2.5)
loglog(h_vec, h_vec.^4, '--', 'DisplayName', '$h^4$', 'Color', Colors(4, :), 'LineWidth', 2.5)
legend('Location', 'north', 'Orientation', 'horizontal', 'FontSize', 20)
xlabel('time step [s]', 'FontSize', 20); ylabel('error [m/s]', 'FontSize', 20)
title('Numerical error vs time step', 'FontSize', 20)

% Test with high fluid density (Point 5, Page 8)
rho = 900;  % High fluid density [kg/m^3]
F = @(v) -0.5 * rho * Cd * Am * v.^2;  % Updated drag force [N]
dvdt = @(t, v) (F(v) + f(t) - alpha * v) ./ m(t);  % Updated ODE [m/s^2]
v0 = 0; tspan = [0 160];  % Reset initial conditions and time span

h = 1;  % Step size [s]
[tt_heun, vv_heun] = heun_method(dvdt, tspan, v0, h);  % RK2 solution
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);   % High precision for ode45
[tt, vv] = ode45(dvdt, tspan, v0, options);           % Reference solution

% Plot RK2 vs ode45 for h = 1 s
figure
plot(tt_heun, vv_heun, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'RK2 (h = 1 s)')
hold on
plot(tt, vv, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'ode45')
xlabel('time [s]', 'FontSize', 20); ylabel('velocity [m/s]', 'FontSize', 20)
legend('Location', 'east', 'Orientation', 'vertical', 'FontSize', 20)

h = 0.1;  % Smaller step size [s]
[tt_heun, vv_heun] = heun_method(dvdt, tspan, v0, h);  % RK2 solution

% Plot RK2 vs ode45 for h = 0.1 s
figure
plot(tt_heun, vv_heun, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', 'RK2 (h = 0.1 s)')
hold on
plot(tt, vv, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'ode45')
xlabel('time [s]', 'FontSize', 20); ylabel('velocity [m/s]', 'FontSize', 20)
legend('Location', 'east', 'Orientation', 'vertical', 'FontSize', 20)
%% Ex 3
% Script for Exercise 3: Stability Analysis of Numerical Methods
% This script computes and visualizes the stability regions of RK2, RK4, and BI2 methods for a
% linear test system dx/dt = A*x, where A has eigenvalues dependent on alpha. It is designed for
% Exercise 3 of the assignment (Page 9-11).

clearvars; close all; clc;  % Clear variables, close figures, and clear command window
plotStyle;                  % Apply custom plotting style (assumed user-defined function)

% Set options for fsolve (used to find maximum stable step size)
opts = optimoptions('fsolve', 'Display', 'none', 'FiniteDifferenceType', 'central', ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-14, 'Algorithm', 'levenberg-marquardt');

% Define color palette for plotting
Colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% Define alpha range and initialize arrays
alpha_vec = linspace(0, pi, 5e2).';  % Alpha values from 0 to π [rad]
h_vec = nan(size(alpha_vec));        % Maximum stable step sizes
eig_A = nan(length(alpha_vec), 2);   % Eigenvalues of system matrix A

% Compute eigenvalues of A (Page 9)
eig_A(:, 1) = exp(1i * alpha_vec);   % λ1 = e^(iα)
eig_A(:, 2) = exp(-1i * alpha_vec);  % λ2 = e^(-iα)

% Define system matrix A as a function of alpha
A = @(alpha) [0, 1; -1, 2*cos(alpha)];  % A from Equation 18 (Page 9)

% Define RK2 stability function (Page 10)
h0 = 2;  % Initial guess for step size
F_RK2 = @(h, alpha) 1 + h .* exp(1i * alpha) + h.^2 / 2 .* exp(2i * alpha);  % Stability polynomial

% Compute maximum stable step size for RK2
for k = length(alpha_vec) : -1 : 1
    alpha = alpha_vec(k);  % Current alpha [rad]
    fun = @(h) abs(F_RK2(h, alpha)) - 1;  % |R(hλ)| = 1 defines stability boundary
    h_vec(k) = fsolve(fun, h0, opts);     % Solve for h_max
    h0 = h_vec(k);                        % Update initial guess
end

% Plot h_max vs alpha for RK2
figure
plot(alpha_vec, h_vec, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', '$h_{max}$')
hold on
xline(pi/2, '--', 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', '$\alpha=\pi/2$')
xlabel('$\alpha$ [rad]', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$h$ [-]', 'FontSize', 20, 'Interpreter', 'latex')
legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)

% Plot RK2 stability region in complex plane
figure
hold on; grid on; axis equal
fill(real(h_vec .* eig_A), imag(h_vec .* eig_A), Colors(1, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
plot(real(h_vec .* eig_A), imag(h_vec .* eig_A), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('$Re(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$Im(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
xline(0, '--', 'Color', 'k', 'LineWidth', 1.5); yline(0, '--', 'Color', 'k', 'LineWidth', 1.5)
lgd = legend('RK2 stability region');
lgd.Orientation = 'vertical'; lgd.Location = 'northwest'; lgd.FontSize = 20; box on

% Compute analytical h_max for RK4 at alpha = π (Page 11)
a = 1/24; b = -1/6; c = 1/2; d = -1;  % Coefficients of RK4 stability polynomial
Delta_0 = b^2 - 3*a*c;
Delta_1 = 2*b^3 - 9*a*b*c + 27*a^2*d;
X = (Delta_1 + sqrt(Delta_1^2 - 4*Delta_0^3)) * 0.5;
C = nthroot(X, 3);
h0_pi_rk4 = -1 / (3*a) * (b + C + Delta_0 / C);  % h_max at α = π

% Define RK4 stability function and compute h_max
h_vec = nan(size(alpha_vec));  % Reset h_vec
h0 = h0_pi_rk4;  % Initial guess from analytical solution
F_RK4 = @(h, alpha) F_RK2(h, alpha) + h.^3 / 6 .* exp(3i * alpha) + h.^4 / 24 .* exp(4i * alpha);

for k = length(alpha_vec) : -1 : 1
    alpha = alpha_vec(k);  % Current alpha [rad]
    fun = @(h) abs(F_RK4(h, alpha)) - 1;  % |R(hλ)| = 1 defines stability boundary
    h_vec(k) = fsolve(fun, h0, opts);     % Solve for h_max
    h0 = h_vec(k);                        % Update initial guess
end

% Plot h_max vs alpha for RK4
figure
plot(alpha_vec, h_vec, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', '$h_{max}$')
hold on
xline(pi/2, '--', 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', '$\alpha=\pi/2$')
xlabel('$\alpha$ [rad]', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$h$ [-]', 'FontSize', 20, 'Interpreter', 'latex')
legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)

% Plot RK4 stability region in complex plane
figure
hold on; grid on; axis equal
fill(real(h_vec .* eig_A), imag(h_vec .* eig_A), Colors(1, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
plot(real(h_vec .* eig_A), imag(h_vec .* eig_A), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('$Re(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$Im(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
xline(0, '--', 'Color', 'k', 'LineWidth', 1.5); yline(0, '--', 'Color', 'k', 'LineWidth', 1.5)
lgd = legend('RK4 stability region');
lgd.Orientation = 'vertical'; lgd.Location = 'northwest'; lgd.FontSize = 20; box on

% Define stability functions for contour plotting
rk2 = @(z) 1 + z + z.^2 / 2;  % RK2 stability polynomial
rk4 = @(z) 1 + z + z.^2 / 2 + z.^3 / 6 + z.^4 / 24;  % RK4 stability polynomial

% Create grid for contour plot
x = linspace(-3, 3, 1e3); y = linspace(-3, 3, 1e3);
[X, Y] = meshgrid(x, y); Z = X + 1i * Y;  % Complex plane grid

% Evaluate stability functions
val_RK2 = rk2(Z); val_RK4 = rk4(Z);
abs_RK2 = abs(val_RK2); abs_RK4 = abs(val_RK4);
abs_RK4_masked = abs_RK4; abs_RK4_masked(abs_RK4 > 1) = nan;  % Mask outside stability region
abs_RK2_masked = abs_RK2; abs_RK2_masked(abs_RK2 > 1) = nan;

% Plot combined RK2 and RK4 stability regions
figure
M_rk4 = contour(X, Y, abs_RK4, [1, 1], 'LineWidth', 2.5, 'EdgeColor', Colors(1, :), 'FaceAlpha', 0.1);
axis equal; grid on; hold on
contour(X, Y, abs_RK4_masked, 150, 'LineWidth', 2.5, 'EdgeColor', 'none', 'FaceColor', Colors(1, :), 'FaceAlpha', 0.1)
M_rk2 = contour(X, Y, abs_RK2, [1, 1], 'LineWidth', 2.5, 'EdgeColor', Colors(2, :), 'FaceAlpha', 0.1);
contour(X, Y, abs_RK2_masked, 150, 'LineWidth', 2.5, 'EdgeColor', 'none', 'FaceColor', Colors(2, :), 'FaceAlpha', 0.1)
xlim([-3, 1])
xlabel('$Re(h\lambda)$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$Im(h\lambda)$', 'FontSize', 25, 'Interpreter', 'latex')
xline(0, '--', 'Color', 'k', 'LineWidth', 1.5); yline(0, '--', 'Color', 'k', 'LineWidth', 1.5)
lgd = legend('RK4 stability region', '', 'RK2 stability region', '');
lgd.Location = 'southwest'; lgd.Orientation = 'vertical'; lgd.FontSize = 21;

% Extract and plot h_max vs alpha from contour for RK4
real_rk4 = M_rk4(1, 2:end); imag_rk4 = M_rk4(2, 2:end);
h_vec_rk4 = sqrt(real_rk4.^2 + imag_rk4.^2);  % Magnitude of hλ
alpha_vec_rk4 = atan2(imag_rk4, real_rk4);    % Angle of hλ
figure
plot(alpha_vec_rk4, h_vec_rk4, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', '$h_{max}$')
xlim([0, 3.5]); hold on
xline(pi/2, '--', 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', '$\alpha=\pi/2$')
xlabel('$\alpha$ [rad]', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$h$ [-]', 'FontSize', 20, 'Interpreter', 'latex')
legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)

% Extract and plot h_max vs alpha from contour for RK2
real_rk2 = M_rk2(1, 2:end); imag_rk2 = M_rk2(2, 2:end);
h_vec_rk2 = sqrt(real_rk2.^2 + imag_rk2.^2);  % Magnitude of hλ
alpha_vec_rk2 = atan2(imag_rk2, real_rk2);    % Angle of hλ
figure
plot(alpha_vec_rk2, h_vec_rk2, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', '$h_{max}$')
xlim([0, 3.5]); hold on
xline(pi/2, '--', 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', '$\alpha=\pi/2$')
xlabel('$\alpha$ [rad]', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$h$ [-]', 'FontSize', 20, 'Interpreter', 'latex')
legend('Orientation', 'vertical', 'Location', 'northwest', 'FontSize', 20)

% Define BI2 stability function (Page 11)
h_vec = nan(size(alpha_vec));  % Reset h_vec
h0 = -10;  % Initial guess
Num = @(h, alpha, theta) 1 + h .* theta .* exp(1i * alpha) + (h .* theta).^2 / 2 .* exp(2i * alpha);
Den = @(h, alpha, theta) 1 - (1 - theta) .* h .* exp(1i * alpha) + ((1 - theta) .* h).^2 / 2 .* exp(2i * alpha);
F_BI2 = @(h, alpha, theta) Num(h, alpha, theta) ./ Den(h, alpha, theta);  % BI2 stability function

% Compute h_max for BI2 with theta = 0.3
for k = length(alpha_vec) : -1 : 1
    alpha = alpha_vec(k);  % Current alpha [rad]
    fun = @(h) abs(F_BI2(h, alpha, 0.3)) - 1;  % |R(hλ)| = 1 defines stability boundary
    h_vec(k) = fsolve(fun, h0, opts);          % Solve for h_max
    h0 = h_vec(k);                             % Update initial guess
end

% Plot BI2 stability region for theta = 0.3
figure
fill([-2 6 6 -2], [-4 -4 4 4], Colors(1, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Background
hold on; grid on; axis equal
fill(real(h_vec .* eig_A), imag(h_vec .* eig_A), [1 1 1], 'FaceAlpha', 0.9, 'EdgeColor', 'none')  % White out unstable region
plot(real(h_vec .* eig_A), imag(h_vec .* eig_A), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('$Re(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$Im(h\lambda)$', 'FontSize', 20, 'Interpreter', 'latex')
xline(0, '--', 'Color', 'k', 'LineWidth', 1.5); yline(0, '--', 'Color', 'k', 'LineWidth', 1.5)
lgd = legend('$BI2_{0.3}$ stability region');
lgd.Location = 'north'; lgd.FontSize = 20; box on
xlim([-2, 6]); ylim([-4, 4]);
%%
% Compute BI2 stability regions for multiple theta values
theta_vec = [0.2, 0.4, 0.6, 0.8];  % Theta values to test
h_vec = nan(length(alpha_vec), length(theta_vec));  % h_max for each alpha and theta
h0_vec = [-3, -10, 10, 3];  % Initial guesses for each theta

for ii = 1 : length(theta_vec)
    theta = theta_vec(ii);  % Current theta
    h0 = h0_vec(ii);        % Initial guess for this theta
    for k = length(alpha_vec) : -1 : 1
        alpha = alpha_vec(k);  % Current alpha [rad]
        fun = @(h) abs(F_BI2(h, alpha, theta)) - 1;  % |R(hλ)| = 1
        h_vec(k, ii) = fsolve(fun, h0, opts);        % Solve for h_max
        h0 = h_vec(k, ii);                           % Update initial guess
    end
end

% Plot BI2 stability regions for all theta values
h_fill = gobjects(1, length(theta_vec));
figure; hold on; grid on; axis equal
for ii = 1 : length(theta_vec)
    theta = theta_vec(ii);  % Current theta
    if theta > 0.5
        name = sprintf('$BI2_{%.1f}$ stability region', theta);  % Stable for theta > 0.5
    else
        name = sprintf('$BI2_{%.1f}$ instability region', theta); % Unstable for theta < 0.5
    end
    graf = fill(real(h_vec(:, ii) .* eig_A), imag(h_vec(:, ii) .* eig_A), Colors(ii, :), ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', name);
    set(graf(1), 'DisplayName', name); set(graf(2), 'HandleVisibility', 'off');  % Handle legend duplication
    plot(real(h_vec(:, ii) .* eig_A), imag(h_vec(:, ii) .* eig_A), 'Color', Colors(ii, :), ...
        'LineWidth', 3.5, 'HandleVisibility', 'off')
end
xlabel('$Re(h\lambda)$', 'FontSize', 34, 'Interpreter', 'latex')
ylabel('$Im(h\lambda)$', 'FontSize', 34, 'Interpreter', 'latex')
xline(0, '--', 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off')
yline(0, '--', 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off')
legend('Location', 'north', 'Orientation', 'vertical', 'FontSize', 28)
ax1 = gca; ax1.FontSize = 34;
%box on; xlim([-11 11]); ylim([-8 8])
%% Ex 4
% Script for Exercise 4: Temperature Evolution of a Cooling Body
% This script solves the ODE dT/dt = -Kc/C*T - Kr/C*T^4 + 1/C*(Kr*Ta^4 + Kc*Ta) using RK2 and RK4
% methods with different step sizes, compares results against ode89, and evaluates accuracy and
% function evaluations. It is designed for Exercise 4 of the assignment (Page 12-14).

clear; close all; clc;  % Clear workspace, close figures, and clear command window
plotStyle;              % Apply custom plotting style (assumed user-defined function)

% Define color palette for plotting
Colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% Define physical parameters (Page 12)
Kc = 0.0042;     % Convection coefficient [W/K]
Kr = 6.15e-11;   % Radiation coefficient [W/K^4]
C = 45;          % Thermal capacitance [J/K]
Ta = 277;        % Ambient temperature [K]
T0 = 555;        % Initial temperature [K]

% Define ODE for temperature evolution
dTdt = @(t, T) -Kc/C * T - Kr/C * T.^4 + 1/C * (Kr * Ta^4 + Kc * Ta);  % dT/dt [K/s]

% Define time span and step sizes
tspan = [0, 1440*10];  % Time span: 0 to 10 days [s]
h_1 = 720;             % Step size 1: 12 minutes [s]
h_2 = 1440;            % Step size 2: 24 minutes [s]

% Solve ODE using RK2 (Heun's method) with both step sizes
[t_heun_1, T_heun_1, fun_eval_heun_1] = heun_method(dTdt, tspan, T0, h_1);
[t_heun_2, T_heun_2, fun_eval_heun_2] = heun_method(dTdt, tspan, T0, h_2);

% Solve ODE using RK4 with both step sizes
[t_RK4_1, T_RK4_1, fun_eval_rk4_1] = RK4(dTdt, tspan, T0, h_1);
[t_RK4_2, T_RK4_2, fun_eval_rk4_2] = RK4(dTdt, tspan, T0, h_2);

% Solve ODE using ode89 as reference solution (Point 1, Page 13)
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);  % High precision settings
[t_ex, T_ex] = ode89(dTdt, tspan, T0, options);          % Exact solution [s, K]

% Plot temperature evolution
figure
plot(t_ex/3600, T_ex, 'Color', Colors(3, :), 'DisplayName', 'ode89', 'LineWidth', 3.5)
hold on; box on
plot(t_heun_1/3600, T_heun_1, 'o', 'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :), ...
    'MarkerSize', 5, 'DisplayName', 'RK2')
plot(t_RK4_2/3600, T_RK4_2, 'o', 'MarkerEdgeColor', Colors(2, :), 'MarkerFaceColor', Colors(2, :), ...
    'MarkerSize', 5, 'DisplayName', 'RK4')
yline(Ta, 'LineWidth', 3.5, 'Color', 'k', 'DisplayName', 'Air temperature')
yline(0.95*Ta, '--', 'LineWidth', 3.5, 'Color', 'k', 'DisplayName', 'Settling temperature')
yline(1.05*Ta, '--', 'LineWidth', 3.5, 'Color', 'k', 'HandleVisibility', 'off')
xlabel('t [h]', 'FontSize', 28); ylabel('T [K]', 'FontSize', 28)
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 25)

% Customize main plot axes
ax1 = gca; ax1.FontSize = 28;
ax1.XTick = linspace(0, max(t_ex/3600), 6);  % 6 evenly spaced ticks [h]
ax1.YTick = linspace(250, 550, 6);           % 6 evenly spaced ticks [K]

% Add inset plot for initial transient
axes('position', [0.53 0.3 0.33 0.44])  % Position: [left bottom width height]
box on
plot(t_ex(1:77)/3600, T_ex(1:77), 'Color', Colors(3, :), 'DisplayName', 'ode89', 'LineWidth', 3.5)
hold on
plot(t_heun_1(1:3)/3600, T_heun_1(1:3), 'o', 'MarkerEdgeColor', Colors(1, :), ...
    'MarkerFaceColor', Colors(1, :), 'MarkerSize', 5, 'DisplayName', 'RK2')
plot(t_RK4_2(1:2)/3600, T_RK4_2(1:2), 'o', 'MarkerEdgeColor', Colors(2, :), ...
    'MarkerFaceColor', Colors(2, :), 'MarkerSize', 5, 'DisplayName', 'RK4')
xlim([0, t_ex(77)/3600]); ylim([400 555])  % Zoom to initial 77 points

% Customize inset axes
ax2 = gca; ax2.FontSize = 28;
ax2.XTick = linspace(0, t_ex(77)/3600, 5);  % 5 ticks [h]
ax2.YTick = linspace(400, 555, 5);          % 5 ticks [K]

% Compute reference solutions at numerical time points for error analysis
[t_ex_heun_1, T_ex_heun_1] = ode89(dTdt, t_heun_1, T0, options);
[t_ex_heun_2, T_ex_heun_2] = ode89(dTdt, t_heun_2, T0, options);
[t_ex_RK4_1, T_ex_RK4_1] = ode89(dTdt, t_RK4_1, T0, options);
[t_ex_RK4_2, T_ex_RK4_2] = ode89(dTdt, t_RK4_2, T0, options);

% Compute relative errors (Point 2, Page 13)
err_heun_1 = abs(T_ex_heun_1 - T_heun_1.') ./ T_ex_heun_1;  % RK2, h = 720 s
err_heun_2 = abs(T_ex_heun_2 - T_heun_2.') ./ T_ex_heun_2;  % RK2, h = 1440 s
err_RK4_1 = abs(T_ex_RK4_1 - T_RK4_1.') ./ T_ex_RK4_1;     % RK4, h = 720 s
err_RK4_2 = abs(T_ex_RK4_2 - T_RK4_2.') ./ T_ex_RK4_2;     % RK4, h = 1440 s

% Plot relative errors
figure
plot(t_ex_heun_1/3600, err_heun_1, 'o', 'MarkerEdgeColor', Colors(1, :), ...
    'MarkerFaceColor', Colors(1, :), 'DisplayName', 'RK2')
hold on
plot(t_ex_RK4_2/3600, err_RK4_2, 'o', 'MarkerEdgeColor', Colors(2, :), ...
    'MarkerFaceColor', Colors(2, :), 'DisplayName', 'RK4')
xlabel('t [h]', 'FontSize', 20); ylabel('Relative Error [-]', 'FontSize', 20)
legend('Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 20)

% Display maximum relative errors and function evaluations (Point 3, Page 14)
fprintf('Relative error with RK2 (h = %d): %e\n', h_1, max(abs(err_heun_1)))
fprintf('Function Evaluations with RK2 (h = %d): %d\n\n', h_1, fun_eval_heun_1)
fprintf('Relative error with RK2 (h = %d): %e\n', h_2, max(abs(err_heun_2)))
fprintf('Function Evaluations with RK2 (h = %d): %d\n\n', h_2, fun_eval_heun_2)
fprintf('Relative error with RK4 (h = %d): %e\n', h_1, max(abs(err_RK4_1)))
fprintf('Function Evaluations with RK4 (h = %d): %d\n\n', h_1, fun_eval_rk4_1)
fprintf('Relative error with RK4 (h = %d): %e\n', h_2, max(abs(err_RK4_2)))
fprintf('Function Evaluations with RK4 (h = %d): %d\n\n', h_2, fun_eval_rk4_2)
%% Ex 5
% Script for Exercise 5: RLC Circuit Simulation and Stability Analysis
% This script solves the ODE system dx/dt = A*x for an RLC circuit using RK2 and IEX4 methods,
% compares stability regions, and evaluates accuracy and computational effort against an exact
% solution (expmv). It is designed for Exercise 5 of the assignment (Page 15-17).

clear; close all; clc;  % Clear workspace, close figures, and clear command window

% Set options for fsolve (used to find stability boundaries)
opts = optimoptions('fsolve', 'Display', 'none', 'FiniteDifferenceType', 'central', ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-14, 'Algorithm', 'levenberg-marquardt');

plotStyle;  % Apply custom plotting style (assumed user-defined function)

% Define color palette for plotting
Colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% Define RLC circuit parameters (Page 15)
R = 25;        % Resistance [Ohm]
L = 20e-3;     % Inductance [H]
C = 200e-3;    % Capacitance [F]
V0 = 12;       % Initial voltage [V]

% Define system matrix A
A = [0, 1; -1/(L*C), -R/L];  % A = [0 1; -1/LC -R/L] from Equation 25

% Compute eigenvalues of A
eigA = eig(A);  % Eigenvalues λ1 and λ2 [s^-1]

% Define initial conditions and ODE
xx0 = [V0*C; 0];  % x0 = [Q0; i0], Q0 = C*V0, i0 = 0 [Coulomb, Ampere]
dxdt = @(t, x) A * x;  % ODE: dx/dt = A*x
tmax = 40;  % Maximum time [s]
tspan = [0, tmax];  % Time span [s]

% Define test matrix AA for stability analysis (from Exercise 3)
AA = @(alpha) [0, 1; -1, 2*cos(alpha)];  % Test matrix with complex eigenvalues
alpha_vec = linspace(0, pi, 1e3).';  % Alpha range [rad]
eig_A = nan(length(alpha_vec), 2);   % Eigenvalues array

% Compute eigenvalues of AA
for k = 1 : length(alpha_vec)
    alpha = alpha_vec(k);
    eig_A(k,1) = cos(alpha) + 1i * sin(alpha);  % λ1 = e^(iα)
    eig_A(k,2) = cos(alpha) - 1i * sin(alpha);  % λ2 = e^(-iα)
end

% Compute RK2 stability boundary
h_vec_heun = nan(size(alpha_vec));  % Maximum stable step sizes for RK2
h0 = 2;  % Initial guess
for k = length(alpha_vec) : -1 : 1
    alpha = alpha_vec(k);  % Current alpha [rad]
    fun = @(h) abs(1 + h * exp(1i * alpha) + h.^2 / 2 * exp(2i * alpha)) - 1;  % |R(hλ)| = 1
    h_vec_heun(k) = fsolve(fun, h0, opts);  % Solve for h_max
    h0 = h_vec_heun(k);  % Update initial guess
end

% Compute IEX4 stability boundary
h_vec_IEX4 = nan(size(alpha_vec));  % Maximum stable step sizes for IEX4
h0 = -10;  % Initial guess
for k = length(alpha_vec) : -1 : 1
    alpha = alpha_vec(k);  % Current alpha [rad]
    fun = @(h) max(abs(eig(F_IEX4(h, alpha)))) - 1;  % |R(hλ)| = 1 for IEX4
    h_vec_IEX4(k) = fzero(fun, h0);  % Solve for h_max (using fzero instead of fsolve)
    h0 = h_vec_IEX4(k);  % Update initial guess
end

% Plot stability regions with RLC eigenvalues
figure
fill([-10 15 15 -10], [-8 -8 8 8], Colors(1, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
    'DisplayName', 'IEX4 Stability Region');  % Full IEX4 region
hold on
fill(real(h_vec_IEX4 .* eig_A), imag(h_vec_IEX4 .* eig_A), [1 1 1], 'FaceAlpha', 0.9, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off')  % White out unstable region
fill(real(h_vec_heun .* eig_A), imag(h_vec_heun .* eig_A), [1 1 1], 'FaceAlpha', 0.9, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off')  % White out outside RK2 region
plot(real(h_vec_IEX4 .* eig_A), imag(h_vec_IEX4 .* eig_A), 'Color', Colors(1, :), ...
    'LineWidth', 3.5, 'HandleVisibility', 'off')
graf = fill(real(h_vec_heun .* eig_A), imag(h_vec_heun .* eig_A), Colors(2, :), 'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', 'DisplayName', 'RK2 Stability Region');
set(graf(1), 'DisplayName', 'RK2 Stability Region'); set(graf(2), 'HandleVisibility', 'off');  % Fix legend duplication
plot(real(h_vec_heun .* eig_A), imag(h_vec_heun .* eig_A), 'Color', Colors(2, :), ...
    'LineWidth', 3.5, 'HandleVisibility', 'off')
xline(0, '--', 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off')
yline(0, '--', 'Color', 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off')
eig1 = '$\lambda_1=-0.2$'; eig2 = '$\lambda_2=-1249.8$';
plot(real(eigA(1)), imag(eigA(1)), 'o', 'MarkerEdgeColor', Colors(3, :), 'MarkerFaceColor', Colors(3, :), ...
    'DisplayName', eig1, 'MarkerSize', 10)  % Plot λ1
plot(real(eigA(2)), imag(eigA(2)), 'o', 'MarkerEdgeColor', Colors(4, :), 'MarkerFaceColor', Colors(4, :), ...
    'DisplayName', eig2, 'MarkerSize', 10)  % Plot λ2
grid on; axis equal; xlim([-10, 15]); ylim([-8, 8])
xlabel('$Re \left(h \, \lambda \right)$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$Im \left(h \, \lambda \right)$', 'FontSize', 25, 'Interpreter', 'latex')
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 30)

% Customize main plot axes
ax1 = gca; ax1.FontSize = 28;
ax1.XTick = linspace(-10, 15, 6); ax1.YTick = linspace(-8, 8, 5);

% Add inset for zoomed view of λ2
axes('position', [0.17 0.17 0.25 0.25])  % Position: [left bottom width height]
fill([-1251 -1249 -1249 -1251], [-1 -1 1 1], Colors(1, :), 'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', 'DisplayName', 'IEX4 Stability Region');
box on; hold on
yline(0, '--', 'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
plot(real(eigA(2)), imag(eigA(2)), 'o', 'MarkerEdgeColor', Colors(4, :), 'MarkerFaceColor', Colors(4, :), ...
    'MarkerSize', 10)
xlim([-1251, -1249]); ylim([-1, 1])
ax2 = gca; ax2.FontSize = 25;
ax2.XTick = linspace(-1251, -1249, 3); ax2.YTick = linspace(-1, 1, 3);
axis tight

% Solve ODE with IEX4
tic
[tt_iex, xx_iex, feval_iex] = IEX4_linear(A, tspan, xx0, 0.1);  % h = 0.1 s
t_cpu_iex = toc;  % CPU time [s]
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile; plot(tt_iex, xx_iex(1, :), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('time [s]'); ylabel('Q [C]')
nexttile; plot(tt_iex, xx_iex(2, :), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('time [s]'); ylabel('i [A]')

% Solve ODE with RK2
hmax_heun = -2 / eigA(2);  % Maximum stable step size for RK2 (real eigenvalue)
h_heun = 0.9 * hmax_heun;  % 90% of h_max for stability
tic
[tt_heun, xx_heun, feval_heun] = RK2_linear(A, tspan, xx0, h_heun);
t_cpu_rk2 = toc;  % CPU time [s]
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile; plot(tt_heun, xx_heun(1, :), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('time [s]'); ylabel('Q [C]')
nexttile; plot(tt_heun, xx_heun(2, :), 'Color', Colors(1, :), 'LineWidth', 2.5)
xlabel('time [s]'); ylabel('i [A]')

% Compute exact solution using matrix exponential (expmv)
xx_ex_iex = expmv(A, xx0, tt_iex);   % Exact solution at IEX4 time points
xx_ex_heun = expmv(A, xx0, tt_heun); % Exact solution at RK2 time points

% Compute relative errors
err_iex_mat = (xx_iex - xx_ex_iex) ./ xx_ex_iex;    % Error matrix for IEX4
err_heun_mat = (xx_heun - xx_ex_heun) ./ xx_ex_heun; % Error matrix for RK2
err_iex_q = abs(err_iex_mat(1, :));    % Relative error in charge (Q)
err_heun_q = abs(err_heun_mat(1, :));  % Relative error in charge (Q)
err_iex_i = abs(err_iex_mat(2, :));    % Relative error in current (i)
err_heun_i = abs(err_heun_mat(2, :));  % Relative error in current (i)

% Plot relative errors
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
semilogy(tt_iex(2:end), err_iex_q(2:end), 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'IEX4')
hold on
semilogy(tt_heun(2:end), err_heun_q(2:end), 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'RK2')
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 20)
xlabel('time [s]', 'FontSize', 20); ylabel('Relative error [-]', 'FontSize', 20)
title('Relative error in terms of Q', 'FontSize', 20)
nexttile
semilogy(tt_iex(2:end), err_iex_i(2:end), 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'IEX4')
hold on
semilogy(tt_heun(2:end), err_heun_i(2:end), 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'RK2')
legend('Location', 'northwest', 'Orientation', 'vertical', 'FontSize', 20)
xlabel('time [s]', 'FontSize', 20); ylabel('Relative error [-]', 'FontSize', 20)
title('Relative error in terms of i', 'FontSize', 20)
%% Ex 6
% Script for Exercise 6: Bouncing Ball Simulation and Stability Analysis
% This script simulates the motion of a bouncing ball using integrate_ball and RK4 methods,
% analyzes the effect of air density (rho) on the solution, and evaluates RK4 stability for the
% system. It is designed for Exercise 6 of the assignment (Page 18-20).

clear; close all; clc;  % Clear workspace, close figures, and clear command window
plotStyle;  % Apply custom plotting style (assumed user-defined function)

% Define color palette for plotting
Colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% Define initial conditions and time span
xx0 = [0; 10];  % Initial state: [velocity; position] = [0 m/s; 10 m]
tf = 10;        % Final time [s]

% Load parameters and set initial air density
par = par_ex_6();  % Load parameters from provided function (Page 18)
par.rho = 1.225;   % Air density [kg/m^3], standard atmospheric value

% Integrate with integrate_ball (assumed adaptive solver)
[tt, xx] = integrate_ball(xx0, tf, par);  % Solve ODE
pos = xx(2, :);  % Position [m]
vel = xx(1, :);  % Velocity [m/s]
nanindeces = find(isnan(vel));  % Find NaN indices (bounces or errors)

% Plot position and velocity for rho = 1.225
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(tt, pos, 'LineWidth', 2.5, 'Color', Colors(1, :))
xlabel('time [s]'); ylabel('position [m]'); ylim([0, 10])
nexttile
plot(tt, vel, 'LineWidth', 2.5, 'Color', Colors(2, :))
hold on
yline(0, '--', 'LineWidth', 2, 'Color', 'k')  % Zero velocity line
for k = 1 : length(nanindeces)  % Connect across NaN gaps (bounces)
    ind = nanindeces(k);
    if ind < length(tt)
        plot([tt(ind) tt(ind)], [vel(ind-1), vel(ind+1)], '--', 'LineWidth', 2, 'Color', Colors(2, :))
    end
end
xlabel('time [s]'); ylabel('velocity [m/s]')

% Increase air density to rho = 15
par.rho = 15;  % Higher air density [kg/m^3]
[tt, xx] = integrate_ball(xx0, tf, par);  % Re-solve ODE
pos = xx(2, :); vel = xx(1, :);
nanindeces = find(isnan(vel));

% Plot position and velocity for rho = 15
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(tt, pos, 'LineWidth', 2.5, 'Color', Colors(1, :))
xlabel('time [s]'); ylabel('position [m]'); ylim([0, 10])
nexttile
plot(tt, vel, 'LineWidth', 2.5, 'Color', Colors(2, :))
hold on
yline(0, '--', 'LineWidth', 2, 'Color', 'k')
for k = 1 : length(nanindeces)
    ind = nanindeces(k);
    if ind < length(tt)
        plot([tt(ind) tt(ind)], [vel(ind-1), vel(ind+1)], '--', 'LineWidth', 2, 'Color', Colors(2, :))
    end
end
xlabel('time [s]'); ylabel('velocity [m/s]')

% Increase air density to rho = 60
par.rho = 60;  % Very high air density [kg/m^3]
[tt, xx] = integrate_ball(xx0, tf, par);  % Re-solve ODE
pos = xx(2, :); vel = xx(1, :);
nanindeces = find(isnan(vel));

% Plot position and velocity for rho = 60
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(tt, pos, 'LineWidth', 2.5, 'Color', Colors(1, :))
xlabel('time [s]'); ylabel('position [m]'); ylim([0, 10])
nexttile
plot(tt, vel, 'LineWidth', 2.5, 'Color', Colors(2, :))
hold on
yline(0, '--', 'LineWidth', 2, 'Color', 'k')
for k = 1 : length(nanindeces)
    ind = nanindeces(k);
    if ind < length(tt)
        plot([tt(ind) tt(ind)], [vel(ind-1), vel(ind+1)], '--', 'LineWidth', 2, 'Color', Colors(2, :))
    end
end
xlabel('time [s]'); ylabel('velocity [m/s]'); ylim([-0.9, 0.1])

% Solve with RK4 using h = 3 s (Point 2, Page 19)
tspan = [0, tf];
h = 3;  % Step size [s]
fun = @(t, x) fun_bouncing(t, x, par);  % ODE function from assignment
[tt_rk4, xx_rk4] = RK4(fun, tspan, xx0, h);  % Solve with RK4
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(tt_rk4, xx_rk4(2,:), 'o--', 'LineWidth', 2.5, 'MarkerEdgeColor', Colors(1, :), ...
    'MarkerFaceColor', Colors(1, :), 'Color', Colors(1, :))
xlabel('time [s]'); ylabel('position [m]')
nexttile
plot(tt_rk4, xx_rk4(1,:), 'o--', 'LineWidth', 2.5, 'MarkerEdgeColor', Colors(2, :), ...
    'MarkerFaceColor', Colors(2, :), 'Color', Colors(2, :))
xlabel('time [s]'); ylabel('velocity [m]')

% Solve with RK4 using h = 0.5 s (Adjusted step size, Point 4)
tspan = [0, tf];
fun = @(t, x) fun_bouncing(t, x, par);
[tt_rk4, xx_rk4] = RK4(fun, tspan, xx0, 0.5);  % Solve with smaller h
pos = xx_rk4(2, :); vel = xx_rk4(1, :);
nanindeces = find(isnan(vel));

% Plot position and velocity for RK4 with h = 0.5 s
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(tt_rk4, pos, 'LineWidth', 2.5, 'Color', Colors(1, :))
xlabel('time [s]'); ylabel('position [m]'); ylim([0, 10])
nexttile
plot(tt_rk4, vel, 'LineWidth', 2.5, 'Color', Colors(2, :))
hold on
yline(0, '--', 'LineWidth', 2, 'Color', 'k')
for k = 1 : length(nanindeces)
    ind = nanindeces(k);
    if ind < length(tt)
        plot([tt(ind) tt(ind)], [vel(ind-1), vel(ind+1)], '--', 'LineWidth', 2, 'Color', Colors(2, :))
    end
end
xlabel('time [s]'); ylabel('velocity [m/s]'); ylim([-0.9, 0.1])
%% Functions 

function plotStyle

% Set LaTeX as the default interpreter for text, axes ticks, and legends
set(0, 'defaultTextInterpreter', 'Latex')           % Use LaTeX for general text (e.g., titles)
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')  % Use LaTeX for tick labels on axes
set(0, 'defaultLegendInterpreter', 'Latex')         % Use LaTeX for legend text

% Enable grid lines by default on both axes
set(0, 'defaultAxesXGrid', 'on')                    % Turn on x-axis grid lines
set(0, 'defaultAxesYGrid', 'on')                    % Turn on y-axis grid lines

% Set default line properties
set(0, 'defaultLineLineWidth', 1.5)                % Default line width to 1.5 points
set(0, 'defaultLineMarkerSize', 6)                 % Default marker size to 6 points
set(0, 'defaultLineMarkerEdgeColor', 'k')          % Default marker edge color to black
set(0, 'defaultLineMarkerFaceColor', 'auto')       % Default marker face color matches line color

% Set default legend properties
set(0, 'defaultLegendLocation', 'northoutside')    % Default legend position above plot
set(0, 'defaultLegendOrientation', 'horizontal')   % Default legend layout to horizontal
set(0, 'defaultLegendFontSize', 12)                % Default legend font size to 12 points

% Set default axes font size
set(0, 'defaultAxesFontSize', 16)                  % Default font size for axes labels/ticks to 16 points

end

function [x, iter] = newton(f, df, x0, tol, Nmax)
    % NEWTON Solves a nonlinear equation using Newton's method.
    % This function implements an iterative Newton solver to find a root of the equation f(x) = 0,
    % given its derivative and an initial guess. It is suitable for problems like the Freudenstein
    % equation in a kinematic chain (e.g., Exercise 1 of the assignment).
    %
    % Inputs:
    %   - f: Function handle representing the equation f(x) to solve (e.g., Freudenstein equation)
    %   - df: Function handle for the analytical derivative of f(x)
    %   - x0: Initial guess for the root [scalar, e.g., radians for angles]
    %   - tol: Tolerance for convergence (e.g., 10^-5)
    %   - Nmax: Maximum number of iterations to prevent infinite loops
    %
    % Outputs:
    %   - x: Approximate root of the equation f(x) = 0 [scalar]
    %   - iter: Number of iterations performed

    % Initialize the root estimate with the initial guess
    x = x0;

    % Initialize error as a large value to enter the loop
    err = 1e3;

    % Initialize iteration counter
    iter = 0;

    % Iterate until the error is below tolerance or max iterations are reached
    while err > tol && iter < Nmax
        % save previous value of x
        x_prec = x;
        % Update the estimate using Newton's method: x_{k+1} = x_k - f(x_k) / f'(x_k)
        x = x - f(x) / df(x);

        % Compute the new error as the difference between two successiveiterates
        err = abs(x-x_prec);

        % Increment the iteration counter
        iter = iter + 1;
    end
    
end

function [x, iter] = newton_fd(f, x0, tol, Nmax, scheme)
    % NEWTON_FD Solves a nonlinear equation using Newton's method with finite differences.
    % This function implements an iterative Newton solver to find a root of f(x) = 0, approximating
    % the derivative using finite difference schemes (forward or centered). It is designed for cases
    % where the analytical derivative is unavailable, such as Exercise 1, Point 3 of the assignment.
    %
    % Inputs:
    %   - f: Function handle representing the equation f(x) to solve (e.g., Freudenstein equation)
    %   - x0: Initial guess for the root [scalar, e.g., radians for angles]
    %   - tol: Tolerance for convergence (e.g., 10^-5)
    %   - Nmax: Maximum number of iterations to prevent infinite loops
    %   - scheme: String specifying the finite difference scheme
    %       'fw'  = Forward difference (first-order accurate)
    %       'cen' = Centered difference (second-order accurate)
    %
    % Outputs:
    %   - x: Approximate root of the equation f(x) = 0 [scalar]
    %   - iter: Number of iterations performed

    % Initialize the root estimate with the initial guess
    x = x0;

    % Initialize error as a large value to enter the loop
    err = 1e3;

    % Initialize iteration counter
    iter = 0;

    % Iterate until the error is below tolerance or max iterations are reached
    while err > tol && iter < Nmax
        % save previous value of x
        x_prec = x;
        % Compute step size for finite difference, balancing truncation and round-off errors
        h = sqrt(eps) * max([1, abs(x)]); % eps is machine epsilon (~2.2e-16 for double precision)

        % Evaluate function at current estimate
        a = f(x);

        % Compute approximate derivative based on specified scheme
        switch scheme
            case 'fw'
                % Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
                b = f(x + h);
                c = (b - a) / h;

            case 'bw'
                % Forward difference: f'(x) ≈ (f(h) - f(x-h)) / h
                b = f(x - h);
                c = (a - b) / h;

            case 'cen'
                % Centered difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
                b  = f(x + h);
                bb = f(x - h);
                c  = (b - bb) / (2*h);

            otherwise
                error('Only fw, bw and cen schemes!')
        end

        % Update the estimate using Newton's method: x_{k+1} = x_k - f(x_k) / f'(x_k)
        x = x - a / c;

        % Compute the new error as the difference between two successive iterates
        err = abs(x - x_prec);

        % Increment the iteration counter
        iter = iter + 1;
    end

end

function [tt, xx, fun_eval] = heun_method(f, tspan, xx0, h)
    % HEUN_METHOD Solves an ODE system using Heun's method (RK2).
    % This function implements a fixed-step second-order Runge-Kutta (Heun's) method to numerically
    % integrate a system of ordinary differential equations. It is designed for problems like the
    % moving object with decreasing mass in Exercise 2 of the assignment.
    %
    % Inputs:
    %   - f: Function handle for the right-hand side of the ODE system, f(t, x)
    %   - tspan: Time interval vector [t0, tmax] [s]
    %   - xx0: Initial state vector
    %       For Exercise 2: xx0(1) = v: Initial velocity [m/s]
    %                       xx0(2) = m: Initial mass [kg]
    %   - h: Fixed time step size [s]
    %
    % Outputs:
    %   - tt: Time vector [s]
    %   - xx: State matrix over time
    %       For Exercise 2: xx(1,:) = v: Velocity over time [m/s]
    %                       xx(2,:) = m: Mass over time [kg]
    %   - fun_eval: Number of function evaluations performed

    % Extract initial and final times from tspan
    t0   = tspan(1);   % Start time [s]
    tmax = tspan(2);   % End time [s]

    % Generate time vector with fixed step size
    tt = t0:h:tmax;
    N  = length(tt);   % Number of time points

    % Pre-allocate state matrix with NaNs for efficiency
    xx = nan(size(xx0, 1), N);  % Rows = state variables, Columns = time steps
    xx(:, 1) = xx0;             % Set initial condition

    % Initialize counter for function evaluations
    fun_eval = 0;

    % Iterate over time steps to compute the solution
    for k = 1 : N-1
        % Current state and time
        xk  = xx(:, k);   % State vector at time tk
        tk  = tt(k);      % Current time
        tk1 = tt(k+1);    % Next time

        % Compute RK2 stages (predictor-corrector)
        K1 = f(tk, xk);             % Slope at current point (predictor)
        K2 = f(tk1, xk + h*K1);     % Slope at predicted next point (corrector)

        % Update state using Heun's method: average of K1 and K2
        xx(:, k+1) = xk + h/2 * (K1 + K2);

        % Increment function evaluation counter (two evaluations per step)
        fun_eval = fun_eval + 2;
    end
end

function [tt, xx, fun_eval] = RK4(f, tspan, xx0, h)
    % RK4 Solves an ODE system using the 4th-order Runge-Kutta method.
    % This function implements a fixed-step RK4 method to numerically integrate a system of
    % ordinary differential equations. It is designed for problems like the moving object with
    % decreasing mass in Exercise 2 of the assignment, offering higher accuracy than RK2.
    %
    % Inputs:
    %   - f: Function handle for the right-hand side of the ODE system, f(t, x)
    %   - tspan: Time interval vector [t0, tmax] [s]
    %   - xx0: Initial state vector
    %       For Exercise 2: xx0(1) = v: Initial velocity [m/s]
    %                       xx0(2) = m: Initial mass [kg]
    %   - h: Fixed time step size [s]
    %
    % Outputs:
    %   - tt: Time vector [s]
    %   - xx: State matrix over time
    %       For Exercise 2: xx(1,:) = v: Velocity over time [m/s]
    %                       xx(2,:) = m: Mass over time [kg]
    %   - fun_eval: Number of function evaluations performed

    % Extract initial and final times from tspan
    t0   = tspan(1);   % Start time [s]
    tmax = tspan(2);   % End time [s]

    % Generate time vector with fixed step size
    tt = t0:h:tmax;
    N  = length(tt);   % Number of time points

    % Pre-allocate state matrix with NaNs for efficiency
    xx = nan(size(xx0, 1), N);  % Rows = state variables, Columns = time steps
    xx(:, 1) = xx0;             % Set initial condition

    % Initialize counter for function evaluations
    fun_eval = 0;

    % Iterate over time steps to compute the solution
    for k = 1 : N-1
        % Current state and time
        xk  = xx(:, k);    % State vector at time tk
        tk  = tt(k);       % Current time
        tk1 = tt(k+1);     % Next time

        % Compute RK4 stages (four slope estimates)
        K1 = f(tk, xk);              % Slope at the beginning of the interval
        K2 = f(tk + h/2, xk + h/2 * K1); % Slope at the midpoint, using K1
        K3 = f(tk + h/2, xk + h/2 * K2); % Slope at the midpoint, using K2
        K4 = f(tk1, xk + h * K3);    % Slope at the end of the interval

        % Update state using RK4: weighted average of slopes
        xx(:, k+1) = xk + h/6 * (K1 + 2*K2 + 2*K3 + K4);

        % Increment function evaluation counter (four evaluations per step)
        fun_eval = fun_eval + 4;
    end
end

function F = F_IEX4(h, alpha)
    % F_IEX4 Computes the linear operator for the Implicit Extrapolation Method of order 4 (IEX4).
    % This function constructs the discrete-time transition matrix F for a linear system using IEX4,
    % which combines four Backward Euler predictors to achieve fourth-order accuracy. It is used in
    % Exercise 5, Point 3, to simulate a stiff electrical circuit and analyze stability.
    %
    % Inputs:
    %   - h: Time step size [s]
    %   - alpha: Parameter influencing the system matrix A (not used in Exercise 5 context, but
    %            retained for generality from Exercise 3; here, it could be ignored or fixed)
    %
    % Output:
    %   - F: Linear operator (matrix) mapping x_k to x_{k+1} for the discrete-time system

    % Define the system matrix A from Exercise 5 (or Exercise 3 with alpha dependency)
    A = [0, 1; -1, 2*cos(alpha)]; % State matrix for the 2D system, e.g., electrical circuit dynamics

    % Identity matrix of the same size as A (2x2 in this case)
    I = eye(size(A));

    % Compute Backward Euler predictors for different substep sizes
    M1 = (I - h * A) \ I;       % Predictor with full step h (1 substep)
    M2 = (I - h * A / 2) \ I;   % Predictor with step h/2 (base for 2 substeps)
    M3 = (I - h * A / 3) \ I;   % Predictor with step h/3 (base for 3 substeps)
    M4 = (I - h * A / 4) \ I;   % Predictor with step h/4 (base for 4 substeps)

    % Apply substepping to achieve higher accuracy
    A1 = M1;          % Single step predictor
    A2 = M2^2;        % Two substeps: (M2)^2
    A3 = M3^3;        % Three substeps: (M3)^3
    A4 = M4^4;        % Four substeps: (M4)^4

    % Combine predictors with extrapolation coefficients for fourth-order accuracy
    F = -A1 / 6 + 4 * A2 - 27 * A3 / 2 + 32 * A4 / 3; % IEX4 operator per Exercise 5, Page 24
end

function [tt, xx, feval] = IEX4_linear(A, tspan, xx0, h)
    % IEX4_LINEAR Solves a linear ODE system using the Implicit Extrapolation Method of order 4.
    % This function implements IEX4 to numerically integrate a linear system dx/dt = A*x with a
    % fixed time step. It leverages the linearity to precompute a transition matrix, making it
    % efficient for stiff systems like the electrical circuit in Exercise 5, Point 3.
    %
    % Inputs:
    %   - A: System matrix defining the linear ODE dx/dt = A*x (e.g., 2x2 for the circuit)
    %   - tspan: Time interval vector [t0, tmax] [s]
    %   - xx0: Initial state vector
    %       For Exercise 5: xx0(1) = q: Initial capacitor charge [C]
    %                       xx0(2) = dq/dt: Initial current [A]
    %   - h: Fixed time step size [s]
    %
    % Outputs:
    %   - tt: Time vector [s]
    %   - xx: State matrix over time
    %       For Exercise 5: xx(1,:) = q: Charge over time [C]
    %                       xx(2,:) = dq/dt: Current over time [A]
    %   - feval: Number of function evaluations (approximated for matrix operations)

    % Extract initial and final times from tspan
    t0   = tspan(1);   % Start time [s]
    tmax = tspan(2);   % End time [s]

    % Generate time vector with fixed step size
    tt = t0:h:tmax;
    N  = length(tt);   % Number of time points

    % Pre-allocate state matrix with NaNs for efficiency
    xx = nan(size(xx0, 1), N);  % Rows = state variables, Columns = time steps
    xx(:, 1) = xx0;             % Set initial condition

    % Identity matrix of the same size as A (2x2 for Exercise 5)
    I = eye(size(A));

    % Compute Backward Euler predictors for different substep sizes
    M1 = (I - h * A) \ I;       % Predictor with full step h (1 substep)
    M2 = (I - h * A / 2) \ I;   % Predictor with step h/2 (base for 2 substeps)
    M3 = (I - h * A / 3) \ I;   % Predictor with step h/3 (base for 3 substeps)
    M4 = (I - h * A / 4) \ I;   % Predictor with step h/4 (base for 4 substeps)

    % Apply substepping to achieve higher accuracy
    A1 = M1;          % Single step predictor
    A2 = M2^2;        % Two substeps: (M2)^2
    A3 = M3^3;        % Three substeps: (M3)^3
    A4 = M4^4;        % Four substeps: (M4)^4

    % Compute the IEX4 transition matrix F (constant for LTI systems)
    F = -A1 / 6 + 4 * A2 - 27 * A3 / 2 + 32 * A4 / 3; % Extrapolation coefficients from Page 24

    % Initialize function evaluation counter (10 for initial matrix computations)
    feval = 10;  % Accounts for 4 inverses and 6 operations to form F

    % Iterate over time steps to propagate the solution
    for k = 1 : N-1
        % Current state
        xk = xx(:, k);

        % Update state using the precomputed transition matrix: x_{k+1} = F * x_k
        xx(:, k+1) = F * xk;

        % Increment function evaluation counter (1 matrix-vector multiply per step)
        feval = feval + 1;
    end
end

function [tt, xx, feval] = RK2_linear(A, tspan, xx0, h)
    % RK2_LINEAR Solves a linear ODE system using the 2nd-order Runge-Kutta method.
    % This function implements a fixed-step RK2 (Heun's method) to numerically integrate a linear
    % system dx/dt = A*x. It precomputes a transition matrix for efficiency, suitable for the
    % electrical circuit in Exercise 5, Point 3, to determine the largest stable step size.
    %
    % Inputs:
    %   - A: System matrix defining the linear ODE dx/dt = A*x (e.g., 2x2 for the circuit)
    %   - tspan: Time interval vector [t0, tmax] [s]
    %   - xx0: Initial state vector
    %       For Exercise 5: xx0(1) = q: Initial capacitor charge [C]
    %                       xx0(2) = dq/dt: Initial current [A]
    %   - h: Fixed time step size [s]
    %
    % Outputs:
    %   - tt: Time vector [s]
    %   - xx: State matrix over time
    %       For Exercise 5: xx(1,:) = q: Charge over time [C]
    %                       xx(2,:) = dq/dt: Current over time [A]
    %   - feval: Number of function evaluations (approximated for matrix operations)

    % Extract initial and final times from tspan
    t0   = tspan(1);   % Start time [s]
    tmax = tspan(2);   % End time [s]

    % Generate time vector with fixed step size
    tt = t0:h:tmax;
    N  = length(tt);   % Number of time points

    % Pre-allocate state matrix with NaNs for efficiency
    xx = nan(size(xx0, 1), N);  % Rows = state variables, Columns = time steps
    xx(:, 1) = xx0;             % Set initial condition

    % Identity matrix of the same size as A (2x2 for Exercise 5)
    I = eye(size(A));

    % Compute the squared system matrix for RK2 operator
    A2 = A^2;  % A^2 term for second-order accuracy

    % Compute the RK2 transition matrix F (constant for LTI systems)
    F = I + h * A + h^2 * A2 / 2;  % RK2 operator from Exercise 5, Page 23 (or Exercise 3, Page 14)

    % Initialize function evaluation counter (1 for initial matrix computation)
    feval = 1;  % Accounts for forming F (simplified estimate)

    % Iterate over time steps to propagate the solution
    for k = 1 : N-1
        % Current state
        xk = xx(:, k);

        % Update state using the precomputed transition matrix: x_{k+1} = F * x_k
        xx(:, k+1) = F * xk;

        % Increment function evaluation counter (1 matrix-vector multiply per step)
        feval = feval + 1;
    end
end

function par = par_ex_6()
    % PAR_EX_6 Defines parameters for the bouncing ball model in Exercise 6.
    % This function returns a structure containing physical and numerical parameters for the
    % bouncing ball system, including gravity, fluid density, and impact properties, as specified
    % in Exercise 6 of the assignment (Page 27).
    %
    % Inputs:
    %   - None
    %
    % Output:
    %   - par: Structure containing model parameters
    %       par.g: Gravitational acceleration [m/s^2]
    %       par.rho: Fluid density [kg/m^3]
    %       par.mb: Ball mass [kg]
    %       par.Ab: Ball cross-sectional area [m^2]
    %       par.Vb: Ball volume [m^3]
    %       par.Cd: Drag coefficient [dimensionless]
    %       par.k: Velocity attenuation factor after impact [dimensionless]
    %       par.tol: Numerical tolerance for integration [dimensionless]

    % Gravitational acceleration
    par.g = 9.81;    % Standard gravity on Earth [m/s^2]

    % Fluid density (set to 60 kg/m^3 as per Exercise 6, Point 3)
    par.rho = 60;    % Density of the fluid medium (e.g., dense air or liquid) [kg/m^3]

    % Ball mass
    par.mb = 1;      % Mass of the bouncing ball [kg]

    % Ball cross-sectional area
    par.Ab = 0.07;   % Area exposed to drag force [m^2]

    % Ball volume
    par.Vb = 0.014;  % Volume affecting buoyancy force [m^3]

    % Drag coefficient
    par.Cd = 1.17;   % Aerodynamic drag coefficient [dimensionless]

    % Velocity attenuation factor after ground impact
    par.k = 0.9;     % Coefficient reducing velocity post-bounce (0.9 means 90% retained) [dimensionless]

    % Numerical tolerance for ODE solver
    par.tol = 1e-6;  % Tolerance for convergence or event detection in integration [dimensionless]
end

function dxxdt = fun_bouncing(~, xx, par)
    % FUN_BOUNCING Computes the time derivatives for the bouncing ball system.
    % This function defines the right-hand side (RHS) of the differential equations governing the
    % dynamics of a bouncing ball under gravity, drag, and buoyancy forces, as specified in
    % Exercise 6 of the assignment (Page 27).
    %
    % Inputs:
    %   - ~: Unused time input (for compatibility with ODE solvers like ode45)
    %   - xx: State vector
    %       xx(1) = v: Velocity of the ball [m/s]
    %       xx(2) = x: Position of the ball [m]
    %   - par: Structure containing model parameters (from par_ex_6)
    %       par.g: Gravitational acceleration [m/s^2]
    %       par.rho: Fluid density [kg/m^3]
    %       par.mb: Ball mass [kg]
    %       par.Ab: Ball cross-sectional area [m^2]
    %       par.Vb: Ball volume [m^3]
    %       par.Cd: Drag coefficient [dimensionless]
    %
    % Output:
    %   - dxxdt: Time derivative of the state vector
    %       dxxdt(1) = dv/dt: Acceleration [m/s^2]
    %       dxxdt(2) = dx/dt: Velocity [m/s]

    % Extract velocity from state vector
    v = xx(1);  % Current velocity of the ball [m/s]

    % Compute drag force (proportional to v^2, direction opposite to velocity)
    Fd = -0.5 * par.rho * par.Cd * par.Ab * v .* abs(v);  % Drag force [N]

    % Compute buoyancy force (constant upward force)
    Fa = par.rho * par.Vb * par.g;  % Buoyancy force [N]

    % Compute state derivatives
    dxxdt = [ ...
        -par.g + (Fd + Fa) / par.mb;  % dv/dt: Acceleration from gravity, drag, and buoyancy [m/s^2]
        v                            % dx/dt: Velocity [m/s]
    ];
end

function [T, X] = integrate_ball(xx0, tf, par)
    % INTEGRATE_BALL Integrates the bouncing ball ODE system with event detection.
    % This function uses MATLAB's ode45 to solve the bouncing ball dynamics, restarting integration
    % after each ground impact with updated initial conditions. It is designed for Exercise 6 of the
    % assignment (Page 27) to simulate the ball's motion over a specified time interval.
    %
    % Inputs:
    %   - xx0: Initial state vector
    %       xx0(1) = v: Initial velocity [m/s]
    %       xx0(2) = x: Initial position [m]
    %   - tf: Final time for integration [s]
    %   - par: Structure containing model parameters (from par_ex_6)
    %       par.k: Velocity attenuation factor [dimensionless]
    %       par.tol: Numerical tolerance [dimensionless]
    %       (other fields used by fun_bouncing: g, rho, mb, Ab, Vb, Cd)
    %
    % Outputs:
    %   - T: Time vector across all bounces [s]
    %   - X: State matrix across all bounces
    %       X(1,:) = v: Velocity over time [m/s]
    %       X(2,:) = x: Position over time [m]

    % Set ODE solver options with high precision and event detection for ground crossing
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(x, y) x_axis_crossing(x, y));

    % Initialize integration times
    t0   = 0;      % Start time of current segment [s]
    tend = 0;      % End time of current segment [s]

    % Initialize output arrays
    T = [];        % Cumulative time vector
    X = [];        % Cumulative state matrix

    % Initialize iteration counter
    iter = 0;

    % Define the ODE function handle
    fun = @(t, x) fun_bouncing(t, x, par);

    % Initialize peak height to enter the loop
    peak = 1e3;    % Maximum height in current segment [m]

    % Integrate until final time or bounces become negligible
    while tend < tf && peak > 1e-3
        % Solve ODE from current t0 to tf with event detection
        [tt, xx, ~, ~, event] = ode45(fun, [t0, tf], xx0, options);

        % Transpose for consistency (rows = states, columns = time)
        xx = xx.';  % xx(1,:) = velocity, xx(2,:) = position

        % Update initial conditions after bounce: reverse velocity with attenuation, position = 0
        xx0 = [-par.k * xx(1, end); 0];  % v^+ = -k * v^-, x = 0

        % Transpose time vector
        tt = tt.';

        % Update end time and peak height
        tend = tt(end);         % End time of this segment [s]
        peak = max(xx(2, :));   % Maximum height reached [m]

        % If an event (ground impact) occurred, set final velocity to NaN and position to 0
        if event
            xx(:, end) = [nan; 0];  % Mark bounce point for clarity in plots
        end

        % Append results to cumulative arrays
        T = [T, tt];  % Concatenate time points
        X = [X, xx];  % Concatenate state vectors

        % Set new start time for next segment
        t0 = tend;

        % Increment iteration counter
        iter = iter + 1;
    end
end

function [value, isterminal, direction] = x_axis_crossing(~, xx)
    % X_AXIS_CROSSING Event function to detect ground impact in the bouncing ball model.
    % This function is used with MATLAB's ODE solvers (e.g., ode45) to detect when the ball's
    % position crosses the ground (x = 0) from above, triggering a bounce. It is part of the
    % integration process in Exercise 6 of the assignment (Page 28).
    %
    % Inputs:
    %   - ~: Unused time input (for compatibility with ODE solver event functions)
    %   - xx: State vector
    %       xx(1) = v: Velocity of the ball [m/s]
    %       xx(2) = x: Position of the ball [m]
    %
    % Outputs:
    %   - value: Value whose zero crossing triggers the event (position x)
    %   - isterminal: Flag to stop integration (1 = stop, 0 = continue)
    %   - direction: Direction of zero crossing to detect (-1 = from positive to negative)

    % Define the event condition: ball hits the ground when position x = 0
    value = xx(2);  % Position of the ball [m]; event triggers when this equals 0

    % Set integration to terminate at the event
    isterminal = 1;  % 1 means stop integration when the event occurs (ground impact)

    % Detect only downward crossings (from positive to negative position)
    direction = -1;  % -1 means event triggers when x decreases through 0 (ball falling)
end