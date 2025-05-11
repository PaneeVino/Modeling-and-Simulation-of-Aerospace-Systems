clear; close all; clc; 
plotStyle();
Colors = [
    0, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330  % Light Blue
    ];

% This section performs a simulation of the satellite's thermal control system
% using causal modeling and a proportional control (kp).

% Set the total simulation time
tf = 50 * 3600; % Total time in seconds (50 hours)

% Set proportional gain for the control system
kp = 6e-5; % Value of kp (proportional gain for temperature control)

% Define the initial state vector
%   - Temperatures (T1, T2, T3, T4, T5) initialized to 298.15 K
%   - Radiator angle (theta) initialized to -0.4*pi (fully closed)
%   - Radiator angular velocity (omega) initialized to 0 rad/s
%   - Motor current (i) initialized to 0 A
x0 = [298.15 * ones(5, 1); -0.4 * pi; 0; 0];

% options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
% odefun = @(t, x) rhs(t, x, kp);
% [t, x] = ode15s(odefun, [0, tf], x0, options);
% x=x.';
[t_seconds, x, Ie, Te, Xe] = integrate_system(x0, tf, kp);

% Extract specific state variables from the simulation results
T1 = x(1, :); % Temperature of the main body
T2 = x(2, :); % Temperature of solar panel 1
T3 = x(3, :); % Temperature of solar panel 2
T4 = x(4, :); % Temperature of radiator 1
T5 = x(5, :); % Temperature of radiator 2
theta = x(6, :); % Radiator angle [rad]
omega = x(7, :); % Radiator angular velocity [rad/s]
i = x(8, :); % Motor current [A]
t = t_seconds / 3600; %Time in hours
K=constant();

% Plot main body temperature (T1) over time
figure
yline(K.T1ref, '--', 'Color', Colors(1, :), 'LineWidth', 2,  'DisplayName', '$T_{ref}=294.15 \ K$'); % Reference temperature
hold on
box on
plot(t, T1, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', '$T_1$'); % Plot T1 in blue (time converted to hours)
% yline(K.Tmin, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$T_{min}$'); % Minimum temperature threshold
% yline(K.Tmax, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$T_{max}$'); % Maximum temperature threshold
yline(K.T1ref - 1e-3 * K.T1ref, '--', 'Color', Colors(5, :), 'LineWidth', 2, 'DisplayName', 'Boundaries after 10 h'); % T1ref lower bound
yline(K.T1ref + 1e-3 * K.T1ref, '--', 'Color', Colors(5, :), 'LineWidth', 2, 'HandleVisibility','off'); % T1ref upper bound
xline(10, '--', 'Color', 'k', 'LineWidth', .5, 'HandleVisibility','off'); % Vertical line at t = 10 hours
legend('Location','northeast', 'Orientation', 'vertical', 'FontSize', 20); % Legend with LaTeX formatting
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('T [K]', 'FontSize', 20); % Label for y-axis (temperature in Kelvin)

% Plot radiator angle (theta) over time
figure
plot(t, rad2deg(theta), 'Color', Colors(1, :), 'LineWidth', 2.5,  'DisplayName', '$\theta$'); % Radiator angle in degrees
hold on
box on
yline(rad2deg(K.theta_min), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\theta$ boundaries')
yline(rad2deg(K.theta_max), '--', 'LineWidth', 2, 'Color', 'k', 'HandleVisibility','off')
yline(rad2deg(theta(end)), '--', 'LineWidth', 2, 'Color', Colors(2, :), 'DisplayName', '$\overline{\theta} \approx -36.0525$ deg')
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('$\theta$ [deg]', 'FontSize', 20); % Label for y-axis (angle in degrees)
legend('Location','best', 'Orientation','vertical', 'FontSize',20)

% Plot radiator angular velocity (omega) over time
figure
plot(t, rad2deg(omega), 'Color', Colors(1, :), 'LineWidth', 2.5,  'DisplayName', '$\omega$'); % Angular velocity in degrees per second
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('$\omega$ [deg/s]', 'FontSize', 20); % Label for y-axis (angular velocity in deg/s)

% Plot temperatures of the solar panels  over time
figure
plot(t, T2, 'Color', Colors(1, :), 'LineWidth', 2.5,  'DisplayName', '$T_2$'); % Temperature of solar panels
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('T [K]', 'FontSize', 20);

% Plot temperatures of the radiators  over time
figure
plot(t, T4, 'Color', Colors(2, :), 'LineWidth', 2.5,  'DisplayName', '$T_2$'); % Temperature of radiators
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('T [K]', 'FontSize', 20);

% Plot current over time
figure
plot(t, i, 'Color', Colors(1, :), 'LineWidth', 2.5,  'DisplayName', '$T_2$'); % Current in the circuit
xlabel('time [h]', 'FontSize', 20); % Label for x-axis (time in hours)
ylabel('i [A]', 'FontSize', 20);
ylim([-1e-8, 0.5e-8])
%%
Data_tab = readtable("exportedVariables.csv", 'VariableNamingRule', 'preserve');

t_dymola     = Data_tab.time.';
omega_dymola = rad2deg(Data_tab.("electro_Mechanical_Model.inertia.w")*1e-3).';
theta_dymola = rad2deg(Data_tab.("electro_Mechanical_Model.angleSensor.phi")).';
T1_dymola    = Data_tab.("thermal_model.MainBody.T").';
T2_dymola    = Data_tab.("thermal_model.Panel1.T").';
T4_dymola    = Data_tab.("thermal_model.Radiator1.T").';
i_dymola     = Data_tab.("electro_Mechanical_Model.resistor.i")*1e-9.';

linewidth_dymola = 3;
linewidth_matlab = 5;
lgd_fontsize = 28;
axis_lab_fontsize = 28;
tick_fontsize = 28;  % Define tick label font size

figure
tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact')

% First plot
nexttile
plot(t_dymola, T1_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, T1, '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
legend('Location','best', 'Orientation', 'vertical', 'FontSize', lgd_fontsize);
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('$T_1$ [K]', 'FontSize', axis_lab_fontsize);
ax = gca;
ax.FontSize = tick_fontsize;  % Set the font size of the ticks

% Second plot
nexttile
plot(t_dymola, T2_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, T2, '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
legend('Location','best', 'Orientation', 'vertical', 'FontSize', lgd_fontsize);
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('$T_2$ [K]', 'FontSize', axis_lab_fontsize);
ax = gca;
ax.FontSize = tick_fontsize;

% Third plot
nexttile
plot(t_dymola, T4_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, T4, '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
legend('Location','best', 'Orientation', 'vertical', 'FontSize', lgd_fontsize);
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('$T_4$ [K]', 'FontSize', axis_lab_fontsize);
ax = gca;
ax.FontSize = tick_fontsize;

% Fourth plot
nexttile
plot(t_dymola, theta_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, rad2deg(theta), '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
legend('Location','best', 'Orientation', 'vertical', 'FontSize', lgd_fontsize);
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('$\theta$ [deg]', 'FontSize', axis_lab_fontsize);
ax = gca;
ax.FontSize = tick_fontsize;

% Fifth plot
nexttile
plot(t_dymola, omega_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, rad2deg(omega), '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('$\omega$ [deg/s]', 'FontSize', axis_lab_fontsize);
legend('Location','best', 'Orientation','vertical', 'FontSize', lgd_fontsize);
ax = gca;
ax.FontSize = tick_fontsize;

% Sixth plot
nexttile
plot(t_dymola, i_dymola, 'Color', Colors(1, :), 'LineWidth', linewidth_dymola, 'DisplayName', 'Dymola' )
hold on
plot(t, i, '--', 'Color', Colors(2, :), 'LineWidth', linewidth_matlab, 'DisplayName', 'MATLAB')
box on
xlabel('time [h]', 'FontSize', axis_lab_fontsize);
ylabel('i [A]', 'FontSize', axis_lab_fontsize);
legend('Location','best', 'Orientation','vertical', 'FontSize', lgd_fontsize);
ylim([-1e-8, 0.5e-8]);
ax = gca;
ax.FontSize = tick_fontsize;
%%
% Remove duplicate time points from t_dymola
[t_unique, idx_unique] = unique(t_dymola);  % keep first occurrence of each time
theta_unique = theta_dymola(idx_unique);
omega_unique = omega_dymola(idx_unique);
i_unique     = i_dymola(idx_unique);
T1_unique    = T1_dymola(idx_unique);
T2_unique    = T2_dymola(idx_unique);
T4_unique    = T4_dymola(idx_unique);

% Now perform interpolation at desired time vector `t`
theta_dym = interp1(t_unique, theta_unique, t);
omega_dym = interp1(t_unique, omega_unique, t);
i_dym     = interp1(t_unique, i_unique, t);
T1_dym    = interp1(t_unique, T1_unique, t);
T2_dym    = interp1(t_unique, T2_unique, t);
T4_dym    = interp1(t_unique, T4_unique, t);

err_theta = abs(theta_dym - rad2deg(theta))/max(abs(rad2deg(theta)));
err_omega = abs(omega_dym - rad2deg(omega))/max(abs(rad2deg(omega)));
err_i     = abs(i_dym - i)/max(abs(i));
err_T1    = abs(T1_dym - T1)/max(abs(T1));
err_T2    = abs(T2_dym - T2)/max(abs(T2));
err_T4    = abs(T4_dym - T4)/max(abs(T4));

Err       = [err_theta;
             err_omega;
             err_i; 
             err_T1;
             err_T2;
             err_T4];
err = vecnorm(Err);

linewidth_err = 2.5;
figure
semilogy(t, err_theta, 'LineWidth', linewidth_err, 'Color', Colors(1, :), 'DisplayName', '$\varepsilon_{\theta}$')
hold on
grid minor
semilogy(t, err_omega, 'LineWidth', linewidth_err, 'Color', Colors(2, :), 'DisplayName', '$\varepsilon_{\omega}$')
semilogy(t, err_i, 'LineWidth', linewidth_err, 'Color', Colors(3, :), 'DisplayName', '$\varepsilon_{i}$')
semilogy(t, err_T1, 'LineWidth', linewidth_err, 'Color', Colors(4, :), 'DisplayName', '$\varepsilon_{T_1}$')
semilogy(t, err_T2, 'LineWidth', linewidth_err, 'Color', Colors(5, :), 'DisplayName', '$\varepsilon_{T_2}$')
semilogy(t, err_T4, 'LineWidth', linewidth_err, 'Color', Colors(6, :), 'DisplayName', '$\varepsilon_{T_4}$')
xlabel('time [h]', 'FontSize', 20);
ylabel('Relative Error [-]', 'FontSize', 20); 
legend('Location','best', 'Orientation','vertical', 'FontSize', 20, 'NumColumns', 2)
ylim([1e-15, 1e0])

figure
semilogy(t, err, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Relative Error [-]')
yline(err(end), '--', 'Color', 'k', 'LineWidth', 1, 'DisplayName', '$\overline{\varepsilon} \approx 4.04 \cdot 10^{-4}$')
xlabel('time [h]', 'FontSize', 20);
ylabel('Relative Error [-]', 'FontSize', 20); 
legend('Location','best', 'Orientation','vertical', 'FontSize', 20)
%%
function plotStyle()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set figure properties for better looking plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpreter:
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
% lines:
set(0,'defaultLineLineWidth', 1.5);
set(0,'defaultLineMarkerSize',6) ;
set(0,'defaultLineMarkerEdgeColor','k')
set(0,'defaultLineMarkerFaceColor','auto')
% legend:
set(0, 'defaultLegendLocation','northoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end

function K = constant()

    % CONSTANT Define the constants for the thermal and dynamic model.
    % This function creates and returns a structure `K` containing all the
    % constants used in the thermal control and dynamic system of the satellite.

    % Electrical constants for the motor control
    R = 0.1;       % Resistance of the motor winding [Ohm]
    L = 1e-3;      % Inductance of the motor winding [Henry]
    km = 0.3;      % Motor constant, torque per unit current [Nm/A]

    % Mechanical properties of the radiator
    mr = 0.2;      % Mass of the radiator [kg]
    Lr = 0.5;      % Length of the radiator [m]

    % Geometric dimensions of the satellite body and panels
    L1 = 1.5;      % Height of the main satellite body [m]
    L2 = 0.5;      % Width of the main satellite body [m]
    L3 = L2;       % Depth of the main satellite body [m] (assumed equal to width)
    Lp = 0.95;     % Length of the solar panels [m]

    % Thermal properties
    Ps = 1350;     % Solar flux [W/m^2]
    C1 = 1.5e5;    % Heat capacity of the main satellite body [J/K]
    C2 = 1187.5;   % Heat capacity of solar panel 1 [J/K]
    C3 = C2;       % Heat capacity of solar panel 2 [J/K]
    C4 = 30;       % Heat capacity of radiator 1 [J/K]
    C5 = C4;       % Heat capacity of radiator 2 [J/K]

    % Thermal conductances between nodes
    G12 = 10;      % Conductance between main body and solar panel 1 [W/K]
    G13 = G12;     % Conductance between main body and solar panel 2 [W/K]
    G14 = G13;     % Conductance between main body and radiator 1 [W/K]
    G15 = G14;     % Conductance between main body and radiator 2 [W/K]

    % Surface properties (absorptivity and emissivity)
    aLpha1 = 0.6;  % Absorptivity of the main body surface [-]
    aLpha2 = 0.78; % Absorptivity of solar panel 1 [-]
    aLpha3 = aLpha2; % Absorptivity of solar panel 2 [-]
    eps1 = 0.45;   % Emissivity of the main body surface [-]
    eps2 = 0.75;   % Emissivity of solar panel 1 [-]
    eps3 = eps2;   % Emissivity of solar panel 2 [-]

    % Radiator-specific properties
    K.epsmin = 0.01; % Minimum emissivity of the radiator (closed state) [-]
    K.epsmax = 0.98; % Maximum emissivity of the radiator (open state) [-]

    % Boltzmann constant and external conditions
    K.sigma = 5.67e-8; % Stefan-Boltzmann constant [W/m^2/K^4]
    K.Tds = 3;         % Deep space temperature [K]

    % Reference and boundary temperatures
    K.T1ref = 294.15;  % Reference temperature of the main body [K]
    K.Tmin = 290;      % Minimum allowable temperature [K]
    K.Tmax = 300;      % Maximum allowable temperature [K]

    % Radiator angle constraints
    K.theta_min = -0.4 * pi; % Minimum radiator angle (closed state) [rad]
    K.theta_max = 0;         % Maximum radiator angle (open state) [rad]

    % Package all constants into the output structure
    K.R = R;
    K.L = L;
    K.km = km;
    K.mr = mr;
    K.Lr = Lr;
    K.L1 = L1;
    K.L2 = L2;
    K.L3 = L3;
    K.LP = Lp;
    K.Ps = Ps;
    K.C1 = C1;
    K.C2 = C2;
    K.C3 = C3;
    K.C4 = C4;
    K.C5 = C5;
    K.G12 = G12;
    K.G13 = G13;
    K.G14 = G14;
    K.G15 = G15;
    K.aLpha1 = aLpha1;
    K.aLpha2 = aLpha2;
    K.aLpha3 = aLpha3;
    K.eps1 = eps1;
    K.eps2 = eps2;
    K.eps3 = eps3;

end

function xdot = rhs(~, x, kp)
    % RHS Compute the time derivatives of the system's state variables.
    % This function defines the right-hand side (RHS) of the differential equations
    % governing the satellite's thermal control and radiator dynamics.
    %
    % Inputs:
    %   - ~: Unused time input (for compatibility with ODE solvers)
    %   - x: State vector
    %       x(1) = T1: Temperature of the main body [K]
    %       x(2) = T2: Temperature of solar panel 1 [K]
    %       x(3) = T3: Temperature of solar panel 2 [K]
    %       x(4) = T4: Temperature of radiator 1 [K]
    %       x(5) = T5: Temperature of radiator 2 [K]
    %       x(6) = th: Radiator angle [rad]
    %       x(7) = om: Angular velocity of the radiator [rad/s]
    %       x(8) = i: Current in the motor circuit [A]
    %   - kp: Proportional control gain
    %
    % Output:
    %   - xdot: Time derivative of the state vector

    % Pre-allocate xdot as an 8-element column vector
    xdot = nan(8,1);

    % Retrieve constants from the constant function
    K = constant();

    % Extract state variables for clarity
    T1 = x(1); T2 = x(2); T3 = x(3); T4 = x(4); T5 = x(5);
    th = x(6); om = x(7); i = x(8);

    % % %Enforce radiator angle limits and handle boundary conditions
    % 
    % if (th < K.theta_min && T1 < K.T1ref) || (th > K.theta_max && T1 > K.T1ref)
    % 
    %     kp = 0;
    % 
    % end

    % Compute total area of the main body exposed to radiation
    A1tot = 2 * K.L2 * K.L3 + 4 * K.L1 * K.L2;

    % Calculate radiator emissivity as a function of its angle
    eps4 = K.epsmin + (K.epsmax - K.epsmin) / (0.4 * pi) * (th + 0.4 * pi);
    eps5 = eps4; % Radiators are assumed to have identical emissivity

    % Calculate the moment of inertia for the radiator
    Jr = (1/3) * K.mr * K.Lr^2; % Simplified rotational inertia

    % Thermal dynamics equations for each node

    % Main body temperature (Node 1)
    xdot(1) = 1 / K.C1 * ( ...
        K.aLpha1 * K.Ps * K.L2 * K.L3 ... % Absorbed solar flux
        - K.eps1 * A1tot * K.sigma * (T1^4 - K.Tds^4) ... % Radiative losses
        + K.G12 * (T2 - T1) + K.G13 * (T3 - T1) ... % Conductive heat transfer
        + K.G14 * (T4 - T1) + K.G15 * (T5 - T1)); % From radiators

    % Solar panel 1 temperature (Node 2)
    xdot(2) = 1 / K.C2 * ( ...
        K.aLpha2 * K.Ps * K.LP * K.L2 ... % Absorbed solar flux
        - K.eps2 * K.LP * K.L2 * K.sigma * (T2^4 - K.Tds^4) ... % Radiative losses
        - K.G12 * (T2 - T1)); % Conductive heat transfer to the main body

    % Solar panel 2 temperature (Node 3)
    xdot(3) = 1 / K.C3 * ( ...
        K.aLpha3 * K.Ps * K.LP * K.L3 ... % Absorbed solar flux
        - K.eps3 * K.LP * K.L3 * K.sigma * (T3^4 - K.Tds^4) ... % Radiative losses
        - K.G13 * (T3 - T1)); % Conductive heat transfer to the main body

    % Radiator 1 temperature (Node 4)
    xdot(4) = 1 / K.C4 * ( ...
        - eps4 * K.Lr * K.L2 * K.sigma * (T4^4 - K.Tds^4) ... % Radiative losses
        - K.G14 * (T4 - T1)); % Conductive heat transfer to the main body

    % Radiator 2 temperature (Node 5)
    xdot(5) = 1 / K.C5 * ( ...
        - eps5 * K.Lr * K.L2 * K.sigma * (T5^4 - K.Tds^4) ... % Radiative losses
        - K.G15 * (T5 - T1)); % Conductive heat transfer to the main body

    % Radiator angular position (theta, Node 6)
    xdot(6) = om; % Angular velocity (derivative of position)

    % Radiator angular velocity (omega, Node 7)
    xdot(7) = K.km / Jr * i; % Torque induced by motor current

    % Motor current (Node 8)
    xdot(8) = 1 / K.L * kp * (T1 - K.T1ref) ... % Control input proportional to temperature error
            - K.R / K.L * i ... % Ohmic losses
            - K.km / K.L * om; % Back-emf from radiator motion

end

function [T, X, Ie, Te, Xe] = integrate_system(xx0, tf, Kp)
    % INTEGRATE_SYSTEM Numerically integrates the satellite thermal control system ODEs.
    % This function solves the non-linear ODEs defined in rhs() over a specified time
    % span, handling events such as temperature or angle limits, and stitches together
    % solution segments when events occur.
    %
    % Inputs:
    %   - xx0: Initial state vector (8x1) containing:
    %          xx0(1) = T1: Initial main body temperature [K]
    %          xx0(2) = T2: Initial solar panel 1 temperature [K]
    %          xx0(3) = T3: Initial solar panel 2 temperature [K]
    %          xx0(4) = T4: Initial radiator 1 temperature [K]
    %          xx0(5) = T5: Initial radiator 2 temperature [K]
    %          xx0(6) = th: Initial radiator angle [rad]
    %          xx0(7) = om: Initial angular velocity [rad/s]
    %          xx0(8) = i: Initial motor current [A]
    %   - tf: Final simulation time [s]
    %   - Kp: Proportional control gain [V/K] for the motor control law
    %
    % Outputs:
    %   - T: Time vector of the full simulation [s]
    %   - X: State matrix (8xN) of the solution at each time step
    %   - Ie: Event indices (not fully returned in this implementation)
    %   - Te: Event times (not fully returned in this implementation)
    %   - Xe: Event states (not fully returned in this implementation)
    %
    % Note: Ie, Te, Xe are initialized but not populated due to missing event storage.

    % Load system constants from the constant() function
    K = constant();

    % Initialize output arrays as empty
    T = [];  % Time points
    X = [];  % State trajectories
    Ie = []; % Event indices (unused in output)
    Te = []; % Event times (unused in output)
    Xe = []; % Event states (unused in output)

    % Set ODE solver options with high precision and event detection
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t, x) EventFunction(t, x));

    % Initialize starting time
    t0 = 0;

    % Define the ODE function with the given proportional gain Kp
    odefun = @(t, x) rhs(t, x, Kp);

    % Integrate the system in segments, handling events until final time is reached
    while t0 < tf
        % Solve ODEs using ode15s (stiff solver) from t0 to tf with event detection
        [tt, xx, ~, ~, ie] = ode15s(odefun, [t0, tf], xx0, options);

        % Transpose outputs for consistency (rows to columns)
        xx = xx.'; % State matrix: 8xN
        tt = tt.'; % Time vector: 1xN

        % Handle events based on event indices (ie) from EventFunction
        if any(ie == 1) && xx(1, end) - K.T1ref <= 0
            % Event 1: T1 below T1ref and radiator at minimum angle
            xx(6, end) = K.theta_min; % Force radiator to closed position
            xx(7, end) = 0;           % Stop angular velocity
            %xx(8, end) = 0;
            odefun = @(t, x) rhs(t, x, 0); % Disable control (Kp = 0)

        elseif any(ie == 2) && xx(1, end) - K.T1ref >= 0
            % Event 2: T1 above T1ref and radiator at maximum angle
            xx(6, end) = K.theta_max; % Force radiator to fully open position
            xx(7, end) = 0;           % Stop angular velocity
            %xx(8, end) = 0;
            odefun = @(t, x) rhs(t, x, 0); % Disable control (Kp = 0)

        elseif any(ie == 3) && xx(6, end) - K.theta_min >= 0
            % Event 3: Radiator angle exceeds minimum (re-enable control if needed)
            odefun = @(t, x) rhs(t, x, Kp); % Restore original control gain

        elseif any(ie == 4) && xx(6, end) - K.theta_max <= 0
            % Event 4: Radiator angle below maximum (re-enable control if needed)
            odefun = @(t, x) rhs(t, x, Kp); % Restore original control gain
        end

        % Update initial conditions for the next segment
        xx0 = xx(:, end); % Last state becomes new initial condition
        t0 = tt(end);     % Last time becomes new starting time

        % Append current segment results to output arrays
        T = [T, tt]; % Concatenate time points
        X = [X, xx]; % Concatenate state trajectories
    end
end

function [value, isterminal, direction] = EventFunction(~, x)
    % EVENTFUNCTION Defines events for stopping ODE integration.
    % This function specifies conditions under which the ODE solver (e.g., ode15s)
    % should halt integration, such as when the radiator angle or main body temperature
    % reaches specific limits. It is used with the 'Events' option in odeset.
    %
    % Inputs:
    %   - ~: Time input (unused, included for compatibility with ODE solvers)
    %   - x: State vector (8x1) containing:
    %        x(1) = T1: Main body temperature [K]
    %        x(6) = th: Radiator angle [rad]
    %        (Other states like T2, T3, etc., are not used here)
    %
    % Outputs:
    %   - value: Vector of event functions; integration stops when any value = 0
    %   - isterminal: Vector indicating whether each event stops integration (1 = yes)
    %   - direction: Vector specifying the direction of zero-crossing to detect:
    %                -1 = decreasing, 1 = increasing, 0 = either
    %
    % Events:
    %   1. Radiator angle reaches minimum (theta_min, closed position)
    %   2. Radiator angle reaches maximum (theta_max, fully open)
    %   3. Main body temperature exceeds reference (T1 > T1ref)
    %   4. Main body temperature drops below reference (T1 < T1ref)

    % Load system constants from the constant() function
    K = constant();

    % Define event conditions (zero-crossing triggers)
    value = [x(6) - K.theta_min;          % Event 1: theta reaches minimum (closed)
             x(6) - K.theta_max;          % Event 2: theta reaches maximum (open)
             x(1) - K.T1ref;              % Event 3: T1 exceeds T1ref
             x(1) - K.T1ref];             % Event 4: T1 drops below T1ref

    % Specify that all events terminate integration
    isterminal = [1; 1; 1; 1]; % 1 = stop integration when event occurs

    % Specify direction of zero-crossing for each event
    direction = [-1;               % Event 1: Detect when theta decreases to theta_min
                  1;               % Event 2: Detect when theta increases to theta_max
                  1;               % Event 3: Detect when T1 increases past T1ref
                 -1];              % Event 4: Detect when T1 decreases past T1ref
end
