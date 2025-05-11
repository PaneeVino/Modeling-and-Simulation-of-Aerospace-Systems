% Modeling and Simulation of Aerospace systems (2024/2025)
% Assignment # 1
% Author: Marcello Pareschi 
%% Exercise 2 - Causal Modeling 
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

K=constant();

rho = K.P_tank / (K.R * K.T_tank);
C   = sqrt(K.k * rho * K.P_tank * (2 / (K.k + 1))^((K.k + 1) / (K.k - 1)));

xv0 = K.A0/K.xv_max +K.xv_max -1.2*1e-3/(K.xv_max*C*K.ve); %m
dLdx = -K.beta  / (K.alpha + K.beta * xv0)^2;
i0 = -sqrt(2*K.kv*xv0/dLdx);
x0 = [0; 0; xv0; 0; i0];
var = [K.ba_min; K.ka_min; K.Rin_min; K.Rf_max];

options = odeset('RelTol',1e-10, 'AbsTol',1e-12);
odefun = @(t, x) rhs(t, x, var);
T = 2*pi/K.om0;
tf = 3 * T;
tspan = [0, tf];
[t,x] = ode15s(odefun,tspan,x0,options);
[Th, D] = get_forces(t, x);

figure
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Corresponding data and labels
y_data = {x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), Th};
y_labels = {'$x_a$ [m]', '$v_a$ [m/s]', '$x_v$ [m]', ...
            '$v_v$ [m/s]', '$I$ [A]', '$T$ [N]'};

linewidth = 4;
tick_fontsize = 28;

for i = 1:6
    % Create the tile and get its axes
    main_ax = nexttile;

    % Main plot
    plot(t/3600, y_data{i}, 'Color', Colors(1, :), 'LineWidth', linewidth)
    xlabel('time [h]')
    ylabel(y_labels{i}, 'Interpreter', 'latex')
    main_ax.FontSize = tick_fontsize;

end


figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(t/T, Th*1e3, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Thrust')
hold on
plot(t/T, D*1e3, '--', 'Color', Colors(2, :), 'LineWidth', 2.5, 'LineWidth', 2.5, 'DisplayName', 'Drag')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Force [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

err = abs(Th - D);
err_max = max(err);
nexttile
semilogy(t/T, err*1e3, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', '$\varepsilon$')
hold on
yline(err_max*1e3, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\varepsilon_{max} \approx 3.7029 \cdot 10^{-3}$')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Abs. Err. [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 20)
ylim([1e-6, 1e-2])
yticks([1e-6, 1e-4, 1e-2])

disp(err_max)
%%
x0 = [0; 0; K.xv_max; 0; 0];

[t,x] = ode15s(odefun,tspan,x0,options);
[Th, D] = get_forces(t, x);
indx = t <= 30;

linewidth = 4;
tick_fontsize = 28;


figure
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Limits for the inset axes
inset_limits = {
    [0, 2.5e-8], ...
    [-1e-10, 0.5e-11], ...
    [0.3e-5, 1e-5], ...
    [-3e-8, 1e-8], ...
    [-4e-4, 0], ...
    [0, 7e-3]
};

% Corresponding data and labels
y_data = {x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), Th};
y_labels = {'$x_a$ [m]', '$v_a$ [m/s]', '$x_v$ [m]', ...
            '$v_v$ [m/s]', '$I$ [A]', '$T$ [N]'};

% Custom inset positions: [x, y, width, height]
inset_positions = {
    [0.13, 0.65, 0.15, 0.17],   % top-left of first tile
    [0.45, 0.66, 0.15, 0.17],
    [0.78, 0.63, 0.15, 0.17],
    [0.12, 0.13, 0.15, 0.17],
    [0.45, 0.23, 0.15, 0.17],
    [0.78, 0.25, 0.15, 0.17]
};

for i = 1:6
    % Create the tile and get its axes
    main_ax = nexttile;

    % Main plot
    plot(t/3600, y_data{i}, 'Color', Colors(1, :), 'LineWidth', linewidth)
    xlabel('time [h]')
    ylabel(y_labels{i}, 'Interpreter', 'latex')
    main_ax.FontSize = tick_fontsize;
    ylim(inset_limits{i})

    % Create inset using custom position
    inset_ax = axes('Position', inset_positions{i});
    inset_ax.FontSize = tick_fontsize;
    box on; hold on;
    plot(t(indx)/3600, y_data{i}(indx), 'Color', Colors(1, :), 'LineWidth', linewidth)
    set(inset_ax, 'FontSize', 8)

    % Return focus to main axes
    axes(main_ax);
end

figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(t/T, Th*1e3, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Thrust')
hold on
plot(t/T, D*1e3, '--', 'Color', Colors(2, :), 'LineWidth', 2.5, 'LineWidth', 2.5, 'DisplayName', 'Drag')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Force [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

err = abs(Th - D);
err_max = max(err(400:end));
nexttile
semilogy(t/T, err*1e3, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', '$\varepsilon$')
hold on
yline(err_max*1e3, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\varepsilon_{max} \approx 3.7034 \cdot 10^{-6}$')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Abs. Err. [N]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 20)
ylim([1e-5, 1e-1])
yticks([1e-5, 1e-3, 1e-1])
disp(err_max)
%%
modelName = 'Ass2_Ex2_Part2_sim';
load_system(modelName);

K = constant();
T = 2*pi/K.om0;

set_param(modelName, ... 
    'StartTime', '0', ...
    'StopTime', num2str(3*T), ...  
    'SolverType', 'Variable-step', ...
    'Solver', 'ode15s', ...
    'MaxStep', '1e-2', ...
    'RelTol', '1e-10', ...
    'ABsTol', '1e-12');

simOut = sim(modelName);

t   = simOut.tout;
D   = simOut.D;
Th  = simOut.T*1e3;
xa  = simOut.xa;
va  = simOut.va;
xv  = simOut.xv;
vv  = simOut.vv;
I   = simOut.I;

indx = t >= 30;

linewidth = 4;
tick_fontsize = 28;


figure
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Limits for the inset axes
inset_limits = {
    [0, inf], ...
    [-3e-12, 1e-11], ...
    [0.7e-5, 1e-5], ...
    [-2.5e-9, 2.5e-9], ...
    [-6e-4, 1e-4], ...
    [-5e-3, 15e-3]
};

y_data = {xa, va, xv, vv, I, Th-D};
y_labels = {'$x_a$ [m]', '$v_a$ [m/s]', '$x_v$ [m]', ...
            '$v_v$ [m/s]', '$I$ [A]', '$T - D$ [mN]'};

inset_positions = {
    [0.13, 0.65, 0.15, 0.17],  
    [0.45, 0.66, 0.15, 0.17],
    [0.78, 0.63, 0.15, 0.17],
    [0.12, 0.13, 0.15, 0.17],
    [0.45, 0.23, 0.15, 0.17],
    [0.78, 0.25, 0.15, 0.17]
};

for i = 1:6
    main_ax = nexttile;

    plot(t/3600, y_data{i}, 'Color', Colors(1, :), 'LineWidth', linewidth)
    xlabel('time [h]')
    ylabel(y_labels{i}, 'Interpreter', 'latex')
    main_ax.FontSize = tick_fontsize;
    ylim(inset_limits{i})

end

figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(t/T, Th, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Thrust')
hold on
plot(t/T, D, '--', 'Color', Colors(2, :), 'LineWidth', 2.5, 'LineWidth', 2.5, 'DisplayName', 'Drag')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Force [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

err = abs(Th - D);
err_max = max(err(indx));
nexttile
semilogy(t/T, err, 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', '$\varepsilon$')
hold on
yline(err_max, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\varepsilon_{max} \approx 0.0115$')
xlabel('Orbital Periods [-]', 'FontSize', 20)
ylabel('Abs. Err. [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 20)
% ylim([1e-8, 1e-4])
% yticks([1e-8, 1e-6, 1e-4])
%%
x0 = [0; 0; 0.75*K.xv_max; 0; 0];

odefun = @(t, x) rhs(t, x, var);
T = 2*pi/K.om0;
tf = 3 * T;
tspan = [0, tf];
[t,x] = ode15s(odefun,t,x0,options);
[Th, D_matlab] = get_forces(t, x);

yy_data = {x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), Th};

err_matlab = yy_data{6} - D_matlab;
e_max_matlab = max(err_matlab(indx))*1e3;
disp(e_max_matlab)
err_simscape = y_data{6};
e_max_simscape = max(err_simscape(indx));
disp(e_max_simscape)
T_matlab = yy_data{6};
T_simscape = y_data{6} + D;
%%
figure
tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')
nexttile
plot(t(indx)/3600, T_simscape(indx), 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Simscape')
hold on
plot(t(indx)/3600, T_matlab(indx)*1e3, '--', 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'Matlab')
xlabel('Time [h]', 'FontSize', 20)
ylabel('Thrust [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'vertical', 'FontSize', 20)

nexttile
plot(t(indx)/3600, err_simscape(indx), 'Color', Colors(1, :), 'LineWidth', 2.5, 'DisplayName', 'Simscape')
hold on
yline(e_max_simscape, '--', ...
    'Color', Colors(1, :), ...
    'LineWidth', 2.5, ...
    'DisplayName', '0.0115', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left')
plot(t(indx)/3600, err_matlab(indx)*1e3, 'Color', Colors(2, :), 'LineWidth', 2.5, 'DisplayName', 'Matlab')
yline(e_max_matlab, '--', ...
    'Color', Colors(2, :), ...
    'LineWidth', 2.5, ...
    'DisplayName', '0.0036', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left')
xlabel('Time [h]', 'FontSize', 20)
ylabel('Abs. Err. [mN]', 'FontSize', 20)
legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 20, 'NumColumns', 2)
%%

% function K = constant()
% % This function defines and returns a struct `K` containing all the 
% % constants and parameter ranges used in the ACS simulation.
% % Each parameter is named based on its physical significance.
% 
% % Satellite and accelerometer properties
% K.Msc = 300; % Spacecraft mass (kg)
% K.ma  = 0.32; % Seismic mass inside the accelerometer (kg)
% 
% % Accelerometer dynamics parameters
% K.ba_min = 1.5e3; % Minimum damper coefficient for accelerometer (Ns/m)
% K.ba_max = 2e4; % Maximum damper coefficient for accelerometer (Ns/m)
% K.ka_min = 5e-5; % Minimum spring stiffness for accelerometer (N/m)
% K.ka_max = 3e-3; % Maximum spring stiffness for accelerometer (N/m)
% K.Kacc = 1; % Proportionality constant relating seismic velocity to output voltage (Vs/m)
% K.xa_min = 1.95; 
% K.xa_max = 2.19;
% K.va_min = -2.57; 
% K.va_max = 2.57;
% 
% % Amplifier resistance ranges
% K.Rin_min = 0.1; % Minimum inverting resistance for amplifier (Ohms)
% K.Rin_max = 10; % Maximum inverting resistance for amplifier (Ohms)
% K.Rf_min = 1e4; % Minimum feedback resistance for amplifier (Ohms)
% K.Rf_max = 8e4; % Maximum feedback resistance for amplifier (Ohms)
% 
% % Solenoidal valve properties
% K.mv = 0.1; % Mass of the spool in the solenoidal valve (kg)
% K.kv = 1e3; % Spring stiffness for the solenoidal valve (N/m)
% K.bv = 1e3; % Damper coefficient for the solenoidal valve (Ns/m)
% K.alpha = 2.1e-2; % Solenoid constant (1/H)
% K.beta = -60; % Solenoid gain dependent on spool position (1/H·m)
% K.A0 = 4.7e-12; % Minimum cross-sectional area of the valve (m^2)
% K.xv_min = 0;
% K.xv_max_ic = 1; 
% K.xv_max= 1e-5; % Maximum stroke of the valve spool (m)
% K.vv_min = -2.19;
% K.vv_max = 2.19;
% K.I_min = -3.65;
% K.I_max = -3.25;
% 
% % Thruster and Xenon tank properties
% K.k = 1.66; % Heat capacity ratio (dimensionless)
% K.P_tank = 2e5; % Tank pressure (Pa)
% K.T_tank = 240; % Tank temperature (K)
% K.R = 63.32754; % Specific gas constant for Xenon (J/kg·K)
% K.q = 1.6e-19; % Charge of Xenon ions (C)
% K.DV = 2000; % Voltage across the acceleration grid (V)
% K.mi = 2.188e-25; % Mass of a single Xenon ion (kg)
% K.vext = sqrt(2*K.q*K.DV/K.mi); % Ions exit velocity (m/s)
% 
% % Orbital and drag properties
% K.oms = 1.658226e-6; % Secular pulsation due to orbital motion (rad/s)
% K.om0 = 1.160758e-3; % Orbital angular velocity (rad/s)
% 
% end


function K = constant()
% This function defines and returns a struct `K` containing all the 
% constants and parameter ranges used in the ACS simulation.
% Each parameter is named based on its physical significance.

% Satellite and accelerometer properties
K.Msc = 300; % Spacecraft mass (kg)
K.ma  = 0.32; % Seismic mass inside the accelerometer (kg)

% Accelerometer dynamics parameters
K.ba_min = 1.5e3; % Minimum damper coefficient for accelerometer (Ns/m)
K.ba_max = 2e4; % Maximum damper coefficient for accelerometer (Ns/m)
K.ka_min = 5e-5; % Minimum spring stiffness for accelerometer (N/m)
K.ka_max = 3e-3; % Maximum spring stiffness for accelerometer (N/m)
K.Kacc = 1; % Proportionality constant relating seismic velocity to output voltage (Vs/m)

% Amplifier resistance ranges
K.Rin_min = 0.1; % Minimum inverting resistance for amplifier (Ohms)
K.Rin_max = 10; % Maximum inverting resistance for amplifier (Ohms)
K.Rf_min = 1e4; % Minimum feedback resistance for amplifier (Ohms)
K.Rf_max = 8e4; % Maximum feedback resistance for amplifier (Ohms)

% Solenoidal valve properties
K.mv = 0.1; % Mass of the spool in the solenoidal valve (kg)
K.kv = 1e3; % Spring stiffness for the solenoidal valve (N/m)
K.bv = 1e3; % Damper coefficient for the solenoidal valve (Ns/m)
K.alpha = 2.1e-2; % Solenoid constant (1/H)
K.beta = -60; % Solenoid gain dependent on spool position (1/H·m)
K.A0 = 4.7e-12; % Minimum cross-sectional area of the valve (m^2)
K.xv_max = 1e-5; % Maximum stroke of the valve spool (m)

% Thruster and Xenon tank properties
K.k = 1.66; % Heat capacity ratio (dimensionless)
K.P_tank = 2e5; % Tank pressure (Pa)
K.T_tank = 240; % Tank temperature (K)
K.R = 63.32754; % Specific gas constant for Xenon (J/kg·K)
K.q = 1.6e-19; % Charge of Xenon ions (C)
K.DV = 2000; % Voltage across the acceleration grid (V)
K.mi = 2.188e-25; % Mass of a single Xenon ion (kg)
K.ve = sqrt(2*K.q*K.DV/K.mi); % Ions exit velocity (m/s)
rho = K.P_tank / (K.R * K.T_tank); %Xe density (kg/m3)
K.M_flux = sqrt(K.k * rho * K.P_tank ...
    * (2 / (K.k + 1))^((K.k + 1) / (K.k - 1))); %mass flux kg/(m2*s)

% Orbital and drag properties
K.oms = 1.658226e-6; % Secular pulsation due to orbital motion (rad/s)
K.om0 = 1.160758e-3; % Orbital angular velocity (rad/s)

end

function xdot = rhs(t, x, var)
% This function calculates the right-hand side of the system's ODEs.
% It models the dynamics of the ACS components: accelerometer, solenoidal valve, and ion thruster.

% Initialize the derivative vector with NaN values
xdot = nan(size(x));

% State variables
xa = x(1); % Accelerometer seismic mass position (m)
va = x(2); % Accelerometer seismic mass velocity (m/s)
xv = x(3); % Solenoidal valve spool position (m)
vv = x(4); % Solenoidal valve spool velocity (m/s)
I  = x(5); % Solenoid current (A)

% Input parameters from the `var` array
ba  = var(1); % Accelerometer damping coefficient (Ns/m)
ka  = var(2); % Accelerometer spring stiffness (N/m)
Rin = var(3); % Amplifier input resistance (Ohm)
Rf  = var(4); % Amplifier feedback resistance (Ohm)

% Load constant values
K = constant();

% Valve area as a function of spool position
xv = clip(xv, 0, K.xv_max);
Av = K.A0 + K.xv_max * (K.xv_max - xv); % Effective valve cross-sectional area (m^2)

% Xenon properties
rho = K.P_tank / (K.R * K.T_tank); % Xenon density in the tank (kg/m^3)
mdot = Av * sqrt(K.k * rho * K.P_tank * (2 / (K.k + 1))^((K.k + 1) / (K.k - 1))); 
% Mass flow rate of Xenon through the valve (kg/s)

% Solenoid inductance and its derivative with respect to spool position
dLdx = -K.beta  / (K.alpha + K.beta * xv)^2; % Derivative of inductance w.r.t. position (H/m)
L = 1 / (K.alpha + K.beta * xv); % Inductance as a function of spool position (H)

% Thrust generated by the ion thruster
T = mdot * sqrt(2 * K.q * K.DV / K.mi); % Thrust (N)

% Drag acting on the satellite
D = (2.2 - cos(K.oms * t) + 1.2 * sin(K.om0 * t) .* cos(K.om0 * t))*1e-3; % Drag force (N)

% Accelerometer output voltage
Vout = K.Kacc * va; % Output voltage proportional to velocity (V)

% Amplified output voltage
Vout_ampl = -Rf / Rin * Vout; % Amplifier output voltage (V)

% Dynamics equations
xdot(1) = va; % Time derivative of seismic mass position (velocity)
xdot(2) = (T - D) / K.Msc - ba / K.ma * va - ka / K.ma * xa; % Accelerometer seismic mass dynamics
xdot(3) = vv; % Time derivative of spool position (velocity)
xdot(4) = 1 / (2 * K.mv) * I^2 * dLdx - K.bv / K.mv * vv - K.kv / K.mv * xv; % Spool dynamics
xdot(5) = 1 / L * Vout_ampl; % Solenoid current dynamics

end

function plotStyle
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

function [T, D] = get_forces(tt, xx)

K = constant();
T = nan(size(tt));
D = T;

for k = 1 : length(T)

    xv = xx(k, 3);
    t  = tt(k);

    Av = K.A0 + K.xv_max * (K.xv_max - xv);
    rho = K.P_tank / (K.R * K.T_tank);
    mdot = Av * sqrt(K.k * rho * K.P_tank * (2 / (K.k + 1))^((K.k + 1) / (K.k - 1)));

    T(k) = mdot * sqrt(2*K.q * K.DV / K.mi);
    D(k) = (2.2 - cos(K.oms * t) + 1.2 * sin(K.om0 * t) .* cos(K.om0 * t))*1e-3;

end

end

