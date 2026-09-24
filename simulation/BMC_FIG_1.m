%% BMC_FIG_1.m
%
% This script generates panels B, C, and D of Figure 1 in the associated
% manuscript. Panel B shows a histogram of example latent period data
% (mean = 5.5 days) with fitted Exponential (K=1) and Erlang (K=4)
% distributions overlaid, where the Erlang shape and scale parameters
% were previously estimated using the GenErlangFit R package/Shiny app.
% Panels C and D show simulated log daily incidence and cumulative
% incidence trajectories, respectively, over 52 weeks, comparing a
% standard SEIR model (single-stage exposed and infectious periods)
% against an SE4IR model with four Erlang-distributed exposed substates.
% Both models are parameterized to produce identical early exponential
% growth rates, and both are simulated with and without a 50% reduction
% in transmission rate implemented at week 6, to illustrate differences
% in projected intervention impact between the two model structures.
%
% The ODE systems for both models are defined by the local function
% SEIR_ERLANG_ODE_v2 found at the end of this script, which implements a 
% generalized SEIR model with an arbitrary number of Erlang-distributed 
% substates in the exposed (E) and infectious (I) compartments, and is 
% solved numerically using MATLAB's ode45 solver.
%
% Output: a multi-panel figure (exported as workflow_figure.png)
% corresponding to Figure 1, panels B-D.

clear all;
clc;


%% Example latent period data and GenErlangFit parameters
% Example data: Generated from COVID-19 latent period distribution (Xin et al.)
% Time from exposure to detectable viral shedding, n=177
data = [4, 3, 11, 4, 4, 7, 5, 5, 11, 6, 4, 11, 4, 4, 5, 5, 5, 10, 1, 5, ...
        10, 4, 8, 9, 6, 8, 2, 2, 12, 5, 7, 6, 3, 4, 8, 5, 4, 3, 3, 4, ...
        2, 4, 7, 9, 4, 4, 6, 5, 4, 5, 7, 3, 3, 4, 7, 3, 6, 4, 6, 8, ...
        7, 8, 6, 4, 2, 7, 7, 5, 7, 3, 4, 9, 5, 5, 2, 5, 3, 8, 6, 5, ...
        3, 4, 3, 11, 3, 8, 13, 7, 4, 11, 3, 9, 5, 7, 5, 8, 7, 5, 2, 1, ...
        6, 2, 6, 3, 6, 3, 9, 3, 2, 6, 2, 9, 5, 3, 6, 5, 6, 5, 5, 7, ...
        5, 9, 6, 2, 7, 8, 2, 6, 6, 6, 6, 4, 9, 9, 8, 8, 6, 3, 6, 11, ...
        3, 6, 6, 2, 7, 1, 11, 3, 8, 6, 5, 8, 8, 2, 6, 7, 5, 7, 4, 5, ...
        5, 8, 1, 10, 4, 6, 7, 7, 8, 4, 1, 6, 8, 2, 7, 3, 4];

Mu = mean(data); % Mean waiting time in E
OG_sigma = 1/Mu; % 1 / Mean waiting time in E


% Gen Erlang Fit 
% FitR = GenErlang_Fit('erlang',data);

% FitR = [shape, scale] parameters for the fitted Erlang distribution,
% as identified by the GenErlangFit R package/Shiny app (K=4 substates
% was selected as optimal for this dataset). Replace the values below
% with the actual output from the app for your data.
% FitR = [K, λEr];   % <-- UPDATE with values from GenErlangFit app output

FitR = [4, 0.721];   % <-- UPDATE with values from GenErlangFit app output

% Histogram with proper bin width
n = length(data);
binWidth = 2 * iqr(data) / (n^(1/3));

%% Run SEIR Simulations 

% Pop Params 
N = 1e6;
T = 400;
tspan = 0:T;

% SEIR Params
Ro_1 = 2.1;
gamma = 0.5;
% sigma = 0.1;
sigma = OG_sigma;
beta_1 = Ro_1*gamma;

% Solve for early exp growth rate
r_eq = [1/(gamma*sigma),1/gamma + 1/sigma,1-Ro_1];
x = roots(r_eq);
r_1 = max(x);

nE = 4; % Would have been obtained from GenErlangFit
nI = 1;

% S(E4)IR Params for equivalent exp growth rate
Ro_n = (r_1/gamma) * ((1 + r_1/(sigma*nE))^nE) / (1 - (1 + r_1/(gamma*nI))^(-nI));
beta_n = Ro_n*gamma;

%% Simulate K=1 SEIR
p1.N = N;
p1.beta = beta_1;
p1.sigma = sigma;
p1.gamma = gamma;
p1.Esubstates = 1;
p1.Isubstates = 1;
p1.eps = 0;

y0_1 = [0; 1; 0];

[t1, y1] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p1), tspan, y0_1);

E1 = y1(:, 1);
I1 = y1(:, 2);
R1 = y1(:, 3);
S1 = N - sum(y1, 2);

%% Simulate K=4 SEEEEIIR
p4.N = N;
p4.beta = beta_n;
p4.sigma = sigma*nE; % Need to correct
p4.gamma = gamma;
p4.Esubstates = nE;
p4.Isubstates = nI;
p4.eps = 0;

y0_n = [zeros(nE,1); 1; 0];
[t4, y4] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p4), tspan, y0_n);
E4 = sum(y4(:, 1:nE), 2);
I4 = y4(:, nE+1);
R4 = y4(:, end);
S4 = N - sum(y4, 2);

%% Verify growth rates match
[~, peak_idx_1] = max(I1);
[~, peak_idx_4] = max(I4);

early_idx_1 = 1:round(peak_idx_1/3);
early_idx_4 = 1:round(peak_idx_4/3);

p_fit_1 = polyfit(t1(early_idx_1), log(I1(early_idx_1)), 1);
p_fit_4 = polyfit(t4(early_idx_4), log(I4(early_idx_4)), 1);

fprintf('=== Growth Rate Verification ===\n');
fprintf('Target r: %.4f\n', r_1);
fprintf('K=1 fitted r: %.4f\n', p_fit_1(1));
fprintf('K=4 fitted r: %.4f\n\n', p_fit_4(1));


%% Intervention Analysis: Effect of Beta Reduction on K=1 vs K=4

intervention_day = 45;
beta_reduction = 0.5;  % 50% reduction

%% Simulate K=1 with intervention
% Before intervention
tspan_before = 0:intervention_day;
[t1_before, y1_before] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p1), tspan_before, y0_1);

% After intervention
p1_after = p1;
p1_after.beta = beta_1 * beta_reduction;
y0_1_after = y1_before(end, :)';
tspan_after = intervention_day:T;
[t1_after, y1_after] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p1_after), tspan_after, y0_1_after);

% Combine
t1_int = [t1_before; t1_after(2:end)];
y1_int = [y1_before; y1_after(2:end, :)];

% CORRECTED: Cumulative infections = (I + R) / N * 100
I1_int = y1_int(:, 2);
R1_int = y1_int(:, 3);
cum_inf_1_int = (I1_int + R1_int) / N * 100;

%% Simulate K=4 with intervention
% Before intervention
[t4_before, y4_before] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p4), tspan_before, y0_n);

% After intervention
p4_after = p4;
p4_after.beta = beta_n * beta_reduction;
y0_4_after = y4_before(end, :)';
[t4_after, y4_after] = ode45(@(t,y) SEIR_ERLANG_ODE_v2(t,y,p4_after), tspan_after, y0_4_after);

% Combine
t4_int = [t4_before; t4_after(2:end)];
y4_int = [y4_before; y4_after(2:end, :)];

% CORRECTED: Cumulative infections = (I + R) / N * 100
I4_int = y4_int(:, nE+1);
R4_int = y4_int(:, end);
cum_inf_4_int = (I4_int + R4_int) / N * 100;

%% CORRECTED: Calculate cumulative infections for no intervention cases
cum_inf_1 = (I1 + R1) / N * 100;
cum_inf_4 = (I4 + R4) / N * 100;


%% Updated Figures for publication
% PUB Ready Figure - MATLAB Ver (2x2 Layout)
figure('Position', [100 100 1800 1200]);
FS = 24;          % Tick labels
FS_label = 20;    % Axis labels
FS_legend = 22;   % Legend text
FS_title = 16;    % Titles

% Color scheme (colorblind-friendly)
col_k1 = [0.8500, 0.3250, 0.0980]; % Red-orange
col_k4 = [0, 0.4470, 0.7410];      % Blue

% Subplot (1,2): Histogram with fitted distributions
subplot(2,2,2);

% Histogram with proper bin width
n = length(data);
binWidth = 2 * iqr(data) / (n^(1/3));
histogram(data, 'Normalization', 'pdf', 'FaceColor', "#A3B899", ...
    'FaceAlpha', 0.6, 'EdgeAlpha', 0.1, 'DisplayName', 'Example Data');

hold on;

% X-axis for PDFs
x = linspace(0, max(data), 500);

% Exponential PDF (rate = 1/mean)
h1_hist = plot(x, exppdf(x, mean(data)), '-', 'Color', col_k1, ...
    'LineWidth', 5, 'DisplayName', 'Exp Dist (K=1)');

% Fitted Erlang (Gamma) PDF
h2_hist = plot(x, gampdf(x, FitR(1), 1/FitR(2)), '-', 'Color', col_k4, ...
    'LineWidth', 5, 'DisplayName', sprintf('Erlang Dist (K=%d)', round(FitR(1))));

hold off;

% Labels and formatting
xlabel('Latent Period (days)', 'FontSize', FS_label, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', FS_label, 'FontWeight', 'bold');
% title('Erlang Distribution Fit to Latent Period Data', ...
%     'FontSize', FS_title, 'FontWeight', 'bold');

% Legend for histogram only
leg_hist = legend('Location', 'northeast', 'FontSize', FS_legend-2, 'Box', 'on');
leg_hist.ItemTokenSize = [50, 15];

% Grid and axes
% grid on;
% grid minor;
set(gca, 'FontSize', FS_label, 'LineWidth', 1.2, 'Box', 'on');
xlim([0, max(data)]);

% Subplot (2,1): Incidence plot (log scale)
subplot(2,2,3);

% CORRECTED: Daily incidence = diff(I + R)
inc_1 = [0; diff(I1 + R1)];
inc_4 = [0; diff(I4 + R4)];

% For intervention scenarios
inc_1_int = [0; diff(I1_int + R1_int)];
inc_4_int = [0; diff(I4_int + R4_int)];

% Plot curves (save handles for shared legend)
hold on;
h1 = plot(t1/7, inc_1, '-', 'Color', col_k1, 'LineWidth', 5);
h2 = plot(t4/7, inc_4, '-', 'Color', col_k4, 'LineWidth', 5);
h3 = plot(t1_int/7, inc_1_int, '--', 'Color', col_k1, 'LineWidth', 5);
h4 = plot(t4_int/7, inc_4_int, '--', 'Color', col_k4, 'LineWidth', 5);

% Intervention line
xline(intervention_day/7, '--k', 'LineWidth', 2, 'Alpha', 0.7);

hold off;

% Labels and formatting
xlabel('Time (weeks)', 'FontSize', FS_label, 'FontWeight', 'bold');
ylabel('Daily Incidence (per million)', 'FontSize', FS_label, 'FontWeight', 'bold');
% title('Simulation Incidence Trajectories', ...
%     'FontSize', FS_title, 'FontWeight', 'bold');

% Grid and axes
% grid on;
% grid minor;
set(gca, 'FontSize', FS_label, 'LineWidth', 1.2, 'Box', 'on');
set(gca, 'YScale', 'log');

xlim([0 T/7]);
ylim([1e-6, 1e-1]);
xticks([0 10 20 30 40 50]);
% Y-axis: counts from 1 to 100,000
ylim([1, 1e5]);
yticks([1, 10, 100, 1000, 10000, 100000]);
yticklabels({'1', '10', '100', '1,000', '10,000', '100,000'});
% Subplot (2,2): Cumulative % Infected plot
subplot(2,2,4);

% Plot curves
hold on;
plot(t1/7, cum_inf_1, '-', 'Color', col_k1, 'LineWidth', 5);
plot(t4/7, cum_inf_4, '-', 'Color', col_k4, 'LineWidth', 5);
plot(t1_int/7, cum_inf_1_int, '--', 'Color', col_k1, 'LineWidth', 5);
plot(t4_int/7, cum_inf_4_int, '--', 'Color', col_k4, 'LineWidth', 5);

% Intervention line
hv = xline(intervention_day/7, '--k', 'LineWidth', 2, 'Alpha', 0.7);

% Intervention annotation
text(intervention_day/7 + 0.7, 98, sprintf('Intervention (Week %d)', round(intervention_day/7)), ...
    'FontSize', FS_legend-5, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

hold off;

% Labels and formatting
xlabel('Time (weeks)', 'FontSize', FS_label, 'FontWeight', 'bold');
ylabel('Cumulative Infections (%)', 'FontSize', FS_label, 'FontWeight', 'bold');
% title('Simulation Cumulative Infection Trajectories', ...
%     'FontSize', FS_title, 'FontWeight', 'bold');

% Grid and axes - extend xlim to accommodate annotations
% grid on;
% grid minor;
set(gca, 'FontSize', FS_label, 'LineWidth', 1.2, 'Box', 'on', 'Layer', 'top');
xlim([0 T/7*1.1]);
xticks([0 10 20 30 40 50]);
ylim([0 100]);

% Final attack rate annotations (inside plot area)
text(T/7 + 0.9, cum_inf_1(end)-1, sprintf('%.1f%%', cum_inf_1(end)), ...
    'Color', col_k1, 'FontSize', FS_legend-3, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(T/7 + 0.9, cum_inf_4(end)+1, sprintf('%.1f%%', cum_inf_4(end)), ...
    'Color', col_k4, 'FontSize', FS_legend-3, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(T/7 + 0.9, cum_inf_1_int(end), sprintf('%.1f%%', cum_inf_1_int(end)), ...
    'Color', col_k1, 'FontSize', FS_legend-3, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(T/7 + 0.9, cum_inf_4_int(end), sprintf('%.1f%%', cum_inf_4_int(end)), ...
    'Color', col_k4, 'FontSize', FS_legend-3, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

xlim([0 T/7*1.15])

% Reposition subplots to make room for shared legend at bottom
sp1 = subplot(2,2,2);
sp2 = subplot(2,2,3);
sp3 = subplot(2,2,4);

% Get current positions
pos1 = sp1.Position;
pos2 = sp2.Position;
pos3 = sp3.Position;

% Shift all up and compress slightly to make room at bottom
sp1.Position = [pos1(1), pos1(2)+0.08, pos1(3), pos1(4)-0.04];
sp2.Position = [pos2(1), pos2(2)+0.08, pos2(3), pos2(4)-0.04];
sp3.Position = [pos3(1), pos3(2)+0.08, pos3(3), pos3(4)-0.04];

% Create shared legend at bottom for subplots 3 and 4
leg_shared = legend([h1, h2, h3, h4], ...
    {'SEIR_ (no interv.)', ...
     sprintf('SE_%dIR (no interv.)', nE), ...
     sprintf('SEIR_ (50%% \\beta Reduction, week %d)', round(intervention_day/7)), ...
     sprintf('SE_%dIR (50%% \\beta Reduction, week %d)', nE, round(intervention_day/7))}, ...
    'Orientation', 'horizontal', ...
    'FontSize', FS_legend, 'Box', 'on');

% Position legend in cleared space at bottom (positive y-value)
leg_shared.Position = [0.15, 0.01, 0.7, 0.05];
leg_shared.ItemTokenSize = [70, 15];

% Export for publication
exportgraphics(gcf, 'workflow_figure.png', 'Resolution', 300);
% exportgraphics(gcf, 'workflow_figure.pdf', 'ContentType', 'vector');

%%
function dydt = SEIR_ERLANG_ODE_v2(t,y,p)
    
    

    kE = p.Esubstates; % Number of Erlang E Substates 
    kI = p.Isubstates; % Number of Erlang I Substates

  

    dydt = zeros(kE+kI+1,1);
    
    % Equations for E compartment:
    dydt(1) = (1-p.eps) * (p.beta / p.N) * (sum(y(kE+1:kE+kI))) * ((p.N - sum(y))) - p.sigma * y(1);

    % y(2) to y(kE) exists as long as kE > 1
    Ecount = 1;
    while Ecount < kE
        dydt(Ecount + 1) = p.sigma * y(Ecount) - p.sigma * y(Ecount + 1);
        Ecount = Ecount + 1;
    end
    
    dydt(kE+1) = p.sigma*y(kE) - p.gamma*y(kE+1);

    % y(kE+2) to y(kE+kI) exists as long as kI > 1 :: DOUBLE CHECK
    Icount = 1;
    while Icount < kI
        dydt(Icount + 1 + kE) = p.gamma * y(Icount + kE) - p.gamma * y(Icount + 1 + kE);
        Icount = Icount + 1;
    end
    
    dydt(kI + kE + 1) = p.gamma*y(kI + kE);

end