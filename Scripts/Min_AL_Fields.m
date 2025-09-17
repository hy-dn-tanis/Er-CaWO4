clc
close all
clear all

vari = {"xx", "yy", "zz", "zzph"};

%% Load and Process Multiple Er:CWO Datasets at Different Magnetic Fields

% Zero Field (PS Off)
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-5dBm_7mK_26_4th_Decay_Meas_test_30-Apr-2025_@_21-02-19\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-5dBm_7mK_26_4th_Decay_Meas_test_@21-02-19.mat", ...
    vari{:});
freq_0mT = yy/1e9;                    % Frequency in GHz
xx_0mT = xx;                          % Time array
zz_mag_0mT = db2mag(zz);             % Transmission magnitude (linear)
zzph_rad_0mT = unwrap(deg2rad(zzph)); % Phase in radians (unwrapped)
zz_0mT = zz;                         % Original transmission in dB
clear yy xx zz zzph

% 0.5 mT Field
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2654.06_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_10-51-16\mw_time_spectrum_r19s2_nat_Er_CWO_2654.06_MHz_-5dBm_7mK_FieldBatch_@10-51-16.mat", ...
    vari{:});
freq_0p5mT = yy/1e9;
xx_0p5mT = xx;
zz_mag_0p5mT = db2mag(zz);
zzph_rad_0p5mT = unwrap(deg2rad(zzph));
zz_0p5mT = zz;
clear yy xx zz zzph

% 0.75 mT Field
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2653.62_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_14-07-49\mw_time_spectrum_r19s2_nat_Er_CWO_2653.62_MHz_-5dBm_7mK_FieldBatch_@14-07-49.mat", ...
    vari{:});
freq_0p75mT = yy/1e9;
xx_0p75mT = xx;
zz_mag_0p75mT = db2mag(zz);
zzph_rad_0p75mT = unwrap(deg2rad(zzph));
zz_0p75mT = zz;
clear yy xx zz zzph

% 1.0 mT Field
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2652.6_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_17-23-49\mw_time_spectrum_r19s2_nat_Er_CWO_2652.6_MHz_-5dBm_7mK_FieldBatch_@17-23-49.mat", ...
    vari{:});
freq_1mT = yy/1e9;
xx_1mT = xx;
zz_mag_1mT = db2mag(zz);
zzph_rad_1mT = unwrap(deg2rad(zzph));
zz_1mT = zz;
clear yy xx zz zzph

% 1.25 mT Field
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2651.3_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_20-39-43\mw_time_spectrum_r19s2_nat_Er_CWO_2651.3_MHz_-5dBm_7mK_FieldBatch_@20-39-43.mat", ...
    vari{:});
freq_1p25mT = yy/1e9;
xx_1p25mT = xx;
zz_mag_1p25mT = db2mag(zz);
zzph_rad_1p25mT = unwrap(deg2rad(zzph));
zz_1p25mT = zz;
clear yy xx zz zzph

% 1.5 mT Field
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2649.5_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_23-55-42\mw_time_spectrum_r19s2_nat_Er_CWO_2649.5_MHz_-5dBm_7mK_FieldBatch_@23-55-42.mat", ...
    vari{:});
freq_1p5mT = yy/1e9;
xx_1p5mT = xx;
zz_mag_1p5mT = db2mag(zz);
zzph_rad_1p5mT = unwrap(deg2rad(zzph));
zz_1p5mT = zz;
clear yy xx zz zzph

%% Summary of loaded datasets
fprintf('Loaded %d datasets with magnetic fields: 0, 0.5, 0.75, 1.0, 1.25, 1.5 mT\n', 6);
%% Create magnetic field array with offset correction
B = [0, 0.3, 0.55, 0.8, 1.05, 1.3];  % mT

fprintf('Magnetic field values (corrected): %.2f, %.2f, %.2f, %.2f, %.2f, %.2f mT\n', B);

%% Apply Circle Correction to All Datasets

% Zero Field (PS Off)
[zz_corrected_0mT, zzph_corrected_0mT] = circle_correction(zz_mag_0mT, zzph_rad_0mT, 3);

% 0.5 mT Field
[zz_corrected_0p5mT, zzph_corrected_0p5mT] = circle_correction(zz_mag_0p5mT, zzph_rad_0p5mT, 3);

% 0.75 mT Field
[zz_corrected_0p75mT, zzph_corrected_0p75mT] = circle_correction(zz_mag_0p75mT, zzph_rad_0p75mT, 3);

% 1.0 mT Field
[zz_corrected_1mT, zzph_corrected_1mT] = circle_correction(zz_mag_1mT, zzph_rad_1mT, 3);

% 1.25 mT Field
[zz_corrected_1p25mT, zzph_corrected_1p25mT] = circle_correction(zz_mag_1p25mT, zzph_rad_1p25mT, 3);

% 1.5 mT Field
[zz_corrected_1p5mT, zzph_corrected_1p5mT] = circle_correction(zz_mag_1p5mT, zzph_rad_1p5mT, 3);

clear zz_mag_0mT zz_mag_0p5mT zz_mag_0p75mT zz_mag_1mT zz_mag_1p25mT zz_mag_1p5mT
clear zzph_rad_0mT zzph_rad_0p5mT zzph_rad_0p75mT zzph_rad_1mT zzph_rad_1p25mT zzph_rad_1p5mT

fprintf('Circle correction applied to all %d datasets and uncorrected variables cleared\n', 6);

%% Apply min_aL routine to all corrected datasets (updated to get areas)

% Zero Field (PS Off)
[min_alpha_L_0mT, areas_0mT] = min_aL(zz_corrected_0mT, freq_0mT);

% 0.5 mT Field
[min_alpha_L_0p5mT, areas_0p5mT] = min_aL(zz_corrected_0p5mT, freq_0p5mT);

% 0.75 mT Field
[min_alpha_L_0p75mT, areas_0p75mT] = min_aL(zz_corrected_0p75mT, freq_0p75mT);

% 1.0 mT Field
[min_alpha_L_1mT, areas_1mT] = min_aL(zz_corrected_1mT, freq_1mT);

% 1.25 mT Field
[min_alpha_L_1p25mT, areas_1p25mT] = min_aL(zz_corrected_1p25mT, freq_1p25mT);

% 1.5 mT Field
[min_alpha_L_1p5mT, areas_1p5mT] = min_aL(zz_corrected_1p5mT, freq_1p5mT);

fprintf('min_aL routine applied to all %d corrected datasets - absorption and areas calculated\n', 6);


%%

[fit_obj_0mT, gof_0mT, params_0mT] = fit_triple_exponential(xx_0mT, areas_0mT);
[fit_obj_0p5mT, gof_0p5mT, params_0p5mT] = fit_triple_exponential(xx_0p5mT, areas_0p5mT);
[fit_obj_0p75mT, gof_0p75mT, params_0p75mT] = fit_triple_exponential(xx_0p75mT, areas_0p75mT);
[fit_obj_1mT, gof_1mT, params_1mT] = fit_triple_exponential(xx_1mT, areas_1mT);
[fit_obj_1p25mT, gof_1p25mT, params_1p25mT] = fit_triple_exponential(xx_1p25mT, areas_1p25mT);
[fit_obj_1p5mT, gof_1p5mT, params_1p5mT] = fit_triple_exponential(xx_1p5mT, areas_1p5mT);

% %% Create multi-panel figure showing field dependence
% figure('Position', [100, 100, 800, 1000]);
% 
% % Define subplot order (bottom to top: lowest to highest field)
% field_order = [1, 2, 3, 4, 5, 6]; % indices for B array
% datasets = {min_alpha_L_0mT, min_alpha_L_0p5mT, min_alpha_L_0p75mT, ...
%            min_alpha_L_1mT, min_alpha_L_1p25mT, min_alpha_L_1p5mT};
% freq_arrays = {freq_0mT, freq_0p5mT, freq_0p75mT, freq_1mT, freq_1p25mT, freq_1p5mT};
% 
% for i = 1:6
%     subplot(6, 1, 7-i); % 7-i to put lowest field at bottom
%     idx = field_order(i);
% 
%     % Plot the data
%     plot_traces(freq_arrays{idx}, datasets{idx}, 1);
% 
%     % Set x-axis limits
%     xlim([2.646, 2.658]);
% 
%     % Add B-field value in bottom left corner
%     text(0.05, 0.15, sprintf('B = %.2f mT', B(idx)), ...
%          'Units', 'normalized', 'HorizontalAlignment', 'left', ...
%          'VerticalAlignment', 'bottom', 'FontSize', 10, ...
%          'BackgroundColor', 'white', 'EdgeColor', 'black');
% 
%     % Y-axis label only on middle subplot - larger and bold
%     if i == 3
%         ylabel('-$\alpha L$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
%     end
% 
%     % X-axis label only on bottom subplot
%     if i == 1
%         xlabel('Frequency (GHz)', 'FontSize', 12);
%     end
% end
% 
% % Link x-axes for synchronized zooming/panning
% linkaxes(findall(gcf, 'type', 'axes'), 'x');

%% Plot all fit parameters as function of magnetic field
figure('Position', [100, 100, 1200, 600]);

% Collect all parameter structures
all_params = {params_0mT, params_0p5mT, params_0p75mT, params_1mT, params_1p25mT, params_1p5mT};

% Extract parameter values and confidence intervals
param_names = {'a', 'tau1', 'b', 'tau2', 'c', 'tau3', 'offset'};
param_labels = {'Amplitude a', 'Time constant \tau_1 (μs)', 'Amplitude b', ...
               'Time constant \tau_2 (μs)', 'Amplitude c', 'Time constant \tau_3 (μs)', 'Offset'};

% Initialize arrays
n_fields = length(B);
param_values = zeros(n_fields, length(param_names));
param_ci_lower = zeros(n_fields, length(param_names));
param_ci_upper = zeros(n_fields, length(param_names));

% Extract values
for i = 1:n_fields
    for j = 1:length(param_names)
        param_values(i, j) = all_params{i}.(param_names{j});
        param_ci_lower(i, j) = all_params{i}.([param_names{j} '_ci_lower']);
        param_ci_upper(i, j) = all_params{i}.([param_names{j} '_ci_upper']);
    end
end

% Reorganize parameters: amplitudes + offset in top row, time constants in bottom row
amplitude_params = {'a', 'b', 'c', 'offset'};
amplitude_labels = {'Amplitude a', 'Amplitude b', 'Amplitude c', 'd'};
time_params = {'tau1', 'tau2', 'tau3'};
time_labels = {'Time constant \tau_1 (\mus)', 'Time constant \tau_2 (\mus)', 'Time constant \tau_3 (\mus)'};

% Top row: Amplitudes and d (offset)
for i = 1:4
    subplot(2, 4, i);
    
    param_idx = find(strcmp(param_names, amplitude_params{i}));
    
    % Plot without error bars
    plot(B, param_values(:, param_idx), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    % Formatting
    xlabel('Magnetic Field (mT)');
    ylabel(amplitude_labels{i});
    grid on;
    title(sprintf('%s vs B-field', amplitude_labels{i}));
    
    % Set reasonable y-limits for positive parameters (except offset which can be negative)
    if ~strcmp(amplitude_params{i}, 'offset')
        ylim_current = ylim;
        if ylim_current(1) < 0
            ylim([0, ylim_current(2)]);
        end
    end
end

% Bottom row: Time constants
for i = 1:3
    subplot(2, 4, 4 + i);
    
    param_idx = find(strcmp(param_names, time_params{i}));
    
    % Plot without error bars
    plot(B, param_values(:, param_idx), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    % Formatting
    xlabel('Magnetic Field (mT)');
    ylabel(time_labels{i});
    grid on;
    title(sprintf('%s vs B-field', time_labels{i}));
    
    % Set reasonable y-limits (avoid negative for time constants)
    ylim_current = ylim;
    if ylim_current(1) < 0
        ylim([0, ylim_current(2)]);
    end
end

sgtitle('Triple Exponential Fit Parameters vs Magnetic Field', 'FontSize', 16);
%% Create individual figures for each dataset: Area vs Time with fits (clean version)

% Dataset arrays
xx_datasets = {xx_0mT, xx_0p5mT, xx_0p75mT, xx_1mT, xx_1p25mT, xx_1p5mT};
areas_datasets = {areas_0mT, areas_0p5mT, areas_0p75mT, areas_1mT, areas_1p25mT, areas_1p5mT};
fit_datasets = {fit_obj_0mT, fit_obj_0p5mT, fit_obj_0p75mT, fit_obj_1mT, fit_obj_1p25mT, fit_obj_1p5mT};

for i = 1:6
    figure('Position', [100 + i*50, 100 + i*30, 800, 600]);
    
    % Get current dataset
    xx = xx_datasets{i};
    areas = areas_datasets{i};
    f_area = fit_datasets{i};
    
    % Plot area data
    plot(xx, areas, 'm.')
    set(gca, 'XScale', 'log') % Base-10 logarithmic x-axis
    grid on; grid minor
    
    % Larger font sizes for axis labels
    xlabel('$\log_{10}(t)$', 'Interpreter', 'latex', 'FontSize', 18)
    ylabel('$\int (-\alpha L)\, df$', 'Interpreter', 'latex', 'FontSize', 18)
    
    hold on
    % Plot fitted function - create denser time array for smooth curve on log scale
    xx_dense = logspace(log10(min(xx)), log10(max(xx)), 1000);
    plot(xx_dense, f_area(xx_dense), 'Color', 'Black', 'LineWidth', 2)
    
    % Adjust tick label font size
    set(gca, 'FontSize', 16)
    
    % Larger legend font
    legend({"Data", "Tri-Exponential Fit"}, 'FontSize', 16, 'FontName', 'Times New Roman')
    
    % Add title with B-field value
    title(sprintf('B = %.2f mT', B(i)), 'FontSize', 18)
    
    hold off
end

fprintf('Created %d clean decay plots with fitted curves\n', 6);
    %% Helping functions

function [ax] = plot_traces(freq, input, steps)
    cc = turbo(size(input,2));
    for i = 1:steps:size(input,2)
        plot(freq, input(:,i), "Color",cc(i,:));
        hold on
    end
end

function [ax] = plot_traces2(freq, input, steps, xx, tickinterval)
    cc = turbo(size(input,2));
    for i = 1:steps:size(input,2)
        plot(freq, input(:,i), "Color",cc(i,:));
        hold on
        xlabel("Frequency (GHz)",'FontName','Helvetica');
    end
    colormap(cc);
    cb = colorbar;
    cb.Box = 'on';
    cb.Label.String = 'Time (s)';
    cb.Label.FontSize =13 ;
    cb.Label.FontName = 'Helvetica';

    xx_norm = (xx - min(xx)) / (max(xx) - min(xx));

    % Choose desired ticks (e.g., every 500 seconds)
    tick_interval = tickinterval;
    desired_ticks = min(xx)-23:tick_interval:max(xx)-23;
    [~, tick_indices] = min(abs(xx - desired_ticks'), [], 2);

    % Round tick labels to nearest integer
    cb.Ticks = xx_norm(tick_indices);
    cb.TickLabels = string(round(desired_ticks));
end

function [min_alpha_L,areas] = min_aL(zz_corrected,freq)
    min_alpha_L = zeros(size(zz_corrected));
    intensity = zeros(size(zz_corrected));
    idx_l = 1:40;
    idx_r = size(zz_corrected,1)-39:size(zz_corrected,1);

    for i = 1:size(min_alpha_L, 2)
        
        %construct baseline
        mean_l = mean(zz_corrected(idx_l,i));
        mean_r = mean(zz_corrected(idx_r,i));
        
        m = (mean_r - mean_l)./(freq(end)-freq(1));
        c = mean_l - m.*freq(1);
        y_val = m .*freq + c;
        
        %calculate intensity and absorption
        intens_temp = zz_corrected(:,i)./y_val';
        min_alpha_L_temp = log(intens_temp);
    
        %save values in big matrix & iterate
        intensity(:,i) = intens_temp;
        min_alpha_L(:,i) = min_alpha_L_temp;
    
    end

    areas_temp = trapz(min_alpha_L);
    areas = smoothdata(areas_temp, 'movmean', 3);
end

function [fit_obj, gof, fit_params] = fit_triple_exponential(xx, areas)
% FIT_TRIPLE_EXPONENTIAL Fits a triple exponential decay to time-resolved data
%
% Inputs:
%   xx    - time array
%   areas - integrated absorption areas vs time
%
% Outputs:
%   fit_obj    - fit object from MATLAB's fit function
%   gof        - goodness of fit structure
%   fit_params - structure containing fit parameters and confidence intervals
%
% Fit function: a*exp(-x/tau1) + b*exp(-x/tau2) + c*exp(-x/tau3) + offset

    % Define the fit function
    fit_func = fittype('a * exp(-1/tau1*x) + b * exp(-1/tau2*x) + c * exp(-1/tau3*x) + offset', ...
        'dependent', 'y', 'independent', 'x', ...
        'coefficients', {'a', 'tau1', 'b', 'tau2', 'c', 'tau3', 'offset'});
    
    % Fit parameters
    confinterval = 0.95;
    lower_bounds = [0, 0, 0, 0, 0, 0, -Inf];
    upper_bounds = [5, Inf, 5, Inf, 5, Inf, Inf];
    start_point = [0.1, 600, 0.1, 90, 0.2, 50, 0];
    
    % Perform the fit
    [fit_obj, gof] = fit(xx(:), areas(:), fit_func, ...
        'StartPoint', start_point, ...
        'Lower', lower_bounds, ...
        'Upper', upper_bounds, ...
        'Robust', 'LAR');
    
    % Extract fit parameters and confidence intervals
    coeff_values = coeffvalues(fit_obj);
    coeff_names = coeffnames(fit_obj);
    conf_int = confint(fit_obj, confinterval);
    
    % Create output structure with parameters
    fit_params = struct();
    for i = 1:length(coeff_names)
        fit_params.(coeff_names{i}) = coeff_values(i);
        fit_params.([coeff_names{i} '_ci_lower']) = conf_int(1, i);
        fit_params.([coeff_names{i} '_ci_upper']) = conf_int(2, i);
    end
    
    % Add goodness of fit metrics
    fit_params.rsquare = gof.rsquare;
    fit_params.adjrsquare = gof.adjrsquare;
    fit_params.rmse = gof.rmse;
    
    % Display results
    fprintf('Triple Exponential Fit Results:\n');
    fprintf('R-square: %.4f\n', gof.rsquare);
    fprintf('Adjusted R-square: %.4f\n', gof.adjrsquare);
    fprintf('RMSE: %.6f\n', gof.rmse);
    fprintf('\nFit Parameters (95%% confidence intervals):\n');
    for i = 1:length(coeff_names)
        fprintf('%s = %.4f [%.4f, %.4f]\n', coeff_names{i}, ...
                coeff_values(i), conf_int(1,i), conf_int(2,i));
    end
    
end

