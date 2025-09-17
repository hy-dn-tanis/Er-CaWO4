clc
close all

clearvars -except freqzf freq05 freq075 freq1 freq125 freq15 min_alpha_LZF min_alpha_L05 min_alpha_L075 min_alpha_L1 min_alpha_L125 min_alpha_L15

vari = {"bg_RF", "bg_PH", "xx", "yy", "zz", "zzph","zzed"...
    "pumping_zz", "pumping_zzph", "vna_f0", "vna_f1", "vna_fc", "vna_pow", "n_zz", "n_zzph","n_zzed"};




%Zero Field PS Off I
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-10dBm_7mK_26_3rd_Decay_Meas_test_16-Apr-2025_@_13-24-18\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-10dBm_7mK_26_3rd_Decay_Meas_test_@13-24-18.mat", ...
%     vari{:});

%Zero Field PS Off II
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-5dBm_7mK_26_4th_Decay_Meas_test_30-Apr-2025_@_21-02-19\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-5dBm_7mK_26_4th_Decay_Meas_test_@21-02-19.mat", ...
    vari{:});
% % 0.5mT
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2654.06_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_10-51-16\mw_time_spectrum_r19s2_nat_Er_CWO_2654.06_MHz_-5dBm_7mK_FieldBatch_@10-51-16.mat", ...
%     vari{:});
% % 0.75 mT
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2653.62_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_14-07-49\mw_time_spectrum_r19s2_nat_Er_CWO_2653.62_MHz_-5dBm_7mK_FieldBatch_@14-07-49.mat", ...
%     vari{:});
% % 1 mT
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2652.6_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_17-23-49\mw_time_spectrum_r19s2_nat_Er_CWO_2652.6_MHz_-5dBm_7mK_FieldBatch_@17-23-49.mat", ...
%     vari{:});
% % 1.25mT
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2651.3_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_20-39-43\mw_time_spectrum_r19s2_nat_Er_CWO_2651.3_MHz_-5dBm_7mK_FieldBatch_@20-39-43.mat", ...
%     vari{:});
% % 1.5mT
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\lowfieldbatch_time_spectrum\mw_time_spectrum_r19s2_nat_Er_CWO_2649.5_MHz_-5dBm_7mK_FieldBatch_01-May-2025_@_23-55-42\mw_time_spectrum_r19s2_nat_Er_CWO_2649.5_MHz_-5dBm_7mK_FieldBatch_@23-55-42.mat", ...
%     vari{:});
freq = yy/1e9;
clear yy

% Unit conversions                           

%Unit conversions
% Convert Transmission from dB to linear scale
zz_mag = db2mag(zz);
n_zz_mag = db2mag(n_zz);


% Convert phase to radian and unwrap
zzph_rad = unwrap(deg2rad(zzph));
n_zzph_rad = unwrap(deg2rad(n_zzph));
pumping_zzph_rad = unwrap(deg2rad(pumping_zzph));


%% Raw Data Plotting
% figure;
%     subplot(2,3,1); plot_traces(freq, zz_mag, 5);
%     title("(i) Transmission"); ylabel("S_{21} (dB)");
% 
%     subplot(2,3,2); plot_traces(freq, rad2deg(zzph_rad), 5);
%     title("(ii) Phase"); ylabel("Phase (deg)");
% 
%     subplot(2,3,3); plot_traces(freq, zzed, 5);
%     title("(iii) Group Delay"); ylabel("Group Delay (s)");
% 
%     subplot(2,3,4); plot_traces(freq, n_zz_mag, 5);
%     title("(iv) Normalized Transmission"); ylabel("S_{21} (dB)");
% 
%     subplot(2,3,5); plot_traces(freq, rad2deg(n_zzph_rad), 5);
%     title("(v) Normalized Phase"); ylabel("Phase (deg)");
% 
%     subplot(2,3,6); plot_traces(freq, n_zzed, 5);
%     title("(vi) Normalized Group Delay"); ylabel("Group Delay (s)");


%% Circle Correction


%
[zz_corrected, zzph_corrected] = circle_correction(zz_mag, zzph_rad, 3);
% 
% figure;
%     % subplot(2,2,1);
%     plot_traces(freq, zz_mag, 5); 
%     xlabel("freq (GHz)")
%     ylabel("|S_{21}|")
%     % title("Original Transmission")
%     axis tight
% figure;
%     % subplot(2,2,2);
%     plot_traces(freq, zzph_rad, 5);
%     xlabel("Frequency (GHz)", 'FontSize',13)
%     ylabel("Phase (deg)", 'FontSize',13)
%     % title("Original Phase")
%     axis tight
% figure;
%     % subplot(2,2,3);
%     plot_traces(freq, zz_corrected,2);
%     xlabel("Frequency (GHz)", 'FontSize',13)
%     ylabel("|S_{21}|", 'FontSize',13)
%     % title("Corrected Transmission")
%     axis tight
% 
% figure;
%     % subplot(2,2,4);
%     plot_traces(freq, zzph_corrected,2);
%     xlabel("Frequency (GHz)", 'FontSize',13)
%     ylabel("Phase (deg)", 'FontSize',13)
%     % title("Corrected Phase")
%     axis tight
% 
% % sgtitle("Decay Traces - 2.6 GHz transition", 'FontName', 'Helvetica','FontWeight', 'bold')

%%

min_alpha_L = zeros(size(zz));
intensity = zeros(size(zz));
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

    clear mean_l mean_r m c y_val intens_temp min_alpha_L_temp

end


areas = trapz(min_alpha_L);
areas_smooth = smoothdata(areas, 'movmean', 3);
%%
figure;
    

    plot_traces2(freq, min_alpha_L(:,2:end),3,xx, 1000);
    ylabel("- \alpha L",'FontName','Helvetica', 'FontSize',13);
    xlabel("Frequency (GHz)", 'FontSize',13,'FontName','Helvetica');
    axis tight

figure;

% Create surface plot
surf(xx, freq, min_alpha_L);

    shading interp;
    cb = colorbar;
    cb.Label.String = '$-\alpha L$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 13;
    cb.Label.FontWeight = 'bold';
    cb.Box = 'on';          % Adds a box around the colorbar
    cb.Ticks = -0.04:0.005:0;



    colormap(flipud(turbo));

    view(0, 90);
    axis tight;

xlabel('Time (s)', 'FontSize',13);
ylabel('Frequency','FontSize',13);
zlabel('$-\alpha L');

figure;

    plot(xx, areas_smooth, '.');
    set(gca, 'XScale', 'log')  % Apply decadic logarithmic scale
    xlabel("$\log_{10}(t)$", 'Interpreter','latex', 'FontWeight','bold')
    ylabel('$\int -\alpha L\, df$', 'Interpreter', 'latex', 'FontWeight','bold')
    title("Logarithmic Time Scale", 'Interpreter', 'latex', 'FontWeight', 'bold')


    % sgtitle("Integral of $- \alpha L$ over time", 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 15)
%%
%custom tri-exp fitting function


fit_func = fittype('a * exp(-1/tau1*x) + b * exp(-1/tau2*x) + c * exp(-1/tau3*x) + d', ...
    'dependent', 'y', 'independent', 'x', ...
    'coefficients', {'a', 'tau1', 'b', 'tau2', 'c', 'tau3', 'offset'});

confinterval = 0.95;
lower_bounds = [0, 0, 0, 0, 0, 0, 0];
upper_bounds = [5, Inf,5,Inf, 5,Inf, Inf];
%fit
[f_area, gof_area] = fit(xx', areas_smooth', fit_func, 'StartPoint',[0.1, 600, 0.1, 90, 0.2, 50, 1], ...
    'Lower', lower_bounds, 'Upper', upper_bounds, 'Robust', 'LAR');

%parameters
    fit_coeff_area = coeffvalues(f_area);
    
    ci_area = confint(f_area, confinterval);
    paramValues_area = coeffvalues(f_area);
    paramNames_area = coeffnames(f_area);

    amplitude_area = [fit_coeff_area(1,1), fit_coeff_area(1,3), fit_coeff_area(1,5)];
    ci_ampl_area = ci_area(:, [1,3,5]);
    
    T1_area = [fit_coeff_area(1,2), fit_coeff_area(1,4), fit_coeff_area(1,6)];
    ci_Time_area = ci_area(:, [2 4 6]); 

    T_area = table(fit_coeff_area(1,1), fit_coeff_area(1,3), fit_coeff_area(1,5), fit_coeff_area(1,2), fit_coeff_area(1,4), fit_coeff_area(1,6),fit_coeff_area(1,7), ...
    'VariableNames', {'a', 'b', 'c', 'tau1', 'tau2', 'tau3', 'offset'});

%%

% Plotting
figure;

% Plot area data using original xx values
    plot(xx, areas_smooth, 'm.')
    set(gca, 'XScale', 'log')  % Base-10 logarithmic x-axis
    grid on; grid minor
    % Larger font sizes for axis labels
    xlabel('$\log_{10}(t)$', 'Interpreter', 'latex', 'FontSize', 18)
    ylabel('$\int (-\alpha L)\, df$', 'Interpreter', 'latex', 'FontSize', 18)

hold on

% Plot fitted function
    plot(xx, f_area(xx), 'Color', 'Black', 'LineWidth', 2)
    % Adjust tick label font size
    set(gca, 'FontSize', 16)
    % Larger legend font
    legend({"Data", "Tri-Exponential Fit"}, 'FontSize', 16, 'FontName', 'Times New Roman')

 disp(T_area)  
 disp(gof_area)

% figure;
% 
%     subplot(2,1,1)
%     e = errorbar(1:length(T1_area), T1_area, abs(T1_area - ci_Time_area(1, :)), abs(ci_Time_area(2, :) - T1_area), 'o');
%     e.Marker = 'o';              
%     e.MarkerEdgeColor = 'b';   
%     e.MarkerSize = 5;  
% 
%     xticks(1:length(paramValues_area));
%     xticklabels(paramNames_area([2,4,6],:));
%     xlabel('Parameters');
%     ylabel('Value');
%     title(strcat("Time Constants with ", num2str(confinterval * 100), "% Confidence Interval"));
%     grid on; grid minor;
% 
%     subplot(2,1,2)
% 
%     e = errorbar(1:length(amplitude_area), amplitude_area, abs(amplitude_area - ci_ampl_area(1, :)), abs(ci_ampl_area(2, :) - amplitude_area), 'o');
%     e.Marker = 'o';              
%     e.MarkerEdgeColor = 'b';   
%     e.MarkerSize = 5;  
% 
%     xticks(1:length(paramValues_area));
%     xticklabels(paramNames_area([1,3,5],:));
%     xlabel('Parameters');
%     ylabel('Value');
%     title(strcat("Amplitude Parameters with ", num2str(confinterval * 100), "% Confidence Interval"));
%     grid on; grid minor;

    % legend("Amplitude Fitting", "Area Fitting")

    %% Helping functions

function [ax] = plot_traces(freq, input, steps)
    cc = turbo(size(input,2));
    for i = 1:steps:size(input,2)
        plot(freq, input(:,i), "Color",cc(i,:));
        hold on
        xlabel("Frequency (GHz)");
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


%%
