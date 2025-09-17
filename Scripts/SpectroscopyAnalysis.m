clc
clear all
close all

vari = {"xx", "yy", "zz", "zzph", "zzed", "n_zz", "n_zzph", "n_zzed", "vna_pow"};

%Broadband Spectroscopy
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_100_4000MHz_7mK_Broadband_06-Mar-2025_@_15-05-52\power_sweep_r19s1_2_nat_Er_CWO_100_4000MHz_7mK_Broadband@_15-05-52.mat", ...
%     vari{:}) ;

%3022_3034MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_3022_3034MHz_7mK_Broadband_06-Mar-2025_@_15-35-29\power_sweep_r19s1_2_nat_Er_CWO_3022_3034MHz_7mK_Broadband@_15-35-29.mat", ...
%     vari{:});

% % 2700_2710MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_2700_2710MHz_7mK_Broadband_06-Mar-2025_@_15-49-36\power_sweep_r19s1_2_nat_Er_CWO_2700_2710MHz_7mK_Broadband@_15-49-36.mat", ...
%      vari{:});

% %2649_2659MHz
load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_2649_2659MHz_7mK_Broadband_06-Mar-2025_@_15-54-36\power_sweep_r19s1_2_nat_Er_CWO_2649_2659MHz_7mK_Broadband@_15-54-36.mat", ...
     vari{:});

% 2649_2659MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_2649_2659MHz_7mK_Broadband_06-Mar-2025_@_15-54-36\power_sweep_r19s1_2_nat_Er_CWO_2649_2659MHz_7mK_Broadband@_15-54-36.mat", ...
%      vari{:});

% 368_378MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_368_378MHz_7mK_Broadband_07-Mar-2025_@_10-30-04\power_sweep_r19s1_2_nat_Er_CWO_368_378MHz_7mK_Broadband@_10-30-04.mat", vari{:});

%3020_3034MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_3020_3034MHz_7mK_Broadband_07-Mar-2025_@_10-35-54\power_sweep_r19s1_2_nat_Er_CWO_3020_3034MHz_7mK_Broadband@_10-35-54.mat", vari{:});

%3160_3240MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_3160_3240MHz_7mK_Broadband_07-Mar-2025_@_10-56-57\power_sweep_r19s1_2_nat_Er_CWO_3160_3240MHz_7mK_Broadband@_10-56-57.mat", vari{:});
% 
% 2695_2715MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_2695_2715MHz_7mK_Broadband_07-Mar-2025_@_11-01-03\power_sweep_r19s1_2_nat_Er_CWO_2695_2715MHz_7mK_Broadband@_11-01-03.mat", vari{:});

%3475_3495MHz
% load("\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\power_sweep_r19s1_2_nat_Er_CWO_3475_3495MHz_7mK_Broadband_07-Mar-2025_@_10-42-19\power_sweep_r19s1_2_nat_Er_CWO_3475_3495MHz_7mK_Broadband@_10-42-19.mat", vari{:});
 

%%
%Unit Conversions

freq = yy/1e9;
% 
% id_left = find(freq == 3.18);
% id_right = find(freq == 3.19);
% 
id_left = 1;
id_right = length(freq);

zz_mag = db2mag(zz);
zzph_rad = unwrap(deg2rad(zzph));

%% Plotting

figure; 
    subplot(3,1,1)
        plot(freq(1,id_left:id_right),zz(id_left:id_right,:),"Color","k", "LineWidth",1.05)
        xlabel("Frequency (GHz)",'FontSize', ...
        12)
        ylabel("|S_{21}| (dB)",'FontSize', ...
        12)
        % xlim([3.181 3.19])
        % ylim([-25 -15])
        title("(a)",'FontSize', ...
        12, 'FontWeight', 'bold')
    subplot(3,1,2)
        plot(freq(1,id_left:id_right), zzph(id_left:id_right,:), "Color","k", "LineWidth",0.8)
        xlabel("Frequency (GHz)",'FontSize', ...
        12)
        % xlim([3.181 3.19])
        ylabel("Phase (deg)",'FontSize', ...
        12)
        
        title("(b)",'FontSize', ...
        12, 'FontWeight', 'bold')
    subplot(3,1,3)
        plot(freq(1,id_left:id_right),zzed(id_left:id_right,:), "Color","k", "LineWidth",0.8)
        xlabel("Frequency (GHz)",'FontSize', ...
        12)
        ylabel("Group Delay (s)",'FontSize', ...
        12)
        % xlim([3.181 3.19])
        title("(c)",'FontSize', ...
        12, 'FontWeight', 'bold')
    % sgtitle(strcat("Broadband Zero Field Spectroscopy ", num2str(freq(1,id_left)), " - ", num2str(freq(1,id_right)), " GHz"), 'FontSize', ...
    %     12, 'FontWeight', 'bold', 'FontName', 'Helvetica');

%%


    zz_corrected27 = zeros(size(zz_mag));
    zzph_corrected27 = zeros(size(zzph_rad));
    in_point27 = 100;

    %Final Trace
    
    i = size(zz_corrected27,2);
        trace27 = zz_mag(:,i);
        trace_ph27 = zzph_rad(:,i);

    
        % trace = n_zz(:,i);
        % trace_ph = n_zzph(:,i);
    
        data_real27 = trace27 .* cos(trace_ph27);
        data_imag27 = trace27 .* sin(trace_ph27);
        dataset27 = [data_real27 data_imag27];
    

        %27
        Pars27 = CircleFitByPratt(dataset27(400:2000,:));
        dataset_tr27 = [data_real27-Pars27(1) data_imag27-Pars27(2)];
        theta = linspace(0,2*pi, 360);
        x = Pars27(:,1) + Pars27(:,3) * cos(theta);
        y = Pars27(:,2) + Pars27(:,3) * sin(theta);


        phi_rad27 = atan2(dataset_tr27(in_point27,2), dataset_tr27(in_point27,1));
        Rot27 = [cos(-phi_rad27) -sin(-phi_rad27); sin(-phi_rad27) cos(-phi_rad27)];

        dataset_rot27 = dataset_tr27;
        for k = 1:size(dataset_rot27,1)
            points = transpose(dataset_tr27(k, :));
            rot_point = Rot27 * points;
            dataset_rot27(k,:) = transpose(rot_point);
        end

        dist27 = 1 - dataset_rot27(in_point27,1);
        dataset_norm27 = [dataset_rot27(:,1) + dist27, dataset_rot27(:,2)];


        zz_corrected27(:,i) = sqrt(dataset_norm27(:,1).^2 + dataset_norm27(:,2).^2);
        zzph_corrected27(:,i) = atan2(dataset_norm27(:,2), dataset_norm27(:,1));

       
        figure;
            plot(data_real27, data_imag27,'.');
            xlabel("$\mathrm{Re}|S_{21}|$",'Interpreter', 'latex')
            ylabel("$\mathrm{Imag}|S_{21}|$",'Interpreter', 'latex')
                    grid on;grid minor

            hold on
            plot(x, y, 'r', 'LineWidth', 2);  % 'r' for red, 'LineWidth' controls thickness            
            axis equal

    
        figure;
            plot(dataset_norm27(:,1), dataset_norm27(:,2),'.');
                    grid on;grid minor

            xlabel("$\mathrm{Re}|S_{21}|$",'Interpreter', 'latex')
            ylabel("$\mathrm{Imag}|S_{21}|$",'Interpreter', 'latex')
            axis equal

       


        % sgtitle("Spectroscopy Data in Complex Plane", 'FontName', 'Helvetica','FontWeight', 'bold')
        % 
        % figure;
        % 
        % subplot(2,2,1)
        % plot(freq, zz27,'.');
        % xlabel("freq (GHz)")
        % ylabel("|S_{21}| (Norm)")
        % title("Transmission (Original)")
        % subplot(2,2,2)
        % plot(freq, zzph27,'.');
        % xlabel('Re(zz)')
        % ylabel("Imag(zz)")
        % title("Phase (Original)")
        % 
        % subplot(2,2,3)
        % plot(freq, zz_corrected27,'.');
        % xlabel("freq (GHz)")
        % ylabel("|S_{21}| (Norm)")
        % title("Transmission (Corrected)")
        % subplot(2,2,4)
        % plot(freq, zzph_corrected27,'.');
        % xlabel("freq (GHz)")
        % ylabel("|S_{21}| (Norm)")
        % title("Phase (Corrected)")


%         figure;
% subplot(2,2,1)
% plot(freq, zz26);
% xlabel("freq (GHz)")
% ylabel("|S_{21}| (Norm)")
% title("Transmission (Original)")
% subplot(2,2,2)
% plot(freq, zzph26);
% xlabel('Re(zz)')
% ylabel("Imag(zz)")
% title("Phase (Original)")
% subplot(2,2,3)
% plot(freq, zz_corrected26);
% xlabel("freq (GHz)")
% ylabel("|S_{21}| (Norm)")
% title("Transmission (Corrected)")
% subplot(2,2,4)
% plot(freq, zzph_corrected26);
% xlabel("freq (GHz)")
% ylabel("|S_{21}| (Norm)")
% title("Phase (Corrected)")