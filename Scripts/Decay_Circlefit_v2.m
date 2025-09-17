%% Circle Fit for Decay Measurement

clc;
clear all;
close all;


%Load Data

load('\\badwwmi-k04-qtm\data_qtm\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-10dBm_7mK_2.6GHz_2ndMeas_09-Apr-2025_@_10-35-08\mw_time_spectrum_r19s2_nat_Er_CWO_2654_MHz_-10dBm_7mK_2.6GHz_2ndMeas_@10-35-08.mat')
% % load('\\badwwmi-k04-qtm\Data_QTM\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2705_MHz_-10dBm_7mK_2.6GHz_2ndMeas_09-Apr-2025_@_13-46-27\mw_time_spectrum_r19s2_nat_Er_CWO_2705_MHz_-10dBm_7mK_2.6GHz_2ndMeas_@13-46-27.mat')

% load('\\badwwmi-k04-qtm\data_qtm\_data\run19\Er_CWO\mw_time_spectrum_r19s2_nat_Er_CWO_2653.5_MHz_-10dBm_7mK_LowFieldBatch_1_17-Apr-2025_@_10-56-13\mw_time_spectrum_r19s2_nat_Er_CWO_2653.5_MHz_-10dBm_7mK_LowFieldBatch_1_@10-56-13.mat')

freq = yy/1e9;
fr_idx = find(freq == vna_fc/1e9); %Take index for center frequency

%% Cut Range

% only if necessary
% fr_left = find(freq == 2.65);
% fr_right = find(freq == 2.658);

fr_left = 1;
fr_right = length(freq);

freq = freq(fr_left:fr_right);
zz = zz(fr_left:fr_right,:);
zzph = zzph(fr_left:fr_right, :);
% 
% zz_new = bg_RF(:,fr_left:fr_right);
% zzph_new = bg_PH(:,fr_left:fr_right);

%Normalized Data
% zz_new = n_zz(fr_left:fr_right,:);
% zzph_new = n_zzph(fr_left:fr_right, :);



%% Plotting
% 
cc = turbo(size(zz,2));
figure; 

for i = 1:size(zz,2)
    subplot(1,2,1)
    plot(freq, zz(:,i), 'Color',cc(i,:));
    hold on
    grid on
end
title('Transmission Spectrum')
xlabel('f (GHz)')
ylabel('S21 (dBm)')

for i = 1:size(zz,2)
    subplot(1,2,2)
    plot(freq, zzph(:,i), 'Color',cc(i,:));
    hold on
    grid on; grid minor
end
title('Phase spectrum')
xlabel('f (GHz)')
ylabel('Phase (deg)')
hold off

%%

cc = turbo(size(zz,2));
zz_corrected = [];
zzph_corrected = [];

for i = 1:size(zz,2)
% for i = 1
    zzmag = db2mag(zz(:,i));
    zzph_new = unwrap(zzph(:,i));
    zzph_new = deg2rad(zzph(:, i));

    % zzmag = db2mag(n_zz(:,i));
    % zzph_new = unwrap(n_zzph(:,i));
    % zzph_new = deg2rad(n_zzph(:, i));
    
    zzreal = zzmag .* cos(zzph_new);
    zzimag = zzmag .* sin(zzph_new);
    
    % figure;
    % plot(freq,zz(:,1))
    % title("trace 1")
    % 
    % figure;
    % plot(zzreal, zzimag, '.')
    % axis equal
    % title("Real and Imaginary Datapoints")
    
    dataset = [zzreal, zzimag];
    
    %Take parameters from spectroscopy correction
    % a = -0.0038;
    % b = 0.0305;
    % R = 0.0238;
    
    Pars = CircleFit(dataset);

    a = Pars(:,1);
    b = Pars(:,2);
    R = Pars(:,3);
    
    theta = linspace(0, 2*pi, 360);
    x = a + R * cos(theta);
    y = b + R * sin(theta);
    % 
    % figure;
    % plot(dataset(:,1), dataset(:,2), '.');
    % hold on
    % plot(x,y);
    % axis equal

    dataset_trans = [zzreal-a zzimag-b];
    
    dataset_rot = dataset_trans;
    phi_rad = atan2(dataset_trans(100,2), dataset_trans(100,1));
    % phi_rad = 2.13;
    Rot = [cos(-phi_rad) -sin(-phi_rad); sin(-phi_rad) cos(-phi_rad)];

    
    for i = 1:size(dataset_trans, 1)
        points = transpose(dataset_trans(i, :));
        rot_point = Rot * points;
        dataset_rot(i,:) = transpose(rot_point);
    end


    % figure
    % subplot(1,2,1)
    % plot(dataset_trans(:,1), dataset_trans(:,2), '.')
    % hold on
    % scatter(dataset_trans(200,1), dataset_trans(200,2), 'square', 'red')
    % grid on
    % grid minor
    % 
    % subplot(1,2,2)
    % plot(dataset_rot(:,1), dataset_rot(:,2), '.')
    % hold on
    % scatter(dataset_rot(200,1), dataset_rot(200,2), 'square', 'red')
    % grid on
    % grid minor
    
    min_re = min(dataset_rot(:,1));
    dataset_rot(:,1) = dataset_rot(:,1) - min_re;
    
    % figure(2)
    % plot(zzreal(:,1), zzimag(:,1), '.');
    % hold on
    % grid on
    % plot(x, y, '.');
    % scatter(a,b, 'square', 'filled', 'k')
    % yline(0, '--');
    % xline(0, '--');
    % plot(dataset_trans(:,1), dataset_trans(:,2), '.');
    % plot(dataset_rot(:,1), dataset_rot(:,2), '.');
    % axis equal
    % title("Datapoints")

    dist = 1 - dataset_rot(1,1);
    fin_data = [dataset_rot(:,1) + dist, dataset_rot(:,2)];
    % 
    % figure(3)
    % plot(fin_data(:,1), fin_data(:,2), '.')
    % grid on
    % grid minor
    % 
    zz_correct = sqrt(fin_data(:,1).^2 + fin_data(:,2).^2);
    zzph_correct = atan2(fin_data(:,2), fin_data(:,1));
    zzph_correct = rad2deg(zzph_correct);

    zz_corrected = [zz_corrected zz_correct];
    zzph_corrected = [zzph_corrected zzph_correct];

end

clear zz_correct zzph_correct;

%%

figure(2)

for i = 1:size(zz_corrected,2)
    subplot(2,2,1)
    plot(freq, zz(:,i), 'Color',cc(i,:));
    hold on
    grid on
    grid minor;
end
title('Original Transmission Spectrum')
xlabel('Frequency (GHz)')
ylabel('S21 (dBm)')

for i = 1:size(zz_corrected,2)
    subplot(2,2,2)
    plot(freq, zzph(:,i), 'Color',cc(i,:));
    hold on
    grid on
    grid minor;
end
title('Original Phase spectrum')
xlabel('Frequency (GHz)')
ylabel('Phase (degrees)')
hold off

for i = 1:size(zz_corrected, 2)
    subplot(2,2,3)
    plot(freq, zz_corrected(:, i), 'Color', cc(i,:));
    hold on
    grid on
    grid minor;
    title("Corrected Transmission")
    xlabel("Frequency (GHz)")
    ylabel("S21 (linear)")

    subplot(2,2,4)
    plot(freq, zzph_corrected(:,i), 'Color', cc(i,:));
    hold on
    grid on
    grid minor;
    title("Corrected Phase")
    xlabel("Frequency (GHz)")
    ylabel("Phase (degrees)")

end

%% Save Data

comment = '2nd_Decay_Meas_2.6GHz_Transition';
directory = "\\badwwmi-k04-qtm\Data_QTM\_data\run19\MeasurementProtocolHT\Circle Corrected Data"

filename_new = fullfile(directory, ['CorrectedData_' comment '.mat']);
save(filename_new)

disp("Data saved")

%% Helping functions

function Par = CircleFit(XY)

%--------------------------------------------------------------------------
%  
%     Circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------

n = size(XY,1);      % number of data points

centroid = mean(XY);   % the centroid of the data set

%     computing moments (note: all moments will be normed, i.e. divided by n)

Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;

for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end
   
Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;

%    computing the coefficients of the characteristic polynomial

Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Mxz2 = Mxz*Mxz;
Myz2 = Myz*Myz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;

epsilon=1e-12; 
ynew=1e+20;
IterMax=20;
xnew = 0;

%    Newton's method starting at x=0

for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
    if (abs(ynew)>abs(yold))
        disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
        xnew = 0;
        break;
    end
    Dy = A1 + xnew*(A22 + 16*xnew*xnew);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Pratt will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end

%    computing the circle parameters

DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];

end    %    CircleFitByPratt

