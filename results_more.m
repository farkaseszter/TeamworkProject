clc
close all
clear all
format long
% --- User settings ---
filename_d = 'proppel_15_down_1460_1240_07.csv';
headerLines = 7;
rpm_d = zeros(1,12);
rpm_d(1,1) = 1460;
for rpm_var = 1:11
    rpm_d(1,rpm_var+1) = rpm_d(1,rpm_var)-20;
end
rpm_plot_d = [rpm_d;rpm_d];

ends_d = [35 65;75 125;135 185;195 245;255 305;315 365;375 425;
    435 485;495 545;555 605;615 665;675 725];
time_diff_d = 50;

% --- Read the CSV file, skipping header lines ---
opts_d = detectImportOptions(filename_d);
opts_d.DataLines = [headerLines+1, Inf];   % start reading after header lines
data_d = readtable(filename_d, opts_d);
time_d = data_d{:, 2};
increment_d = time_d(2);

% --- Extract columns & transform to local coordinate system ---
data_all_d = zeros(3,3, length(time_d));
k_var = [3 2 1];
sens_var = [2 1 3];
for k = 1:3
    for k2 = 1:3
        for k3 = 1:length(time_d)
            if isnan(data_d{k3, 3*sens_var(k)-1+k_var(k2)})
                data_all_d(k,k2,k3)=10;
            else
                data_all_d(k,k2,k3)=data_d{k3, 3*sens_var(k)-1+k_var(k2)}-data_d{1, 5+k_var(k2)};
            end
        end
    end
end

% --- Solution matrices ---
end_no_d = zeros(2,length(ends_d(:,1)));
x_dat_d = 10*ones(3,3,length(ends_d(:,1)),floor(time_diff_d/increment_d)+1);
time_segmented_d = zeros(length(ends_d(:,1)),floor(time_diff_d/increment_d)+1);
mean_values_d = zeros(3,3,3,length(ends_d(:,1)));
normalised_top_bot_d = zeros(3,3,2,length(ends_d(:,1)));
ranges_d = zeros(3,3,2,length(ends_d(:,1)));
g_d = zeros(3,3,3, length(ends_d(:,1)));
ylims_d = zeros(3,3,2);

% --- Solution ---
for j1 = 1:3
    for j2 = 1:3
        for i = 1:length(ends_d(:,1))
            end_no_d(1,i)=ceil(ends_d(i,1)/increment_d);
            end_no_d(2,i)=end_no_d(1,i)+floor((ends_d(i,2)-ends_d(i,1))/increment_d);
            x_dat_d(j1,j2,i,1:floor((ends_d(i,2)-ends_d(i,1))/increment_d)+1) =  data_all_d(j1,j2,end_no_d(1,i):end_no_d(2,i));
            sortedCol = sort(x_dat_d(j1,j2,i,:));
            var_lh = [40 39];
            Lowest  = sortedCol(var_lh(1));
            Highest = sortedCol(end-var_lh(2));
            if Highest==10
                while Highest==10
                    var_lh(2)=var_lh(2)+1;
                    Highest=sortedCol(end-var_lh(2));
                end
                var_lh(2)=var_lh(2)+40;
                Highest=sortedCol(end-var_lh(2));
            end
            mask_mean = (x_dat_d(j1,j2,i,:) <= Highest) & (x_dat_d(j1,j2,i,:) >= Lowest);
            mean_values_d(j1,j2,2,i) = median(x_dat_d(j1,j2,i,mask_mean), 'omitnan');
            ranges_d(j1,j2,:,i) = [0.05*abs(Highest-mean_values_d(j1,j2,2,i)) 0.05*abs(Lowest-mean_values_d(j1,j2,2,i))];
            mask_top = (x_dat_d(j1,j2,i,:) >= Highest-ranges_d(j1,j2,1,i)) & (x_dat_d(j1,j2,i,:) <= Highest+0.5*ranges_d(j1,j2,1,i));
            mask_bot = (x_dat_d(j1,j2,i,:) <= Lowest+ranges_d(j1,j2,2,i)) & (x_dat_d(j1,j2,i,:) >= Lowest-0.5*ranges_d(j1,j2,2,i));
            x1_top = x_dat_d(j1,j2,i,mask_top);
            x1_bot = x_dat_d(j1,j2,i,mask_bot);
            mean_values_d(j1,j2,1,i) = median(x1_top, 'omitnan');
            mean_values_d(j1,j2,3,i) = median(x1_bot, 'omitnan');
            normalised_top_bot_d(j1,j2,:,i)= [mean_values_d(j1,j2,1,i)-mean_values_d(j1,j2,2,i) -(mean_values_d(j1,j2,2,i)-mean_values_d(j1,j2,3,i))];
            time_segmented_d(i,1:floor((ends_d(i,2)-ends_d(i,1))/increment_d)+1) =  time_d(end_no_d(1,i):end_no_d(2,i));
            g_d(j1,j2,:,i) = mean(time_segmented_d(i,:));
        end
        tl = max(mean_values_d(j1,j2,1,:));
        bl = min(mean_values_d(j1,j2,3,:));
        diff_l = 0.1*abs(tl-bl);
        ylims_d(j1,j2,:)=[bl-diff_l tl+diff_l];
    end
end



% --- User settings - 2nd csv ---
filename_u = 'proppel_15_up_1240_1460_07.csv';
rpm_u = zeros(1,12);
rpm_u(1,1) = 1240;
for rpm_var = 1:11
    rpm_u(1,rpm_var+1) = rpm_u(1,rpm_var)+20;
end
rpm_plot_u = [rpm_u;rpm_u];

ends_u = [30 60;75 120;135 180;195 240;255 300;315 360;375 420;
    435 480;495 540;555 600;615 660;675 720];
time_diff_u = 45;

% --- Read the CSV file, skipping header lines ---
opts_u = detectImportOptions(filename_u);
opts_u.DataLines = [headerLines+1, Inf];   % start reading after header lines
data_u = readtable(filename_u, opts_u);

% --- Extract columns ---
time_u = data_u{:, 2};
increment_u = time_u(2);
data_all_u = zeros(3,3, length(time_u));
for k = 1:3
    for k2 = 1:3
        for k3 = 1:length(time_u)
            if isnan(data_u{k3, 3*sens_var(k)-1+k_var(k2)})
                data_all_u(k,k2,k3)=10;
            else
                data_all_u(k,k2,k3)=data_u{k3, 3*sens_var(k)-1+k_var(k2)}-data_u{1, 5+k_var(k2)};
            end
        end
    end
end

% --- Solution matrices ---
end_no_u = zeros(2,length(ends_u(:,1)));
x_dat_u = 10*ones(3,3,length(ends_u(:,1)),floor(time_diff_u/increment_u)+1);
time_segmented_u = zeros(length(ends_u(:,1)),floor(time_diff_u/increment_u)+1);
mean_values_u = zeros(3,3,3,length(ends_u(:,1)));
normalised_top_bot_u = zeros(3,3,2,length(ends_u(:,1)));
ranges_u = zeros(3,3,2,length(ends_u(:,1)));
g_u = zeros(3,3,3, length(ends_u(:,1)));
ylims_u = zeros(3,3,2);

% --- Solution ---
for j1 = 1:3
    for j2 = 1:3
        for i = 1:length(ends_u(:,1))
            end_no_u(1,i)=ceil(ends_u(i,1)/increment_u);
            end_no_u(2,i)=end_no_u(1,i)+floor((ends_u(i,2)-ends_u(i,1))/increment_u);
            x_dat_u(j1,j2,i,1:floor((ends_u(i,2)-ends_u(i,1))/increment_u)+1) =  data_all_u(j1,j2,end_no_u(1,i):end_no_u(2,i));
            sortedCol = sort(x_dat_u(j1,j2,i,:));
            var_lh = [40 39];
            Lowest  = sortedCol(var_lh(1));
            Highest = sortedCol(end-var_lh(2));
            if Highest==10
                while Highest==10
                    var_lh(2)=var_lh(2)+1;
                    Highest=sortedCol(end-var_lh(2));
                end
                var_lh(2)=var_lh(2)+40;
                Highest=sortedCol(end-var_lh(2));
            end
            mask_mean = (x_dat_u(j1,j2,i,:) <= Highest) & (x_dat_u(j1,j2,i,:) >= Lowest);
            mean_values_u(j1,j2,2,i) = median(x_dat_u(j1,j2,i,mask_mean), 'omitnan');
            ranges_u(j1,j2,:,i) = [0.05*abs(Highest-mean_values_u(j1,j2,2,i)) 0.05*abs(Lowest-mean_values_u(j1,j2,2,i))];
            mask_top = (x_dat_u(j1,j2,i,:) >= Highest-ranges_u(j1,j2,1,i)) & (x_dat_u(j1,j2,i,:) <= Highest+0.5*ranges_u(j1,j2,1,i));
            mask_bot = (x_dat_u(j1,j2,i,:) <= Lowest+ranges_u(j1,j2,2,i)) & (x_dat_u(j1,j2,i,:) >= Lowest-0.5*ranges_u(j1,j2,2,i));
            x1_top = x_dat_u(j1,j2,i,mask_top);
            x1_bot = x_dat_u(j1,j2,i,mask_bot);
            mean_values_u(j1,j2,1,i) = median(x1_top, 'omitnan');
            mean_values_u(j1,j2,3,i) = median(x1_bot, 'omitnan');
            normalised_top_bot_u(j1,j2,:,i)= [mean_values_u(j1,j2,1,i)-mean_values_u(j1,j2,2,i) -(mean_values_u(j1,j2,2,i)-mean_values_u(j1,j2,3,i))];
            time_segmented_u(i,1:floor((ends_u(i,2)-ends_u(i,1))/increment_u)+1) =  time_u(end_no_u(1,i):end_no_u(2,i));
            g_u(j1,j2,:,i) = mean(time_segmented_u(i,:));
        end
        tl = max(mean_values_u(j1,j2,1,:));
        bl = min(mean_values_u(j1,j2,3,:));
        diff_l = 0.1*abs(tl-bl);
        ylims_u(j1,j2,:)=[bl-diff_l tl+diff_l];
    end
end

% --- Plot ---
% -- plot matrices --
plot_dat_d = zeros(1,length(time_d));
plot_dat_x_d = zeros(1,length(time_d));
plot_dat_yz_d = zeros(1,length(time_d));
plot_dat_seg_d = zeros(length(ends_d(:,1)),floor((ends_d(i,2)-ends_d(i,1))/increment_d)+1);
plot_time_seg_d = zeros(length(ends_d(:,1)),floor((ends_d(i,2)-ends_d(i,1))/increment_d)+1);
plot_mean_x_d = zeros(3,length(ends_d(:,1)));
plot_all_mean_d = zeros(3,length(ends_d(:,1)));
plot_bif_d = zeros(2,length(ends_d(:,1)));

plot_dat_u = zeros(1,length(time_u));
plot_dat_x_u = zeros(1,length(time_u));
plot_dat_yz_u = zeros(1,length(time_u));
plot_dat_seg_u = zeros(length(ends_u(:,1)),floor((ends_u(i,2)-ends_u(i,1))/increment_u)+1);
plot_time_seg_u = zeros(length(ends_u(:,1)),floor((ends_u(i,2)-ends_u(i,1))/increment_u)+1);
plot_mean_x_u = zeros(3,length(ends_u(:,1)));
plot_all_mean_u = zeros(3,length(ends_u(:,1)));
plot_bif_u = zeros(2,length(ends_u(:,1)));

plot_segment = 1;
titles = ["x","y","z"];

% -- plot downnward data --
figure(1)
t1 = tiledlayout(2,3);
for j1 = 1:2
    for j2 = 1:3
        nexttile;
        plot_dat_d(:) = data_all_d(j1,j2,:);
        plot(time_d,plot_dat_d,'k.');
        xlabel('Time [s]');
        ylabel('Displacement [mm]');
        title('Displacement of measurement point No.'+string(j1)+' in the '+string(titles(j2))+'-direction');
        xlim([(ends_d(1,1)-10) (ends_d(12,2)+10)])
        ylim(ylims_d(j1,j2,:))
        grid on;
        hold off
    end
end
title(t1, 'Raw displacement vs time data - down at l=15%')

% -- plot segmented downward data --
if plot_segment==1
    figure(2)
    t2 = tiledlayout(2,3);
    for j1 = 1:2
        for j2 = 1:3
            nexttile;
            plot_dat_seg_d(:,:) = x_dat_d(j1,j2,:,:);
            plot_time_seg_d(:,:) = time_segmented_d(:,:);
            plot(plot_time_seg_d,plot_dat_seg_d,'k.');
            hold on
            plot_mean_x_d(:,:) = g_d(j1,j2,:,:);
            plot_all_mean_d(:,:) = mean_values_d(j1,j2,:,:);
            scatter(plot_mean_x_d,plot_all_mean_d);
            xlabel('Time [s]');
            ylabel('Displacement [mm]');
            title('Displacement of measurement point No.'+string(j1)+' in the '+string(titles(j2))+'-direction');
            xlim([(ends_d(1,1)-10) (ends_d(12,2)+10)])
            ylim(ylims_d(j1,j2,:))
            grid on;
            hold off
        end
    end
    title(t2, 'Segmented displacement vs time diagrams with mean values - down at l=15%')
end

% -- plot upnward data --
figure(3)
t3 = tiledlayout(2,3);
for j1 = 1:2
    for j2 = 1:3
        nexttile;
        plot_dat_u(:) = data_all_u(j1,j2,:);
        plot(time_u,plot_dat_u,'k.');
        xlabel('Time [s]');
        ylabel('Displacement [mm]');
        title('Displacement of measurement point No.'+string(j1)+' in the '+string(titles(j2))+'-direction');
        xlim([(ends_u(1,1)-10) (ends_u(12,2)+10)])
        ylim(ylims_u(j1,j2,:))
        grid on;
        hold off
    end
end
title(t3, 'Raw displacement vs time data - up at l=15%')

% -- plot segmented upward data --
if plot_segment==1
    figure(4)
    t4 = tiledlayout(2,3);
    for j1 = 1:2
        for j2 = 1:3
            nexttile;
            plot_dat_seg_u(:,:) = x_dat_u(j1,j2,:,:);
            plot_time_seg_u(:,:) = time_segmented_u(:,:);
            %temp_2(:,:) = data_x2(j1,j2,:,:);
            plot(plot_time_seg_u,plot_dat_seg_u,'k.');
            hold on
            plot_mean_x_u(:,:) = g_u(j1,j2,:,:);
            plot_all_mean_u(:,:) = mean_values_u(j1,j2,:,:);
            scatter(plot_mean_x_u,plot_all_mean_u);
            xlabel('Time [s]');
            ylabel('Displacement [mm]');
            title('Displacement of measurement point No.'+string(j1)+' in the '+string(titles(j2))+'-direction');
            xlim([(ends_u(1,1)-10) (ends_u(12,2)+10)])
            ylim(ylims_u(j1,j2,:))
            grid on;
            hold off
        end
    end
    title(t4, 'Segmented displacement vs time diagrams with mean values - up at l=15%')
end

% -- plot bistable area --
figure(5)
t5 = tiledlayout(2,3);
for j1 = 1:2
    for j2 = 1:3
        nexttile;
        hold on
        plot_bif_u(:,:) = normalised_top_bot_u(j1,j2,:,:);
        plot_bif_d(:,:) = normalised_top_bot_d(j1,j2,:,:);
        plot(rpm_plot_d(1,:),plot_bif_d(1,:));
        plot(rpm_plot_d(2,:),plot_bif_d(2,:));
        plot(rpm_plot_u(1,:),plot_bif_u(1,:));
        plot(rpm_plot_u(2,:),plot_bif_u(2,:));
        xlabel('Pulse-width modulation (PWM) [-]');
        ylabel('Displacement [mm]');
        title('Displacement of measurement point No.'+string(j1)+' in the '+string(titles(j2))+'-direction');
        grid on;
        hold off
    end
end
title(t5, 'Maximum displacement values vs rpm at l=15%')
%% 

% -- plot phase downward data --
varfig_d = 10;
for j1 = 1:2
    plot_dat_x_d(:) = data_all_d(j1,1,:);
    for j2 = 2:3
        varfig_d = varfig_d+1;
        figure(varfig_d)
        plot_dat_yz_d(:) = data_all_d(j1,j2,:);
        hold on
        for j3 = 1:length(ends_d(:,1))
            plot(plot_dat_x_d(end_no_d(1,j3):end_no_d(2,j3)),plot_dat_yz_d(end_no_d(1,j3):end_no_d(2,j3)),'.');
        end
        plot(plot_dat_x_d(1:end_no_d(1,1)),plot_dat_yz_d(1:end_no_d(1,1)),'k.','MarkerSize',1);
        for j4 = 1:length(ends_d(:,1))-1
            plot(plot_dat_x_d(end_no_d(2,j4):end_no_d(1,j4+1)),plot_dat_yz_d(end_no_d(2,j4):end_no_d(1,j4+1)),'k.','MarkerSize',1);
        end
        leg = legend(string(rpm_d(1)),string(rpm_d(2)),string(rpm_d(3)),string(rpm_d(4)),string(rpm_d(5)), ...
            string(rpm_d(6)),string(rpm_d(7)),string(rpm_d(8)),string(rpm_d(9)),string(rpm_d(10)),string(rpm_d(11)),string(rpm_d(12)));
        title(leg,'PWM values')
        xlabel('Displacement in x-direction [mm]');
        ylabel('Displacement in '+string(titles(j2))+'-direction [mm]');
        title('Phase diagram of measurement point No.'+string(j1)+' in the x'+string(titles(j2))+'-plane - down at l=15%')
        xlim(ylims_d(j1,1,:))
        ylim(ylims_d(j1,j2,:))
        grid on;
        hold off
    end
    varfig_d = varfig_d+1;
    figure(varfig_d)
    plot_dat_x_d(:) = data_all_d(j1,3,:);
    plot_dat_yz_d(:) = data_all_d(j1,2,:);
    hold on
    for j3 = 1:length(ends_d(:,1))
        plot(plot_dat_x_d(end_no_d(1,j3):end_no_d(2,j3)),plot_dat_yz_d(end_no_d(1,j3):end_no_d(2,j3)),'.');
    end
    plot(plot_dat_x_d(1:end_no_d(1,1)),plot_dat_yz_d(1:end_no_d(1,1)),'k.','MarkerSize',1);
    for j4 = 1:length(ends_d(:,1))-1
        plot(plot_dat_x_d(end_no_d(2,j4):end_no_d(1,j4+1)),plot_dat_yz_d(end_no_d(2,j4):end_no_d(1,j4+1)),'k.','MarkerSize',1);
    end
    leg = legend(string(rpm_d(1)),string(rpm_d(2)),string(rpm_d(3)),string(rpm_d(4)),string(rpm_d(5)), ...
        string(rpm_d(6)),string(rpm_d(7)),string(rpm_d(8)),string(rpm_d(9)),string(rpm_d(10)),string(rpm_d(11)),string(rpm_d(12)));
    title(leg,'PWM values')
    xlabel('Displacement in z-direction [mm]');
    ylabel('Displacement in y-direction [mm]');
    title('Phase diagram of measurement point No.'+string(j1)+' in the zy-plane - down at l=15%')
    xlim(ylims_d(j1,3,:))
    ylim(ylims_d(j1,2,:))
    grid on;
    hold off
end


% -- plot phase upward data --
varfig_u = 20;
for j1 = 1:2
    plot_dat_x_u(:) = data_all_u(j1,1,:);
    for j2 = 2:3
        varfig_u = varfig_u+1;
        figure(varfig_u)
        plot_dat_yz_u(:) = data_all_u(j1,j2,:);
        hold on
        for j3 = 1:length(ends_u(:,1))
            plot(plot_dat_x_u(end_no_u(1,j3):end_no_u(2,j3)),plot_dat_yz_u(end_no_u(1,j3):end_no_u(2,j3)),'.');
        end
        plot(plot_dat_x_u(1:end_no_u(1,1)),plot_dat_yz_u(1:end_no_u(1,1)),'k.','MarkerSize',1);
        for j4 = 1:length(ends_u(:,1))-1
            plot(plot_dat_x_u(end_no_u(2,j4):end_no_u(1,j4+1)),plot_dat_yz_u(end_no_u(2,j4):end_no_u(1,j4+1)),'k.','MarkerSize',1);
        end
        leg = legend(string(rpm_u(1)),string(rpm_u(2)),string(rpm_u(3)),string(rpm_u(4)),string(rpm_u(5)), ...
            string(rpm_u(6)),string(rpm_u(7)),string(rpm_u(8)),string(rpm_u(9)),string(rpm_u(10)),string(rpm_u(11)),string(rpm_u(12)));
        title(leg,'PWM values')
        xlabel('Displacement in x-direction [mm]');
        ylabel('Displacement in '+string(titles(j2))+'-direction [mm]');
        title('Phase diagram of measurement point No.'+string(j1)+' in the x'+string(titles(j2))+'-plane - up at l=15%')
        xlim(ylims_u(j1,1,:))
        ylim(ylims_u(j1,j2,:))
        grid on;
        hold off
    end
    varfig_u = varfig_u+1;
    figure(varfig_u)
    plot_dat_x_u(:) = data_all_u(j1,3,:);
    plot_dat_yz_u(:) = data_all_u(j1,2,:);
    hold on
    for j3 = 1:length(ends_u(:,1))
        plot(plot_dat_x_u(end_no_u(1,j3):end_no_u(2,j3)),plot_dat_yz_u(end_no_u(1,j3):end_no_u(2,j3)),'.');
    end
    plot(plot_dat_x_u(1:end_no_u(1,1)),plot_dat_yz_u(1:end_no_u(1,1)),'k.','MarkerSize',1);
    for j4 = 1:length(ends_u(:,1))-1
        plot(plot_dat_x_u(end_no_u(2,j4):end_no_u(1,j4+1)),plot_dat_yz_u(end_no_u(2,j4):end_no_u(1,j4+1)),'k.','MarkerSize',1);
    end
    leg = legend(string(rpm_u(1)),string(rpm_u(2)),string(rpm_u(3)),string(rpm_u(4)),string(rpm_u(5)), ...
        string(rpm_u(6)),string(rpm_u(7)),string(rpm_u(8)),string(rpm_u(9)),string(rpm_u(10)),string(rpm_u(11)),string(rpm_u(12)));
    title(leg,'PWM values')
    xlabel('Displacement in z-direction [mm]');
    ylabel('Displacement in y-direction [mm]');
    title('Phase diagram of measurement point No.'+string(j1)+' in the zy-plane - up at l=15%')
    xlim(ylims_u(j1,3,:))
    ylim(ylims_u(j1,2,:))
    grid on;
    hold off
end
