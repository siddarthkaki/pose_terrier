%% housekeeping
clear; close all; clc;

%% init

OPT_misc = 1;


data_dir = "../data/sfm_2023_test2/";

run_name = "1672605165";

ptl_data = f_read_filter_logs(data_dir + run_name + "_");
ptl_data.SolvedAngVel = rad2deg(f_read_omegas(data_dir + run_name + "_" + "solved_omegas.csv"));

tru_data = f_read_truth_logs(data_dir);

%% angle wrap-around processing -- per-case basis!
% ddx = 2;
% 
% tru_data.TrueRotation(:,ddx) = wrapTo360(tru_data.TrueRotation(:,ddx));
% 
% ptl_data.SolvedRotation(:,ddx) = wrapTo360(ptl_data.SolvedRotation(:,ddx));
% ptl_data.FilteredRotation(:,ddx) = wrapTo360(ptl_data.FilteredRotation(:,ddx));



%% interpolation
for idx = 1:3,
    v = tru_data.TrueRotation(:,idx);
    vq1 = interp1(tru_data.Time(~isnan(v)), v(~isnan(v)), ptl_data.Time,'linear');
    ptl_data.TrueRotation(:,idx) = vq1;

    v = tru_data.TrueTranslation(:,idx);
    vq1 = interp1(tru_data.Time(~isnan(v)), v(~isnan(v)), ptl_data.Time,'linear');
    ptl_data.TrueTranslation(:,idx) = vq1;

    v = tru_data.TrueAngVel(:,idx);
    vq1 = interp1(tru_data.Time(~isnan(v)), v(~isnan(v)), ptl_data.Time,'linear');
    ptl_data.TrueAngVel(:,idx) = vq1;
end


%% interpolated errors
% iidx = find(ptl_data.Time > 20.0519 & ptl_data.Time <= 20.3802);
% ptl_data.FilteredRotation(iidx,2) = ptl_data.FilteredRotation(iidx,2) + 360;

c_err_est_mat = NaN(length(ptl_data.Time),4);

for tdx = 1:length(ptl_data.Time),

    true_rot = ptl_data.TrueRotation(tdx,:);
    filt_rot = ptl_data.FilteredRotation(tdx,:);
    unfl_rot = ptl_data.SolvedRotation(tdx,:);

    true_quat = dcm2quat(rotation.euler2dcm_312(deg2rad(true_rot)));
    filt_quat = dcm2quat(rotation.euler2dcm_312(deg2rad(filt_rot)));
    unfl_quat = dcm2quat(rotation.euler2dcm_312(deg2rad(unfl_rot)));

    dquat_filt = rotation.quatmult_S(filt_quat, quatinv(true_quat));
    ptl_data.ErrorRotation(tdx) = rad2deg( 2*acos( abs( dquat_filt(1) ) ) );
    ptl_data.ErrorGibbs(tdx,:) = dquat_filt(2:4)/dquat_filt(1);

    dquat_unfl = rotation.quatmult_S(unfl_quat, quatinv(true_quat));
    ptl_data.ErrorRotationNLS(tdx) = rad2deg( 2*acos( abs( dquat_unfl(1) ) ) );
    ptl_data.ErrorGibbsNLS(tdx,:) = dquat_unfl(2:4)/dquat_unfl(1);

    %%%

    true_trans = ptl_data.TrueTranslation(tdx,:);
    filt_trans = ptl_data.FilteredTranslation(tdx,:);
    unfl_trans = ptl_data.SolvedTranslation(tdx,:);

    ptl_data.NormErrorTranslation(tdx) = norm(filt_trans - true_trans);
    ptl_data.NormErrorTranslationNormed(tdx) = norm(filt_trans - true_trans)/norm(true_trans);

    ptl_data.NormErrorTranslationNLS(tdx) = norm(unfl_trans - true_trans);
    ptl_data.NormErrorTranslationNormedNLS(tdx) = norm(unfl_trans - true_trans)/norm(true_trans);

    %%%

    true_angvel = ptl_data.TrueAngVel(tdx,:);
    filt_angvel = ptl_data.FilteredAngVel(tdx,:);

    ptl_data.ErrorAngVel(tdx,:) = filt_angvel - true_angvel;
    
    %%%
    a = true_quat.';
    c = filt_quat;
    b = rotation.quatmult_S(c, quatinv(a.').');

    c_err_est(1,1) = -b(2:4).'*a(2:4);
    c_err_est(2,1) = a(1)*b(2) + a(3)*b(4) - a(4)*b(3);
    c_err_est(3,1) = a(1)*b(3) - a(2)*b(4) + a(4)*b(2);
    c_err_est(4,1) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2);

    c_err_est_mat(tdx,:) = c_err_est;
end


%%

% iidx = find(ptl_data.Time <= 20.0519 | ptl_data.Time > 20.3802);
iidx = find(ptl_data.Time <= 11.4901 | ...
           (ptl_data.Time >= 11.7615 & ptl_data.Time <= 50.0415) | ...
           (ptl_data.Time >= 50.2918));
           %(ptl_data.Time >= 37.7801 & ptl_data.Time <= 45.7921) | ...
%iidx = find(ptl_data.Time <= 20.0402 | ptl_data.Time >= 20.3102);

%%
f1111 = figure(1111);
for idx = 1:4,
    c_err_est_mat_denan(:,idx) = denan(c_err_est_mat(:,idx));
    r = c_err_est_mat(:,idx);
    pd = fitdist(r,'Normal');
    subplot(4,1,idx)
    histfit(r,200)
    grid on
    ylabel("\Delta q_{\eta}(" + string(idx-1) + ")")
    set(gca,'FontWeight','bold')

    %if idx == 1, xlim(1.8*[-1 1]*10^-3); end
    %if idx ~= 1, xlim(0.08*[-1 1]); end
end
fontsize(f1111,scale=2)
set(gcf, 'PaperPositionMode', 'auto');

dqm = mean(c_err_est_mat_denan).'
dqs = std(c_err_est_mat_denan).'


%%
data_table = f_statistics_quatera(ptl_data,iidx)

% writematrix(mobilepose_data.ErrorRotationInterp.','csvs/heavy_unfiltered_ER.csv')
% writematrix(ptl_data.ErrorRotation.','csvs/heavy_filtered_ER.csv')
% writematrix(mobilepose_data.NormErrorTranslationInterp.','csvs/heavy_unfiltered_ET.csv')
% writematrix(ptl_data.NormErrorTranslation.','csvs/heavy_filtered_ET.csv')

%% rotation
f1 = figure(1);

for idx = 1:3,
    subplot(4,1,idx)
    hold on

    %         fl = plot(tru_data.Time, tru_data.PredictedRotation(:,idx),'-','DisplayName','filtered');
    %         tr = plot(tru_data.Time, tru_data.TrueRotation(:,idx),'k--','DisplayName','truth');
    fl = plot(ptl_data.Time, ptl_data.FilteredRotation(:,idx),'Color','#D95319','DisplayName','filtered');
    tr = plot(ptl_data.Time, ptl_data.TrueRotation(:,idx),'k--','DisplayName','truth');
    plot(ptl_data.Time, ptl_data.SolvedRotation(:,idx),'Color','#0072BD','DisplayName','NLS');

    uistack(fl,'top')
    uistack(tr,'top')

    hold off
    grid on

    switch (idx)
        case 1
            ylabel('\phi [deg]')
            legend('Location','best')
            %             title('Rotation')
            xlim([0 60]); xticks(0:10:60)
            ylim([-90 95]); yticks(-90:45:90); %%%%%

        case 2
            xlim([0 60]); xticks(0:10:60)
            ylabel('\theta [deg]')
            ylim([-180 180]); yticks(-180:90:180); %%%%%

        case 3
            xlim([0 60]); xticks(0:10:60)
            ylabel('\psi [deg]')
            ylim([-180 180]); yticks(-180:90:180); %%%%%
    end
end

subplot(4,1,4)
hold on
fl = plot(ptl_data.Time(iidx), ptl_data.ErrorRotation(iidx),'Color','#D95319','DisplayName','filtered');
plot(ptl_data.Time(iidx), ptl_data.ErrorRotationNLS(iidx),'Color','#0072BD','DisplayName','NLS');
uistack(fl,'top')
hold off
grid on
xlim([0 60]); xticks(0:10:60)
% ylim([0 35]); yticks(0:10:30)
ylim([0 100]); yticks(0:25:100)
% ylim([0 180]); yticks(0::30);
% ylim([0 180]); yticks(0:45:180); %%%%%
ylabel('error [deg]')
xlabel('time [s]')
% legend('Location','best')

boldify;
set(gcf, 'PaperPositionMode', 'auto');



% %% translation
% f2 = figure(2);
% 
% for idx = 1:3,
%     subplot(4,1,idx)
%     hold on
%     %         fl = plot(tru_data.Time, tru_data.PredictedTranslation(:,idx),'-','DisplayName','filtered');
%     %         tr = plot(tru_data.Time, tru_data.TrueTranslation(:,idx),'k--','DisplayName','truth');
%     fl = plot(ptl_data.Time, ptl_data.FilteredTranslation(:,idx),'-','DisplayName','filtered');
%     tr = plot(ptl_data.Time, ptl_data.TrueTranslation(:,idx),'k--','DisplayName','truth');
%     %         plot(ptl_data.Time, ptl_data.SolvedTranslation(:,idx),'DisplayName','ptl nls');
% 
%     uistack(fl,'top')
%     uistack(tr,'top')
% 
% 
%     hold off
%     grid on
% 
%     switch (idx)
%         case 1
%             ylabel('x [m]')
%             legend('Location','best')
%             %             title('Translation')
% %             ylim([-15 5]);
%             xlim([0 60]); xticks(0:10:60)
%             %             yticks(-15:5:5); %%%%%
% 
%         case 2
%             xlim([0 60]); xticks(0:10:60)
%             ylabel('y [m]')
% %             ylim([-10 10]);
%             yticks(-10:5:10); %%%%%
% 
%         case 3
%             xlim([0 60]); xticks(0:10:60)
%             ylabel('z [m]')
% %             ylim([40 70]);
%             yticks(40:10:70); %%%%%
% 
%     end
% end
% 
% subplot(4,1,4)
% hold on
% % fl = plot(tru_data.Time, tru_data.NormErrorTranslation,'-','DisplayName','filtered');
% fl = plot(ptl_data.Time, ptl_data.NormErrorTranslation,'-','DisplayName','filtered');
% uistack(fl,'top')
% hold off
% grid on
% % ylim([0 0])
% xlim([0 60]); xticks(0:10:60)
% % ylim([0 4]); yticks(0:1:4); %%%%%
% ylabel('error [m]')
% xlabel('time [s]')
% % legend('Location','best')
% 
% boldify;
% set(gcf, 'PaperPositionMode', 'auto');


%% rotation 3sigma - separate
f111 = figure(111);
colourvecd = ["#7E2F8E","#e2a611","#419ac0"]; % dark
colourvecl = ["#d078e2","#ffd879","#4DBEEE"]; % light
displaynames = ["g_1","g_2","g_3"];

for idx = 1:3,
    subplot(3,1,idx)
    hold on
    plot(ptl_data.Time,  3*ptl_data.FilteredStd(:,idx),'Color','black','DisplayName',displaynames(idx) + ": \pm3\sigma");
    pp = plot(ptl_data.Time, -3*ptl_data.FilteredStd(:,idx),'Color','black');
    set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(ptl_data.Time(iidx), ptl_data.ErrorGibbs(iidx,idx),'-','Color','#419ac0','DisplayName',displaynames(idx) + ": error");
    hold off
    grid on
    ylabel("g_" + idx + " error")
    xlim([0 60]); xticks(0:10:60)
    ylim(0.4*[-1 1]); yticks(-0.4:0.2:0.4)

    if idx == 1, legend("\pm3\sigma","error",'Location','best'); end
end
grid on
xlabel('time [s]')
xlim([0 60]); xticks(0:10:60)
% legend('Location','best')

boldify();
set(gcf, 'PaperPositionMode', 'auto');



%% ang vel 3sigma - separate
% f211 = figure(211);
% colourvecd = ["#7E2F8E","#e2a611","#419ac0"]; % dark
% colourvecl = ["#d078e2","#ffd879","#4DBEEE"]; % light
% displaynames = ["\phi","\theta","\psi"];
% 
% for idx = 1:3,
%     subplot(3,1,idx)
%     hold on
%     plot(ptl_data.Time,  3*rad2deg(ptl_data.FilteredStd(:,idx+3)),'Color','black','DisplayName',displaynames(idx) + ": \pm3\sigma");
%     pp = plot(ptl_data.Time, -3*rad2deg(ptl_data.FilteredStd(:,idx+3)),'Color','black');
%     set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     plot(ptl_data.Time, ptl_data.ErrorAngVel(:,idx),'-','Color','#419ac0','DisplayName',displaynames(idx) + ": error");
%     hold off
%     grid on
%     ylabel("deg/s")
%     xlim([0 60]); xticks(0:10:60)
%     ylim(10*[-1 1]);
% %     yticks(-0.5:0.1:0.5)
% 
%     if idx == 1, legend("\pm3\sigma","error",'Location','best'); end
% end
% grid on
% xlabel('time [s]')
% xlim([0 60]); xticks(0:10:60)
% % legend('Location','best')
% 
% boldify();
% set(gcf, 'PaperPositionMode', 'auto');


%% translation 3sigma - separate
% f121 = figure(121);
% 
% displaynames = ["x","y","z"];
% 
% for idx = 1:3,
%     subplot(3,1,idx)
%     hold on
%     plot(ptl_data.Time,  3*ptl_data.FilteredStd(:,idx+9),'Color','black','DisplayName',displaynames(idx) + ": \pm3\sigma");
%     pp = plot(ptl_data.Time, -3*ptl_data.FilteredStd(:,idx+9),'Color','black');
%     set(get(get(pp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     plot(ptl_data.Time, ptl_data.FilteredTranslation(:,idx) - ptl_data.TrueTranslation(:,idx),'-','Color','#419ac0','DisplayName',displaynames(idx) + ": error");
%     hold off
%     grid on
%     ylabel(displaynames(idx) + " error [m]")
%     xlim([0 60]); xticks(0:10:60)
% 
%     switch (idx)
%         case 1
%             legend("\pm3\sigma","error",'Location','best');
% %             ylim(16*[-1 1])
%             %             yticks(-15:5:5); %%%%%
% 
%         case 2
% %             ylim(16*[-1 1])
%             %             yticks(-10:5:10); %%%%%
% 
%         case 3
% %             ylim(65*[-1 1]);
%             yticks(-60:30:60); %%%%%
%     end
% 
% end
% xlabel('time [s]')
% xlim([0 60]); xticks(0:10:60)
% % legend('Location','best')
% 
% boldify();
% set(gcf, 'PaperPositionMode', 'auto');


%%

for idx = 1:length(tru_data.TrueTranslation(:,1)),
    distances(idx) = norm(tru_data.TrueTranslation(idx,:));
end

% %% 3D trajectory plot
% figure(100)
% % plot3(tru_data.TrueTranslation(1,1), tru_data.TrueTranslation(1,2), tru_data.TrueTranslation(1,3),'o')
% % plot3(tru_data.TrueTranslation(end,1), tru_data.TrueTranslation(end,2), tru_data.TrueTranslation(end,3),'x')
% h = scatter3(tru_data.TrueTranslation(:,1),tru_data.TrueTranslation(:,2),tru_data.TrueTranslation(:,3),40,tru_data.Time);
% hold on
% plot3(tru_data.PredictedTranslation(:,1), tru_data.PredictedTranslation(:,2), tru_data.PredictedTranslation(:,3))
% % plot3(ptl_data.FilteredTranslation(:,1), ptl_data.FilteredTranslation(:,2), ptl_data.FilteredTranslation(:,3))
% hold off
% 
% % h.MarkerFaceColor = 'flat';
% colormap(jet)
% cb = colorbar;
% cb.Label.String = 'time [s]';
% % cb.Ticks = [0:10:tru_data.Time(end)];
% 
% hold off
% axis equal
% grid on
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% % legend('relative trajectory','start','stop')
% 
% boldify;
% set(gcf, 'PaperPositionMode', 'auto');



%% misc
if OPT_misc,

    % att 3sigmas
    f101 = figure(101);
    subplot(3,1,1)
    plot(ptl_data.Time, 3*rad2deg(ptl_data.FilteredStd(:,1:3)));
    grid on
    ylabel('deg')
    ylim([0 15])
    % title('3\sigma')

    subplot(3,1,2)
    plot(ptl_data.Time, 3*rad2deg(ptl_data.FilteredStd(:,4:6)));
    grid on
    ylabel('deg/s')
    ylim([0 10])

    subplot(3,1,3)
    plot(ptl_data.Time, 3*rad2deg(ptl_data.FilteredStd(:,7:9)));
    grid on
    ylabel('deg/s^2')
    ylim([0 10])

    legend('\phi', '\theta', '\psi')

    boldify;
    set(gcf, 'PaperPositionMode', 'auto');


    % plotting angular velocity

    f3 = figure(3);
    for idx = 1:3,
        subplot(3,1,idx)
        hold on
        fl = plot(ptl_data.Time, ptl_data.FilteredAngVel(:,idx),'Color','#D95319','DisplayName','filtered');
        tr = plot(ptl_data.Time, ptl_data.TrueAngVel(:,idx),'k--','DisplayName','truth');
        plot(ptl_data.Time, ptl_data.SolvedAngVel(:,idx),'Color','#0072BD','DisplayName','quatera');
        hold off
        uistack(fl,'top')
        uistack(tr,'top')
        grid on
        xlabel('time [s]')
        xlim([0 60]); xticks(0:10:60)
        ylabel('[deg/s]')
        switch (idx)
            case 1
                legend('Location','best')
                ylim(4+[-25 25]); yticks([-20:15:30])
                %title('Angular Velocity');
            case 2
                ylim(-5+[-25 25]); yticks([-30:15:20])
            case 3
                ylim(-3+[-25 25]); yticks([-25:15:30])
        end
    end

    % ylim([0 10])

    boldify
    set(gcf, 'PaperPositionMode', 'auto');

    %% plotting angular acceleration

    f13 = figure(13);
    for idx = 1:3,
        subplot(3,1,idx)
        plot(ptl_data.Time, ptl_data.FilteredAngAcc(:,idx),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
        hold off
        grid on
        xlabel('time [s]')
        ylabel('deg/s^2')
        if idx == 1, title('Angular Acceleration'); end
    end

    % ylim([0 10])
    legend('filtered')

    boldify
    set(gcf, 'PaperPositionMode', 'auto');

%     %% plotting translational velocity
% 
%     f4 = figure(4);
%     for idx = 1:3,
%         subplot(3,1,idx)
%         plot(ptl_data.Time, ptl_data.FilteredVelocity(:,idx),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
%         hold off
%         grid on
%         xlabel('time [s]')
%         ylabel('m')
%         if idx == 1, title('Translational Velocity'); end
%     end
% 
%     % ylim([0 10])
%     legend('filtered')
% 
%     boldify
%     set(gcf, 'PaperPositionMode', 'auto');

    %% plotting translational acceleration

%     f14 = figure(14);
%     for idx = 1:3,
%         subplot(3,1,idx)
%         plot(ptl_data.Time, ptl_data.FilteredAcceleration(:,idx),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
%         hold off
%         grid on
%         xlabel('time [s]')
%         ylabel('m/s^2')
%         if idx == 1, title('Translational Acceleration'); end
%     end
% 
%     % ylim([0 10])
%     legend('filtered')
% 
%     boldify
%     set(gcf, 'PaperPositionMode', 'auto');

end
