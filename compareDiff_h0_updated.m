clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

tic

%%
disjPress_switch = 'off_on';
axisSymmetry_switch = 'on';

h0_init = [300e-09];                 % initial film height in m
A_vw = 1.5e-20;                  % Hamaker constant
gam = 0.0445;                      % surface tension
Rc = 1.8e-3; 
kappa = pi*h0_init^3*gam/A_vw/Rc;               % dimensionless curvature - the free parameter of the system

mk_postProcess = strcat('axisSymmetry','h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), '_disjPr_',disjPress_switch);
mkdir(mk_postProcess);
cd(mk_postProcess);
copyfile('..\*.m', '.')
copyfile('../dataRadoev1984.xlsx', '.')
% copyfile('../getTimesData.m', '.')
% copyfile('../getCrComparisonData.m', '.')
% copyfile('../Manev1984.m', '.')


%% first comparing ratio of thinning rates with Cr

%{
relFolder = 'axisSymmetricFilm\disjPress_off\repulsion_off\Coarse_dt_h0_2.25500nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_100_25nm_c1_0_c2_0';
% relFolder = 'h0_400nm_Avw_1.25e-21_ST_0.034_Rc_0.01_disjPr_off';
filmSizes = {'20', '50', '100', '150'};
JoyeStartCrit = 0.65;
if isequal(disjPress_switch, 'on') 
    JoyeStopCrit = 0.4*kappa^(-2/7);            % keep it lower to be able to probe even smaller thicknesses if it reaches 
elseif isequal(disjPress_switch, 'off') 
    JoyeStopCrit = 0.8*kappa^(-2/7);             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
end
for i = 1:length(filmSizes)
    [c_r_Joye_r_comb{i} ratio_v_vre_comb{i} ratio_vc_vre_comb{i} c_r_Joye_centre_comb{i}] = getCrComparisonData(relFolder, filmSizes{i}, JoyeStartCrit, JoyeStopCrit, mk_postProcess);
end

hfig3 = figure;
hfig3.Renderer = 'Painters';
ColOrd = get(gca, 'colororder');
for i = 1:length(c_r_Joye_centre_comb)

    loglog(c_r_Joye_r_comb{i}, ratio_v_vre_comb{i}, 'o', 'color', ColOrd(i,:))
    hold on
    loglog(c_r_Joye_centre_comb{i}, ratio_vc_vre_comb{i},'--', 'color', ColOrd(i,:))
    hold on
end


% ylim([0.001 2])
xlabel('$C_r$')
ylabel('$v/v_{re}$')

set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');

lgd = legend(strcat('$R_f =$',filmSizes{1},'$\mu$','m'),strcat('$R_f =$',filmSizes{1},'$\mu$','m'),...
            strcat('$R_f =$',filmSizes{2},'$\mu$','m'),strcat('$R_f =$',filmSizes{2},'$\mu$','m'),...
            strcat('$R_f =$',filmSizes{3},'$\mu$','m'),strcat('$R_f =$',filmSizes{3},'$\mu$','m'), ...
            strcat('$R_f =$',filmSizes{4},'$\mu$','m'),strcat('$R_f =$',filmSizes{4},'$\mu$','m'), ...
            'interpreter','latex');
lgd.FontSize = 12;

set(lgd, 'location', 'northeastoutside');
set(lgd, 'unit', 'inches');
legend_size = get(lgd, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size);
set(gca,'FontSize',16)

fileNameSaved_01 = strcat('Cr02');
set(hfig3,'Units','Inches');
pos = get(hfig3,'Position');
set(hfig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig3,fileNameSaved_01,'-dpdf','-r300')

save('Cr01.mat')

%}
% cd ..

%%
% compare rupture times



relFolder = {'axisSymmetricFilm\disjPress_off\repulsion_off\Coarse_dt_h0_2.25300nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_100_25nm_c1_0_c2_0',...
    'axisSymmetricFilm\disjPress_on\repulsion_off\Coarse_dt_h0_300nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_100_25nm_c1_0_c2_0',...
    'axisSymmetricFilm\disjPress_on\repulsion_adhoc_on\RCoarse_dt_h0_300nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_100_25nm_c1_0_c2_0_hmin_5nm' ...
    'axisSymmetricFilm\disjPress_on\repulsion_edl_on\RCoarse_dt_h0_300nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_100_25nm_c1_80_c2_50'};
h_dimensionless = [0.05:0.05:1]';  % dimensionless film thickness vector --> will be used in the calculation of theoretical velocities
for j = 1:length(relFolder)
    pathFol{j} = strcat('..\..\..\finiteSizeStochasticThinFilms\',relFolder{j});
    cd(pathFol{j})
    load('results_differentFilmSize.mat')
    % we take the opportunity to calculate the theoretical values for the conditions prescribed for different h0 --> kappa
    [v_re_det{j} t_re{j} t_re_withoutvdW{j} v_MTR{j} t_MTR{j} t_MTR_withoutvdW{j} v_MTR_1997Paper{j} v_MTR_Tsekov{j}] = Reynolds_and_MTR(h_dimensionless, kappa, L_flat, R_f, h0_init,...
                                                                    t_scale, h_drain_start, h_drain_end, visc, gam, Rc, A_vw);

    
    h0_film{j} = h0_init*1e9;
    R_film_diffh0{j} = R_film;
    t_rupture{j} = t_rupt*t_scale;
    v_thin_min_det{j} = abs(avg_cr_thinningRate_fit).*h0_init*10^10./t_scale;
    t_drain_r{j} = drainageTime_right*t_scale;
    t_drainFull{j} = drainageTime*t_scale;
    drainageTime_Final{j} = [t_drainFull{j}(1:2) t_drain_r{j}(3:end)];
    cd(strcat('../../../../../postProcess/postProcessFiniteSizedFilms/',mk_postProcess))
end
% 
h1a = figure;
h1a.Renderer = 'Painters';
figureName_tr = strcat('ruptureTimesDifferentFilmsThicknesses_',num2str(h0_film{1}),'to',num2str(h0_film{end}), ' nm');

for i = 1:length(R_film_diffh0)
    loglog(R_film_diffh0{i}, t_rupture{i},'o')
    Legend_tr{i} = strcat('$h_0 = $', num2str(h0_film{i}),' nm');
    hold on
end
% ylim([10 80])

xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t_r$ (s)','Fontsize',14)
set(gca,'FontSize',16)

set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');

Legend_tr = {'drain.', 'drain. + vdW', 'drain. + vdW + repulsion1', 'drain. + vdW + repulsion2'};
h_legend_tr = legend(Legend_tr, 'interpreter','latex','location','best' ); 

set(h_legend_tr, 'location', 'northeastoutside');
set(h_legend_tr, 'unit', 'inches');
legend_size = get(h_legend_tr, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size);

set(h1a,'Units','Inches');
pos = get(h1a,'Position');
set(h1a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1a,figureName_tr,'-dpdf','-r300')

%% 

% first get Manev 1984 data sets
[R_manev t_drain_1stData std_t_drain_1stData t_drain_2ndData std_t_drain_2ndData] = Manev1984();
                                                                
h1b = figure;
h1b.Renderer = 'Painters';
figureName_dr = strcat('drainageTimesDifferentFilmsThicknesses_',num2str(h0_film{1}),'to',num2str(h0_film{end}), ' nm');
ColOrd = get(gca, 'colororder');
ColOrd([8:14],:) = parula(7);  
for i = 1:length(R_film_diffh0)
    p1(i) = loglog(R_film_diffh0{i},  drainageTime_Final{i},'o', 'color', ColOrd(i,:));
    hold on
end
hold on
p2 = loglog(R_film_diffh0{length(R_film_diffh0)-1}, t_re{length(R_film_diffh0)-1}, '--', 'color', ColOrd(length(R_film_diffh0) + 1,:));
hold on
p3 = loglog(R_film_diffh0{length(R_film_diffh0)-1}, t_MTR{length(R_film_diffh0)-1},'-.', 'color', ColOrd(length(R_film_diffh0) + 1,:));
hold on
p2a = loglog(R_film_diffh0{length(R_film_diffh0)-1}, t_re_withoutvdW{length(R_film_diffh0)-1}, '--', 'color', ColOrd(length(R_film_diffh0) + 4,:));
hold on
p3a = loglog(R_film_diffh0{length(R_film_diffh0)-1}, t_MTR_withoutvdW{length(R_film_diffh0)-1},'-.', 'color', ColOrd(length(R_film_diffh0) + 4,:));
hold on
p4 = errorbar(R_manev, t_drain_1stData, std_t_drain_1stData, 'd', 'color', ColOrd(length(R_film_diffh0) + 2,:));
set(gca, 'XScale','log');
set(gca, 'YScale','log');
hold on
p5 = errorbar(R_manev, t_drain_2ndData, std_t_drain_2ndData, 'd', 'color', ColOrd(length(R_film_diffh0) + 3,:));

set(gca, 'XScale','log')
set(gca, 'YScale','log')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t_{drain}$ (s)','Fontsize',14)
set(gca,'FontSize',16)
ylim([0.5*min( drainageTime_Final{1}) 1100])% 5*max(drainageTime_Final{end})])
xlim([10 1100])

set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');


Legend_dr = Legend_tr;
Legend_dr{length(Legend_tr) + 1} = 'Manev data 1 (drain.+vdW)';
Legend_dr{length(Legend_tr) + 2} = 'Manev data 2 (drain.+vdW)';
Legend_dr{length(Legend_tr) + 3} = 'Reynolds theory (drain.+vdW)';
Legend_dr{length(Legend_tr) + 4} = 'MTR theory (drain.+vdW)';
Legend_dr{length(Legend_tr) + 5} = 'Reynolds theory (drain.)';
Legend_dr{length(Legend_tr) + 6} = 'MTR theory (drain.)';
h_legend = legend([p1 p4 p5 p2(1) p3(1) p2a(1) p3a(1)], Legend_dr, 'interpreter','latex'); 
set(h_legend, 'location', 'northeastoutside');
set(h_legend, 'unit', 'inches');
legend_size = get(h_legend, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size);

set(h1b,'Units','Inches');
pos = get(h1b,'Position');
set(h1b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h1b,figureName_dr,'-dpdf','-r300');

save(strcat('drainageTimes','.mat'));

%% compare different thinning velocities at critical thickness

% read Radoev data
[L_film_Radoev L_film_Radoev_sort v_thin_Radoev_sort std_v_thin_Radoev h_cr_Radoev_sort std_h_cr_Radoev] = readRadoevData();

h2b = figure;
h2b.Renderer = 'Painters';
figureName_v_cr = strcat('criticalVel_',num2str(h0_film{1}),'to',num2str(h0_film{end}), ' nm');
ColOrd = get(gca, 'colororder');
ColOrd([8:14],:) = parula(7);  
for i = 1:length(R_film_diffh0)
    p1(i) = loglog(R_film_diffh0{i}, v_thin_min_det{i},'o', 'color', ColOrd(i,:));
    hold on
end
hold on
p2 = loglog(R_film_diffh0{length(R_film_diffh0)-1}, v_re_det{length(R_film_diffh0)-1}(2,:), '--', 'color', ColOrd(length(R_film_diffh0) + 1,:));
hold on
p3 = loglog(R_film_diffh0{length(R_film_diffh0)-1}, v_MTR{length(R_film_diffh0)-1}(:,2),'-.', 'color', ColOrd(length(R_film_diffh0) + 1,:));
hold on
p4 = errorbar(L_film_Radoev_sort*10^6, v_thin_Radoev_sort, std_v_thin_Radoev, 'd', 'color', ColOrd(length(R_film_diffh0) + 2,:));
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
ylim([1 100])



xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$v_{cr}$ (s)','Fontsize',14)
set(gca,'FontSize',16)
% ylim([0.5*min(t_drain_r{1}) 1100])% 5*max(t_drain_r{end})])
xlim([10 1100])

set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');


Legend_v_cr = Legend_tr;
Legend_v_cr{length(Legend_tr) + 1} = 'Reynolds theory (drain.+vdW)';
Legend_v_cr{length(Legend_tr) + 2} = 'MTR theory (drain.+vdW)';
Legend_v_cr{length(Legend_tr) + 3} = 'Radoev data (drain + vdW + repulsion)';
h_legend = legend([p1 p2(1) p3(1) p4], Legend_v_cr, 'interpreter','latex'); 
set(h_legend, 'location', 'northeastoutside');
set(h_legend, 'unit', 'inches');
legend_size = get(h_legend, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size);

set(h2b,'Units','Inches');
pos = get(h2b,'Position');
set(h2b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h2b,figureName_v_cr,'-dpdf','-r300');

cd ..



%% finding min h_avg_probe for a given film parameters for different film radii

% relFolder = {'h0_200nm_Avw_1.25e-21_ST_0.034_Rc_0.0018_disjPr_on',...
%     'h0_300nm_Avw_1.25e-21_ST_0.034_Rc_0.0018_disjPr_on'...
%     'h0_500nm_Avw_1.25e-21_ST_0.034_Rc_0.0018_disjPr_on',...
%     'h0_1000nm_Avw_1.25e-21_ST_0.034_Rc_0.0018_disjPr_on'};
% h2b = figure;
% h2b.Renderer = 'Painters';
% for j = 1:length(relFolder)
%     pathFol{j} = strcat('..\..\..\finiteSizeStochasticThinFilms\',relFolder{j});
%     cd(pathFol{j})
%     folders_R = dir('Rf*');
%     a2 = cellfun(@num2str, struct2cell(folders_R), 'UniformOutput', false);
%     Out2 = sortrows(a2.',6);
%     hAvg_Right_relevant = zeros(length(relFolder), max(size(Out2)));
%     for iter = 1:max(size(Out2))
%         subFolderR = Out2{iter};
%         cd(subFolderR);
%         load('hProbeAvg.mat');
%         hAvg_Right_relevant(j,iter) = h_right_avg_j(end);
%         cd ..
%     end
%     load('results_differentFilmSize.mat');
%     Legend_h_avg{j} = strcat('$h_0 = $', num2str(h0_init*10^9),' nm');
%     cd(strcat('../../postProcess/postProcessFiniteSizedFilms/',mk_postProcess))
%     loglog(R_film, hAvg_Right_relevant(j,:)*h0_init*10^9, 'o')
%     hold on
% end
% ylim([15 80])
% xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
% ylabel('$h_{avg}(t=t_r)$ (nm)','Fontsize',14)
% legend(Legend_h_avg, 'interpreter','latex','location','best')
% set(gca,'FontSize',16)
% 
% set(h2b,'Units','Inches');
% pos = get(h2b,'Position');
% set(h2b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h2b,'averageFilmHeight_atRupture','-dpdf','-r300');
% 
% save('h_avg_diffFilmSizes.mat')
toc