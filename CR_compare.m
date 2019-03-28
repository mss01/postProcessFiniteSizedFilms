clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

%% h0 = 500nm,  kappa = 6000 (no disj press) Rf = 25;
load('..\h0_500nm_Avw_1.25e-21_ST_0.034_Rc_0.001\Rf_25_mu_m\JoyeComparison.mat');
close all

c_r_Joye_r_comb{1} = c_r_Joye_r(2:end);
ratio_v_vre_comb{1} = ratio_v_vre;
ratio_vc_vre_comb{1} = ratio_vc_vre;
c_r_Joye_centre_comb{1} = c_r_Joye_centre(2:end);

%% h0 = 500nm, kappa = 6000, Rf = 50

load('..\h0_500nm_Avw_1.25e-21_ST_0.034_Rc_0.001\Rf_50_mu_m\JoyeComparison.mat');
close all

c_r_Joye_r_comb{2} = c_r_Joye_r(2:end);
ratio_v_vre_comb{2} = ratio_v_vre;
ratio_vc_vre_comb{2} = ratio_vc_vre;
c_r_Joye_centre_comb{2} = c_r_Joye_centre(2:end);


%% h0 = 500nm, kappa = 6000 (no disj press) Rf = 100

load('..\h0_500nm_Avw_1.25e-21_ST_0.034_Rc_0.001\Rf_100_mu_m\JoyeComparison.mat');
close all

c_r_Joye_r_comb{3} = c_r_Joye_r(2:end);
ratio_v_vre_comb{3} = ratio_v_vre;
ratio_vc_vre_comb{3} = ratio_vc_vre;
c_r_Joye_centre_comb{3} = c_r_Joye_centre(2:end);

hfig3 = figure;

loglog(c_r_Joye_r_comb{1}, ratio_v_vre_comb{1}, 'o')
hold on
loglog(c_r_Joye_centre_comb{1}, ratio_vc_vre_comb{1},'--')
hold on
loglog(c_r_Joye_r_comb{2}, ratio_v_vre_comb{2}, 'o')
hold on
loglog(c_r_Joye_centre_comb{2}, ratio_vc_vre_comb{2},'--')
hold on
loglog(c_r_Joye_r_comb{3}, ratio_v_vre_comb{3}, 'o')
hold on
loglog(c_r_Joye_centre_comb{3}, ratio_vc_vre_comb{3},'--')

xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

set(hfig3,'Units','Inches');
pos = get(hfig3,'Position');
set(hfig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig3,'compareJoyeFig8','-dpdf','-r300')


%% h0 = 500nm, kappa = 6000, Rf = 10

hfig4 = figure;
load('Z:\from_shah_ms\simulations\firstPaperSimulations\doubleSided_Curvature\deterministic\Manev_etal_1984a_Malhotra_etal1987_reverseEngg\h0_500nm_Radoev\deterministic\Rf_10_mu_m\workspace_deterministic_t_cr.mat');

loglog(c_r_Joye(2:end), ratio_v_vre.*8/3, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre.*8/3,'--')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

%% h0 = 500nm, kappa = 6000, Rf = 30

load('Z:\from_shah_ms\simulations\firstPaperSimulations\doubleSided_Curvature\deterministic\Manev_etal_1984a_Malhotra_etal1987_reverseEngg\h0_500nm_Radoev\deterministic\Rf_30_mu_m\workspace_deterministic_t_cr.mat');

hold on
loglog(c_r_Joye(2:end), ratio_v_vre.*8/3, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre.*8/3,'--')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

%% h0 = 500nm, kappa = 6000, Rf = 50

load('Z:\from_shah_ms\simulations\firstPaperSimulations\doubleSided_Curvature\deterministic\Manev_etal_1984a_Malhotra_etal1987_reverseEngg\h0_500nm_Radoev\deterministic\Rf_50_mu_m\workspace_deterministic_t_cr.mat');

hold on
loglog(c_r_Joye(2:end), ratio_v_vre.*8/3, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre.*8/3,'--')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

%% h0 = 500nm, kappa = 6000, Rf = 100

load('Z:\from_shah_ms\simulations\firstPaperSimulations\doubleSided_Curvature\deterministic\Manev_etal_1984a_Malhotra_etal1987_reverseEngg\h0_500nm_Radoev\deterministic\Rf_100_mu_m\workspace_deterministic_t_cr.mat');

hold on
loglog(c_r_Joye(2:end), ratio_v_vre.*8/3, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre.*8/3,'--')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)


set(hfig4,'Units','Inches');
pos = get(hfig4,'Position');
set(hfig4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig4,'thinningRates_vs_C_R_Compare_withDisjPress_cart','-dpdf','-r300')