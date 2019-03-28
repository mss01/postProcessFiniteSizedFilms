clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

cd ..
data_fromExcel = xlsread('dataRadoev1984.xlsx','Radoevdata');

%% input parameters
h0 = 200e-9;
h_dimensionless = [0.05:0.05:1];
A = 1.25e-21;
Rc = 10e-3;
gamma = 0.034;
visc = 0.00089;
k = 1.38*10^-23;
Tmp = 300; 
a0 = 1;
h_drain_start = 100e-9/h0;          % start recording  
h_drain_end = 25e-9/h0;             % was used to end recording times

%% dependent dimensionless parameters

kappa = round(pi*h0^3*gamma/(A*Rc));
t_scale = 12*pi^2*visc*gamma*h0^5/A^2;
l_scale = h0^2*sqrt(2*pi*gamma/A);

%% length of the  films studied

L_film_det  = [50 60 65 70 75 80 90 100 115 150 200 300 400 500 600 700 800 900 1000]';
L_film_Radoev = [50 60 65 70 75 80 90 100 115 150 200 500 700 1000].*10^-6;
L_film_scaled = L_film_det.*10^-6./l_scale;

%% extract deterministic data

for i = 1:length(L_film_det);
        [t_rupt_det(i,:) t_drain_det_wholeFilm(i,:) t_drain_det_left(i,:) t_drain_det_right(i,:) v_thin_min_det(i,:) h_cr_final_FullFilmavg(i,:) h_cr_det(i,:) v_thin_centre(i,:) v_thin_rim(i,:)] = ...
        extractData_fromDet(L_film_det(i));
    L_film_det(i)
    close all;
end

t_rupt_det = t_rupt_det.*t_scale;
t_drain_det_wholeFilm = t_drain_det_wholeFilm.*t_scale;
t_drain_det_left = t_drain_det_left.*t_scale;
t_drain_det_right = t_drain_det_right.*t_scale;
drainageTime_det = [t_drain_det_left((1:end),:)];
v_thin_min_det = abs(v_thin_min_det).*h0*10^10./t_scale;
h_cr_det_final = [h_cr_final_FullFilmavg((1:3),:); h_cr_det((4):end,:)].*h0*10^10;
v_thin_avg_det = (h_drain_start - h_drain_end)*h0*10^10./(drainageTime_det);

save('extracted_det_Data01.mat')

% load('extracted_det_Data01.mat');

%% Reynolds theory followed by MTR theory

v_re_det = ((1 + 6.*h_dimensionless.^3*kappa)./(L_film_det.*10^-6/l_scale).^2)*h0/t_scale*10^10;
% t_re = h0*10^10./v_re_det;
preFactorFunc = @(x) 1./(1+6.*kappa.*x.^3);
preFactorRe = integral(preFactorFunc, h_drain_end ,h_drain_start);
t_re = preFactorRe.*L_film_scaled.^2.*t_scale;             %% 0.241 is coming from integrating 1/(1 + 6h^3*kappa) from 0 to 1

% v_MTR = (((1 + 6*kappa)^8.*h_dimensionless.^12)./(3.*L_film_scaled.^4)).^0.2*h0/t_scale*10^10;
for i = 1:length(h_dimensionless)
    v_MTR(:,i) = 1/6.*(((1 + 6*kappa.*h_dimensionless(i).^3).^8)./(108.*h_dimensionless(i).^12*L_film_scaled.^4)).^0.2*h0/t_scale*10^10;
end

preFactorFuncMTR = @(y) 6.*(108.*y.^12).^(1/5)./(1+6.*kappa.*y.^3).^(8/5);
% preFactorFuncMTR = @(y) 1./(1+6.*kappa.*y.^3);

preFactorMTR = integral(preFactorFuncMTR, h_drain_end ,h_drain_start);
t_MTR = preFactorMTR.*(L_film_scaled.^4.).^(1/5).*t_scale;

h_MTR = h0.*h_dimensionless;
for i = 1:length(h_dimensionless)
        v_MTR_1997Paper(:,i) = 8/(3*visc)*((4.*h_MTR(i).^12.*(gamma/Rc + A/(6*pi*h_MTR(i).^3)).^8)./(3.83^12*gamma^3.*(L_film_det.*10^-6).^4)).^0.2.*10^10;%% from Manev 1997 paper
        v_MTR_Tsekov(:,i) = 1/(6*visc).*((h_MTR(i).^12*(gamma/Rc + A/(6*pi*h_MTR(i).^3)).^8)./(4*gamma^3.*(L_film_det.*10^-6).^4)).^0.2.*10^10;
end


%% data from Radoev paper

L_film_Radoev2 = data_fromExcel(:,28)*10^-6;

L_film_Radoev_sort = [L_film_Radoev2(1) L_film_Radoev2(2) mean(L_film_Radoev2(3:4)) ...
    mean(L_film_Radoev2(5:6)) mean(L_film_Radoev2(7:9)) mean(L_film_Radoev2(10:12)) ...
    mean(L_film_Radoev2(13:15)) mean(L_film_Radoev2(16:25)) mean(L_film_Radoev2(26)) mean(L_film_Radoev2(27)) ...
    mean(L_film_Radoev2(28:34)) mean(L_film_Radoev2(35:49)) mean(L_film_Radoev2(50:53)) ...
    mean(L_film_Radoev2(54:61))];

v_thin_Radoev2 = data_fromExcel(:,30);
v_thin_Radoev_sort = [v_thin_Radoev2(1) v_thin_Radoev2(2) mean(v_thin_Radoev2(3:4)) ...
    mean(v_thin_Radoev2(5:6)) mean(v_thin_Radoev2(7:9)) mean(v_thin_Radoev2(10:12)) ...
    mean(v_thin_Radoev2(13:15)) mean(v_thin_Radoev2(16:25)) mean(v_thin_Radoev2(26)) mean(v_thin_Radoev2(27)) ...
    mean(v_thin_Radoev2(28:34)) mean(v_thin_Radoev2(35:49)) mean(v_thin_Radoev2(50:53)) ...
    mean(v_thin_Radoev2(54:61))];

std_v_thin_Radoev = [0 0 std(v_thin_Radoev2(3:4)) ...
    std(v_thin_Radoev2(5:6)) std(v_thin_Radoev2(7:9)) std(v_thin_Radoev2(10:12)) ...
    std(v_thin_Radoev2(13:15)) std(v_thin_Radoev2(16:25)) 0 0 ...
    std(v_thin_Radoev2(28:34)) std(v_thin_Radoev2(35:49)) std(v_thin_Radoev2(50:53)) ...
    std(v_thin_Radoev2(54:61))];

h_cr_Radoev2 = data_fromExcel(:,29);
sizeRadoevsData = [1 1 2 2 3 3 3 10 1 1 7 15 4 8];
h_cr_Radoev_sort = [h_cr_Radoev2(1) h_cr_Radoev2(2) mean(h_cr_Radoev2(3:4)) ...
    mean(h_cr_Radoev2(5:6)) mean(h_cr_Radoev2(7:9)) mean(h_cr_Radoev2(10:12)) ...
    mean(h_cr_Radoev2(13:15)) mean(h_cr_Radoev2(16:25)) mean(h_cr_Radoev2(26)) mean(h_cr_Radoev2(27)) ...
    mean(h_cr_Radoev2(28:34)) mean(h_cr_Radoev2(35:49)) mean(h_cr_Radoev2(50:53)) ...
    mean(h_cr_Radoev2(54:61))];

std_h_cr_Radoev = [0 0 std(h_cr_Radoev2(3:4)) ...
    std(h_cr_Radoev2(5:6)) std(h_cr_Radoev2(7:9)) std(h_cr_Radoev2(10:12)) ...
    std(h_cr_Radoev2(13:15)) std(h_cr_Radoev2(16:25)) 0 0 ...
    std(h_cr_Radoev2(28:34)) std(h_cr_Radoev2(35:49)) std(h_cr_Radoev2(50:53)) ...
    std(h_cr_Radoev2(54:61))];

%% theories on critical thicknesses

% Vrij's thoery

h_cr_vrij = 0.268*(A^2*(L_film_det.*10^-6).^2.*Rc./gamma^2).^(1/7)*10^10;                   % from Vrij's 1966 paper
h_cr_vrij_dim = 0.268*h0*(2*pi^2.*L_film_scaled.^2./kappa).^(1/7).*10^10;                   % from Vrij's 1966 paper nondimensionlized based on the scales
h_cr_MTR = 0.98.*(k*Tmp)^(1/12)*(A/(6*pi))^(1/3)./(visc^(1/6)*gamma^(1/4)*(a0)^(1/6)).*(L_film_det.*10^-6).^(2/15).*10^10;   % from 2005 Manev and Angarska paper
h_cr_Radoev = (9.*sqrt(3).*sqrt(k*Tmp/gamma)*(A/(6*pi))^2./(16.*v_thin_Radoev_sort.*10^-10.*visc*gamma)).^(1/5).*10^10;

%% Reynolds theory

% v_re_det = ((1 + 6.*h_dimensionless.^3*kappa)./(L_film_det.*10^-6/l_scale).^2)*h0/t_scale*10^10;
% t_re = h0*10^10./v_re_det;
% 
% v_re_Radoev = ((1 + 6*h_dimensionless^3*kappa)./(L_film_Radoev/l_scale).^2)*h0/t_scale*10^10;
% v_re_stoc = ((1 + 6*h_dimensionless^3*kappa)./((L_film_stoc*10^-6)/l_scale).^2)*h0/t_scale*10^10;


%% comparison of drainage times

% may be at some point we look at the rupture times from the simulations
h1a = figure;
loglog(L_film_det, t_rupt_det, 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t$ (s)','Fontsize',14)
set(gca,'FontSize',14)

set(h1a,'Units','Inches');
pos = get(h1a,'Position');
set(h1a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1a,'ruptureTime_h0_500nm','-dpdf','-r300')

% next is drainage time deterministic
h1b = figure;
loglog(L_film_det, drainageTime_det, 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t_{drain}$ (s)','Fontsize',14)
set(gca,'FontSize',14)


% this one is for getting a handle on how much Reynolds theory predicts

% loglog(L_film*10^6, t_re, 'o')

% hold on
% errorbar(L_film_stoc, mean_drainageTime, std_drainageTime,'o')
% set(gca, 'Xscale', 'log')
% set(gca, 'Yscale', 'log')

hold on 
loglog(L_film_det, t_re, 'o')
hold on
loglog(L_film_det, t_MTR,'o')
ylim([10 1200])
xlim([8 1100])
legend('$\theta$ = 0','Reynolds theory','MTR theory','Location','best')

set(h1b,'Units','Inches');
pos = get(h1b,'Position');
set(h1b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1b,'drainageTime_Radoe_h0_500nm','-dpdf','-r300')

%% Now on to the thinning rates

% for probe sensititvity check RadoevData_01.m file

h2b = figure;

loglog(L_film_det, v_thin_min_det, 'o')

% hold on
% errorbar(L_film_stoc, mean_v_thin_min_stoc_dim, std_v_thin_min_stoc_dim, 'o')

hold on
errorbar(L_film_Radoev_sort*10^6, v_thin_Radoev_sort, std_v_thin_Radoev, 'o')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')

hold on
loglog(L_film_det, v_re_det(:,7),'o')
hold on
loglog(L_film_det, v_MTR(:,7),'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$V$ ($\AA/s$)','Fontsize',14)
set(gca,'FontSize',14)
ylim([0.5 120])


legend('thinning rate, $\theta = 0$', 'Radoev thin rate', 'Reynolds thinning rate', 'MTR thinning rate')

set(h2b,'Units','Inches');
pos = get(h2b,'Position');
set(h2b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2b,'velThinningRadoev_h0_500nm','-dpdf','-r300')


h2c = figure;

plot(L_film_det, v_thin_min_det, 'o')

% hold on
% errorbar(L_film_stoc, mean_v_thin_min_stoc_dim, std_v_thin_min_stoc_dim, 'o')

hold on
errorbar(L_film_Radoev*10^6, v_thin_Radoev_sort, std_v_thin_Radoev, 'o')

hold on
plot(L_film_det, v_re_det(:,7),'o')
hold on
plot(L_film_det, v_MTR(:,7),'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$V$ ($\AA/s$)','Fontsize',14)
set(gca,'FontSize',14)

legend('thinning rate, $\theta$ = 0','Radoev thinning rate', 'Reynolds thinning rate', 'MTR thinning rate')

set(h2c,'Units','Inches');
pos = get(h2c,'Position');
set(h2c,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2c,'velThinningRadoev_linlin_h0_500nm','-dpdf','-r300')

%% CR

c_r_det = 2*h0*Rc./(L_film_det*10^-6).^2;
c_r_Radoev = 2*h0*Rc./L_film_Radoev.^2;
c_r_stoc = 2*h0*Rc./(L_film_stoc*10^-6).^2;

v_rim = v_thin_rim./v_re_det(:,end);
v_centre = (v_thin_centre./v_re_det(:,end));

h3 = figure;
loglog(c_r_det, v_rim, 'o')
xlabel('$C_R$','Fontsize',14)
ylabel('$V/V_{Re}$','Fontsize',14)
set(gca,'FontSize',14)

hold on
loglog(c_r_det, v_centre, 'o')
% ylim([10^-4 13])


legend('thinning rate rim','thinning rate centre')
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h3,'Vel_By_V_Re_C_r_h0_500nm','-dpdf','-r300')

%% critical thicknesses

% h4a = figure;
% for i = 1:3
%     loglog(L_film_det, h_cr_det(:,i), markIdhar{i}, 'color', colorIdhar{i})
%     hold on
% end
% xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
% ylabel('$h_{cr}$ ($\AA$)','Fontsize',14)
% legend('res = 20 $\mu$m','res = 22 $\mu$m', 'res = 24 $\mu$m','location','southeast')
% set(gca,'FontSize',14)
% 
% set(h4a,'Units','Inches');
% pos = get(h4a,'Position');
% set(h4a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h4a,'h_cr_Radoev_SensAnalysis','-dpdf','-r300')



h4b = figure;
loglog(L_film_det, h_cr_det_final(:,3), 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$h_{cr}$ ($\AA$)','Fontsize',14)
set(gca,'FontSize',14)

% hold on
% errorbar(L_film_stoc, mean_h_cr_stoc_dim, std_h_cr_stoc_dim, 'o')
% set(gca, 'Xscale', 'log')
% set(gca, 'Yscale', 'log')

hold on
errorbar(L_film_Radoev*10^6, h_cr_Radoev_sort, std_h_cr_Radoev, 'o')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')


loglog(L_film_det, h_cr_vrij, 'k--')
hold on
loglog(L_film_det, h_cr_MTR, 'k-.')
hold on
loglog(L_film_Radoev_sort.*10^6, h_cr_Radoev, 'k+') 
set(gca,'FontSize',14)

% ylim([140 450])
legend('$\theta$ = 0', 'Radoev expt. data', 'Vrij','MTR','Radoev','Location','best')

set(h4b,'Units','Inches');
pos = get(h4b,'Position');
set(h4b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h4b,'h_cr_Radoev_withVrij_h0_500nm','-dpdf','-r300')


% mean_h_cr_stoc_dim = [mean_h_cr_full_stoc(1:5) mean_h_cr_stoc(6:end)];
% std_h_cr_stoc_dim = [std_h_cr_full_stoc(1:5) std_h_cr_stoc(6:end)];

% h4b = figure;
% loglog(L_film_det, h_cr_det_final(:,3), 'o')
% xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
% ylabel('$h_{cr}$ ($\AA$)','Fontsize',14)
% set(gca,'FontSize',14)
% 
% hold on
% errorbar(L_film_stoc, mean_h_cr_stoc_dim, std_h_cr_stoc_dim, 'o')
% set(gca, 'Xscale', 'log')
% set(gca, 'Yscale', 'log')
% 
% hold on
% errorbar(L_film_Radoev*10^6, h_cr_Radoev_sort, std_h_cr_Radoev, 'o')
% set(gca, 'Xscale', 'log')
% set(gca, 'Yscale', 'log')
% 
% y_hcr_vrij = 6e3.*(L_film_stoc.*10^-6).^(2/7);
% y_hcr_radoev = 0.88e3*(L_film_det.*10^-6).^0.13;
% hold on
% loglog(L_film_stoc, y_hcr_vrij, 'k')
% hold on
% loglog(L_film_det, y_hcr_radoev, 'k--')
% ylim([140 450])
% legend('$\theta$ = 0','$\theta$ = 0.001', 'Radoev expt. data', 'slope 2/7 - Vrij','slope 0.13 - Radoev','Location','northwest')
% 
% set(h4b,'Units','Inches');
% pos = get(h4b,'Position');
% set(h4b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h4b,'h_cr_Radoev01','-dpdf','-r300')
%% 

stdErMean = std_v_thin_Radoev./sizeRadoevsData;

