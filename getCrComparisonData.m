function [c_r_Joye_r_comb ratio_v_vre_comb ratio_vc_vre_comb c_r_Joye_centre_comb] = getCrComparisonData(relFolder, filmSizes, JoyeStartCrit, JoyeStopCrit, mk)
% 
% str = string(filmSizes);
relSubFolder = strcat('\Rf_',filmSizes,'_mu_m');
pathFol = strcat('..\..\..\finiteSizeStochasticThinFilms\',relFolder, relSubFolder);
pathFol
cd(pathFol)
load('JoyeComparison.mat')
close all
% c_r_Joye_r_comb = c_r_Joye_r(2:end);
% ratio_v_vre_comb = ratio_v_vre;
% ratio_vc_vre_comb = ratio_vc_vre;
% c_r_Joye_centre_comb = c_r_Joye_centre(2:end);
c_r_Joye_r_comb = c_r_Joye_r(JoyeStopCrit < h_min_rim & h_min_rim < JoyeStartCrit);
ratio_v_vre = [ratio_v_vre(1) ratio_v_vre];
ratio_v_vre_comb = ratio_v_vre(JoyeStopCrit < h_min_rim & h_min_rim < JoyeStartCrit);
ratio_vc_vre = [ratio_vc_vre(1) ratio_vc_vre];
ratio_vc_vre_comb = ratio_vc_vre(JoyeStopCrit < h_min_rim & h_min_rim < JoyeStartCrit);
c_r_Joye_centre_comb = c_r_Joye_centre(JoyeStopCrit < h_min_rim & h_min_rim < JoyeStartCrit);
cd(strcat('../../../postProcess/postProcessFiniteSizedFilms/',mk));



% 
% matFilesToRead = dir('*.mat'); 
% for q = 1:length(matFilesToRead)
%     load(matFilesToRead(q).name); 
% end 

end