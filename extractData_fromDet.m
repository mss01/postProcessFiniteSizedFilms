function [t_rupt_det t_drain_det_wholeFilm t_drain_det_left t_drain_det_right v_thin_min_det h_cr_final_FullFilmavg h_cr_det v_thin_centre v_thin_rim] = extractData_fromDet(L_film_det)

str1 = strcat('../Rf_', num2str(L_film_det),'_mu_m'); %num2str(L_film(i));  %

cd(str1)
matFilesToRead = dir('*.mat'); 
for q = 1:length(matFilesToRead)
    load(matFilesToRead(q).name); 
end 
t_drain_det_wholeFilm = drainageTime;
t_drain_det_left = drainageTime_left;
t_drain_det_right = drainageTime_right;
t_rupt_det = t_rupt;

v_thin_min_det = avg_cr_thinningRate_fit;
h_cr_det = h_cr_final;

cd('Z:\from_shah_ms\simulations\firstPaperSimulations\doubleSided_Curvature\deterministic\Manev_etal_1984a_Malhotra_etal1987_reverseEngg\h0_75nm_RadoevData\probeCentredAtDimple\postProcessAll');

end