PRO cross_difference, choose_term = choose_term

;; get constants from the Aug23 run:
hubble_param = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'HUBBLE_PARAM', names = varnames)
kperp_lambda_conv = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'KPERP_LAMBDA_CONV')
delay_params = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'DELAY_PARAMS')
kx_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'KX_MPC')
ky_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'KY_MPC')
kz_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
     'KZ_MPC')
     
;max_kperp_lambda = min([file_struct_arr.kspan/2.,file_struct_arr.max_baseline_lambda])    
kperp_plot_range = [5./kperp_lambda_conv, 300./kperp_lambda_conv]

;; get data
IF SIZE(choose_term, /N_ELEMENTS) LT 1 THEN choose_term = 1
CASE choose_term OF
1: BEGIN
     ;; get data and variances for term 1:
     Aug23_diff = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
          'DATA_DIFF_1')
     Aug23_sigma2 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $ 
          'SIGMA2_1')
     Aug27_diff = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'DATA_DIFF_1')
     Aug27_sigma2 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'SIGMA2_1') 
END
2: BEGIN
     ;; get data and variances for term 2
     Aug23_diff = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'DATA_DIFF_2')
     Aug23_sigma2 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'SIGMA2_2')
     Aug27_diff = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'DATA_DIFF_2')
     Aug27_sigma2 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
          'SIGMA2_2')
END
ENDCASE

;; calculate cross difference and weights for term 1        
diff_cross = Aug23_diff * CONJ(Aug27_diff)
sigma2 = diff_cross^2 * (Aug23_sigma2/Aug23_diff^2 + Aug27_sigma2/Aug27_diff^2)
weights = 1/sigma2
wh_diff0 = WHERE(Aug23_diff EQ 0 OR Aug27_diff EQ 0, count_diff0)
IF count_diff0 GT 0 THEN weights[wh_diff0] = 0
wh_sig0 = WHERE(sigma2 EQ 0, count_sig0)
IF count_sig0 GT 0 THEN weights[wh_sig0] = 0

;; rebin to 2D
diff_cross_2d = kspace_rebinning_2d(diff_cross, kx_mpc, ky_mpc, kz_mpc, $
     kperp_edges_mpc, kpar_edges_mpc, weights = weights, kperp_bin = kperp_bin, kpar_bin = kpar_bin)
     
;; plot:
kpower_2d_plots, power = real_part(diff_cross_2d), kperp_edges = kperp_edges_mpc, kpar_edges = kpar_edges_mpc, $
     kperp_bin = kperp_bin, kpar_bin = kpar_bin, kperp_lambda_conv = kperp_lambda_conv, $
     delay_params = delay_params, hubble_param = hubble_param, window_num = 1, color_profile = 'sym_log', $
     /delay_axis, /hinv, /baseline_axis, kperp_plot_range = kperp_plot_range, full_title = 'Difference Cross - Real Part', $
     note = 'Long Run Aug23 Difference * Aug27 Difference'
     
kpower_2d_plots, power = abs(diff_cross_2d), kperp_edges = kperp_edges_mpc, kpar_edges = kpar_edges_mpc, $
     kperp_bin = kperp_bin, kpar_bin = kpar_bin, kperp_lambda_conv = kperp_lambda_conv, $
     delay_params = delay_params, hubble_param = hubble_param, window_num=2, /delay_axis, /hinv, /baseline_axis, $
     kperp_plot_range = kperp_plot_range, full_title = 'Difference Cross Term' + STRCOMPRESS(choose_term) + '- Abs Value', $
     note = 'Long Run Aug23 Difference * Aug27 Difference'





END