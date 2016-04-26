pro make_plots_for_adams_sem1_paper
  ; This is a script to generate plots for Adam's paper (and record how they were made).
  ; This is to be run _after_ one round of PS has already been run using the proper
  ; cluster management.
  ; Note that I moved the ps data into ps_jan2016, thus the ps_foldername option


  ; Diffuse model
  diffuse_file=filepath('EoR0_diffuse_model_94.sav',root=rootdir('FHD'),subdir='catalog_data')
  restore,diffuse_file

  map = map * float(nside2npix(nside)) / (4.*!pi)
  healpix_quickimage,map,hpx_inds,nside,ra_range=ra_range,dec_range=dec_range,map_out=map_out,/noplot ; get the map_out and ra/dec info
  
  ra_axis=lindgen((size(map_out))[1])*(ra_range[1]-ra_range[0])/(size(map_out))[1]+ra_range[0]
  dec_axis=lindgen((size(map_out))[2])*(dec_range[1]-dec_range[0])/(size(map_out))[2]+dec_range[0]
  dra = ra_axis[n_elements(ra_axis)/2]-ra_axis[n_elements(ra_axis)/2+1]
  ddec = dec_axis[n_elements(dec_axis)/2]-dec_axis[n_elements(dec_axis)/2-1]
  
  fwhm = 1. ; degrees
  fwhm = [fwhm/dra,fwhm/ddec]
  im=gauss_smooth(map_out,fwhm/2.355)
  ra_range = [-13,13]
  dec_range = [-37,-16]
  ra_ind=where((ra_axis le ra_range[1]) and (ra_axis ge ra_range[0]))
  dec_ind=where((dec_axis le dec_range[1]) and (dec_axis ge dec_range[0]))
  im=im[ra_ind,*]
  im=im[*,dec_ind]
  quick_image,im,xtitle='ra (degrees)', ytitle='dec (degrees)',xrange=ra_range,yrange=dec_range,savefile='/nfs/eor-00/h1/beards/temp/diffuse',/pdf,data_range=[-1e4,1e4]
  

  ; ###################
  ; Intermediate PS
  ; ###################
  




  ; ###################
  ; Deep integration
  ; ###################
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1'
  obs_name='wedge_cut_plus_res_cut'
  low=[0,95]
  mid=[96-48,191-48]
  hi=[96,191]

  ; Make a narrow kpar power plot
  ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
    freq_ch_range=low,kperp_range_lambda_1dave=[15,20],$
    kpar_range_1dave=[0,200],ps_foldername='ps_jan2016',$
    /plot_1to2d,/plot_kpar_power,type_inc='res',pol_inc='xx'
  ; Narrow kperp power plot
  ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
    freq_ch_range=low,kperp_range_lambda_1dave=[0,200],$
    kpar_range_1dave=[0.6,.7],ps_foldername='ps_jan2016',$
    /plot_1to2d,/plot_kperp_power,type_inc=['dirty','res'],pol_inc='xx',$
    range_1d=[1.e5,1.e11]

  ; Run the final cut
  ; Low band
  ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
    freq_ch_range=low,kperp_range_lambda_1dave=[10,70],$
    kpar_range_1dave=[0.15,200],ps_foldername='ps_jan2016',$
    /plot_1to2d,coarse_harm_width=3,/refresh_ps,wedge_angles=[120d]
  ; Mid band
  ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
    freq_ch_range=mid,kperp_range_lambda_1dave=[10,70],$
    kpar_range_1dave=[0.15,200],ps_foldername='ps_jan2016',$
    /plot_1to2d,coarse_harm_width=3,/refresh_ps,wedge_angles=[120d]
  ; High band
  ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
    freq_ch_range=hi,kperp_range_lambda_1dave=[10,70],$
    kpar_range_1dave=[0.15,200],ps_foldername='ps_jan2016',$
    /plot_1to2d,coarse_harm_width=3,/refresh_ps,wedge_angles=[120d]

end