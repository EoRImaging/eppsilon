function translate_1d_power,k_edges,k_bin,hubble_param,power,k_out=k_out
  kscale = k_edges[1:*]/hubble_param
  k_mid = kscale[1:*] - k_bin/hubble_param
  ;k_out = [min(kscale), k_mid, max(kscale)]
  k_out = [min(kscale)-k_bin/hubble_param/2, k_mid, max(kscale),max(kscale)+k_bin/hubble_param]
  power_out = hubble_param^3. * [power[1],power[1:-1],power[-1]]  
  
  return,power_out
end

pro make_plots_for_adams_sem1_paper,diffuse_sky=diffuse_sky,reflection_diff=reflection_diff,$
  diffuse_diff=diffuse_diff,kpar_plot=kpar_plot,kperp_plot=kperp_plot,final_cut_2d=final_cut_2d,$
  plot_all=plot_all,include_comparison=include_comparison
  ; This is a script to generate plots for Adam's paper (and record how they were made).
  ; This is to be run _after_ one round of PS has already been run using the proper
  ; cluster management.
  ; Note that I moved the ps data into ps_jan2016, thus the ps_foldername option
  
  if n_elements(include_comparison) eq 0 then include_comparison=1 ; do it by default
    
  ; Directory to put all the plots;
  outdir = '/nfs/eor-00/h1/beards/temp/plots_for_paper/'
  if not file_test(outdir,/directory) then file_mkdir,outdir
  


  if keyword_set(diffuse_sky) or keyword_set(plot_all) then begin
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
    outfile=outdir+'diffuse_sky'
    quick_image,im,xtitle='ra (degrees)', ytitle='dec (degrees)',xrange=ra_range,yrange=dec_range,savefile=outfile,/pdf,data_range=[-1e4,1e4]
  endif

  ; ###################
  ; Intermediate PS
  ; ###################
  
  if keyword_set(reflection_diff) or keyword_set(plot_all) then begin
    ; Reflection diff
    outfile=outdir+'reflections_diff.pdf'
    obs_name=['golden','golden']
    dir=['/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_before_reflections/','/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_after_reflections/']
    ps_diff_wrapper,dir,obs_name,/exact_obsnames,/pdf
    file_copy,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_reflections__before_minus_after/plots/fhd_apb_reflections_bh__before_golden_minus_after_golden_dencorr_2dkpower.pdf',outfile,/overwrite
  endif
  
  if keyword_set(diffuse_diff) then begin
    ; Diffuse diff
    outfile=outdir+'diffuse_diff.pdf'
    obs_name=['golden','golden']
    dir=['/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_jan2016_no_diffuse/','/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_after_diffuse/']
    ps_diff_wrapper,dir,obs_name,/exact_obsnames,/pdf
    file_copy,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb__jan2016_no_diffuse_minus_after_diffuse/plots/fhd_apb_bh__jan2016_no_diffuse_golden_minus_after_diffuse_golden_dencorr_2dkpower.pdf',outfile,/overwrite
  endif
  
 
  ; ###################
  ; Deep integration
  ; ###################
  
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1'
  obs_name='wedge_cut_plus_res_cut'
  low=[0,95]
  mid=[96-48,191-48]
  hi=[96,191]
  
  ;kperp_range_lambda_use=[10,70]
  kperp_range_lambda_use=[10,60]
  kpar_range_1dave_use=[0.15,200]
  
  ; Few plotting options that should be uniform
  linethick=5
  charsize=2.5
  
  kperp_slice=[13,16]
  kpar_slice=[1.05,1.15]

  if keyword_set(kpar_plot) or keyword_set(plot_all) then begin
    ; Make a narrow kpar power plot
        
    ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
      freq_ch_range=mid,kperp_range_lambda_1dave=kperp_slice,$
      kpar_range_1dave=[0,200],ps_foldername='ps_jan2016',$
      /plot_kpar_power,type_inc=['dirty','res'],pol_inc='xx'
    restore,dir+'/ps_jan2016/Combined_obs_wedge_cut_plus_res_cut_cubeXX__even_odd_joint_ch48-143_res_xx_bh_dencorr_kperplambda'+$
          strtrim(kperp_slice[0],2)+'-'+strtrim(kperp_slice[1],2)+'_kpar_power.idlsave'
    coarse_harm_ranges=transpose([[.244,.372],[.552,.680],[.860,.988],[1.168,1.296]])/hubble_param ; These are really hard to get... had to dig into plotting code. These are specifically for mid-band, coarse_harm_width=3
    res_power_out=translate_1d_power(k_edges,k_bin,hubble_param,power,k_out=k_out)
    noise_out=translate_1d_power(k_edges,k_bin,hubble_param,1./sqrt(weights),k_out=k_out)
    lin_delay_kpar_slope = (delay_params[1] - delay_params[0])/(max(k_edges)/hubble_param - k_bin)
    lin_delay_kpar_intercept = delay_params[0] / (lin_delay_kpar_slope * k_bin)
    ;delay_edges = lin_delay_kpar_slope * k_edges + lin_delay_kpar_intercept
    ; plot it
    outfile=outdir+'kpar.ps'
    cgPS_Open,File=outfile
    device,xsize=10,ysize=4,/inches
    xrange=[.01,max(coarse_harm_ranges)]
    yrange=[1e5,1e13]
    axis_range = xrange*lin_delay_kpar_slope + lin_delay_kpar_intercept
    axis_title = 'delay (ns)'
    xtitle=textoidl('k_{||} (h Mpc^{-1})', font = 1)
    ytitle=textoidl('P_k (mK^2 h^{-3} Mpc^3)', font = 1)
    
    cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
      psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=linethick,color='blue',$
      xtitle=xtitle,ytitle=ytitle,charsize=charsize,linestyle=2
    cgcolorfill,[xrange[0],kpar_range_1dave_use[0],kpar_range_1dave_use[0],xrange[0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
    ;cgcolorfill,[kpar_range_1dave[1],xrange[1],xrange[1],kpar_range_1dave[1]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
    for i=0,3 do begin
      cgcolorfill,[coarse_harm_ranges[i,0],coarse_harm_ranges[i,1],coarse_harm_ranges[i,1],coarse_harm_ranges[i,0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
    endfor
    ; replot to fix damage done by polyfill
    ;cgcolorfill,[kpar_slice[0],kpar_slice[1],kpar_slice[1],kpar_slice[0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("plum")
    ;cgoplot,[mean(kpar_slice),mean(kpar_slice)],[yrange[0],yrange[1]],col='blue',thick=3*linethick
    cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
      psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=linethick,color='blue',$
      xtitle=xtitle,ytitle=ytitle,charsize=charsize,/noerase,linestyle=2
    ; draw top axis
    cgaxis, xaxis=1, xrange = axis_range,  xtitle = axis_title, font = 1, xstyle = 1,charsize=charsize
    cgplot,/overplot,k_out,noise_out,psym=10,linestyle=2,thick=linethick,color='red'
    
    cgPS_Close,/PDF,/Delete_PS
    
    ;;;; Comparison version
    if keyword_set(include_comparison) then begin
      ;xrange=[.01,2*max(coarse_harm_ranges)]
      ; Read in Caths data
      textfast,chips_kpar,file_path=dir+'/chips_data/kpar',/read
      ; NEED TO CONFIRM WITH CATH - I think these are edges, not centers
      dkpar = chips_kpar[1]-chips_kpar[0]
      chips_kpar = chips_kpar + dkpar/2
      textfast,chips_ppar,file_path=dir+'/chips_data/power_vs_kpar',/read
      outfile=outdir+'kpar_comparison.ps'
      cgPS_Open,File=outfile
      device,xsize=10,ysize=4,/inches
      cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
        psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=2*linethick,color='blue',$
        xtitle=xtitle,ytitle=ytitle,charsize=charsize,linestyle=2
      cgcolorfill,[xrange[0],kpar_range_1dave_use[0],kpar_range_1dave_use[0],xrange[0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
      for i=0,3 do begin
        cgcolorfill,[coarse_harm_ranges[i,0],coarse_harm_ranges[i,1],coarse_harm_ranges[i,1],coarse_harm_ranges[i,0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
      endfor
      ; replot to fix damage done by polyfill
      ;cgcolorfill,[kpar_slice[0],kpar_slice[1],kpar_slice[1],kpar_slice[0]],[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("plum")
      ;cgoplot,[mean(kpar_slice),mean(kpar_slice)],[yrange[0],yrange[1]],col='blue',thick=3*linethick
      cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
        psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=2*linethick,color='blue',$
        xtitle=xtitle,ytitle=ytitle,charsize=charsize,/noerase,linestyle=2
      cgplot,/overplot,chips_kpar,chips_ppar,psym=10,thick=linethick-2,color=cgcolor('dark gray'),$
        linestyle=0
      ; draw top axis
      cgaxis, xaxis=1, xrange = axis_range,  xtitle = axis_title, font = 1, xstyle = 1,charsize=charsize
      cgplot,/overplot,k_out,noise_out,psym=10,linestyle=0,thick=linethick,color='red'
        
      cgPS_Close,/PDF,/Delete_PS
    endif
    
    
  endif
  
  if keyword_set(kperp_plot) or keyword_set(plot_all) then begin
    ; Narrow kperp power plot
    ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
      freq_ch_range=mid,kperp_range_lambda_1dave=[0,2000],$
      kpar_range_1dave=kpar_slice,ps_foldername='ps_jan2016',$
      /plot_kperp_power,type_inc=['dirty','res'],pol_inc='xx'
    restore,dir+'/ps_jan2016/Combined_obs_wedge_cut_plus_res_cut_cubeXX__even_odd_joint_ch48-143_res_xx_bh_dencorr_kpar'+$
          string(kpar_slice[0],format='(f4.2)')+'-'+string(kpar_slice[1],format='(f4.2)')+'_kperp_power.idlsave'
    deep_file = dir+'/ps_jan2016/Combined_obs_wedge_cut_plus_res_cut_cubeXX__even_odd_joint_ch48-143_res_xx_bh_dencorr_no_120deg_wedge_cbw3_'+$
          'kperplambda'+strtrim(kperp_range_lambda_use[0],2)+'-'+strtrim(kperp_range_lambda_use[1],2)+$
          '_kpar'+string(kpar_range_1dave_use[0],format='(f4.2)')+'-'+strtrim(fix(kpar_range_1dave_use[1]),2)+'_1dkpower.idlsave'
    wedge_amp=getvar_savefile(deep_file,'wedge_amp')
    wedge_line=mean(kpar_slice)/wedge_amp
    res_power_out=translate_1d_power(k_edges,k_bin,hubble_param,power,k_out=k_out)
    noise_out=translate_1d_power(k_edges,k_bin,hubble_param,1./sqrt(weights),k_out=k_out)
    ; Cut out power beyond kperp=0.3 to match 2D plots
    ind = where(k_out lt 0.3)
    k_out = k_out[ind]
    res_power_out = [res_power_out[ind[0:-1]],res_power_out[ind[-1]]]
    noise_out = [noise_out[ind[0:-1]],noise_out[ind[-1]]]
    ; plot it
    outfile=outdir+'kperp.ps'
    cgPS_Open,File=outfile
    device,xsize=10,ysize=4,/inches
    xrange=[.002,.5] 
    yrange=[1e6,1e10]
    axis_range = minmax(xrange * hubble_param * kperp_lambda_conv)
    axis_title = 'baseline length ' + textoidl('(\lambda)', font = 1)
    perp_char = '!9' + string(94B) + '!X'
    xtitle=textoidl('k_{perp} (h Mpc^{-1})', font = 1)
    xtitle = repstr(xtitle, 'perp', perp_char)
    ytitle=textoidl('P_k (mK^2 h^{-3} Mpc^3)', font = 1)
    cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
      psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=linethick,color='blue',$
      xtitle=xtitle,ytitle=ytitle,charsize=charsize
    ; Plot exclusion regions
    cgcolorfill,[xrange[0],kperp_range_lambda_use[0]/(hubble_param*kperp_lambda_conv),kperp_range_lambda_use[0]/(hubble_param*kperp_lambda_conv),xrange[0]],$
          [yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
    cgcolorfill,[xrange[1],kperp_range_lambda_use[1]/(hubble_param*kperp_lambda_conv),kperp_range_lambda_use[1]/(hubble_param*kperp_lambda_conv),xrange[1]],$
          [yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")    
    ;cgcolorfill,[kperp_slice[0],kperp_slice[1],kperp_slice[1],kperp_slice[0]]/(hubble_param*kperp_lambda_conv),[yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("plum")
    ;cgoplot,[mean(kperp_slice),mean(kperp_slice)]/(hubble_param*kperp_lambda_conv),[yrange[0],yrange[1]],col='blue',thick=3*linethick,linestyle=2
    ; Fix damage
    cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
      psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=linethick,color='blue',$
      xtitle=xtitle,ytitle=ytitle,charsize=charsize,/noerase      
    ; draw top axis
    cgaxis, xaxis=1, xrange = axis_range,  xtitle = axis_title, font = 1, xstyle = 1,charsize=charsize
    cgplot,/overplot,k_out,noise_out,psym=10,linestyle=0,thick=linethick,color='red'
    cgoplot,[wedge_line,wedge_line],yrange,linestyle=2,thick=linethick,color='black'
    cgPS_Close,/PDF,/Delete_PS
    
    ;;;; Comparison version
    if keyword_set(include_comparison) then begin
      ; Read in Caths data
      textfast,chips_kperp,file_path=dir+'/chips_data/kperp',/read
      ; NEED TO CONFIRM WITH CATH - I think these are edges, not centers
      dkperp = chips_kperp[1]-chips_kperp[0]
      chips_kperp = chips_kperp + dkperp/2
      textfast,chips_pperp,file_path=dir+'/chips_data/power_vs_kperp',/read
      ind=where(chips_pperp eq 0)
      remove,ind,chips_kperp
      remove,ind,chips_pperp
      outfile=outdir+'kperp_comparison.ps'
      cgPS_Open,File=outfile
      device,xsize=10,ysize=4,/inches
      cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
        psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=2*linethick,color='blue',$
        xtitle=xtitle,ytitle=ytitle,charsize=charsize,linestyle=0
      cgcolorfill,[xrange[0],kperp_range_lambda_use[0]/(hubble_param*kperp_lambda_conv),kperp_range_lambda_use[0]/(hubble_param*kperp_lambda_conv),xrange[0]],$
          [yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")
      cgcolorfill,[xrange[1],kperp_range_lambda_use[1]/(hubble_param*kperp_lambda_conv),kperp_range_lambda_use[1]/(hubble_param*kperp_lambda_conv),xrange[1]],$
          [yrange[0],yrange[0],yrange[1],yrange[1]],col=cgcolor("medium gray")    
      ; Fix damage
      cgplot,k_out,res_power_out,/ylog,/xlog,font=1,xtickformat='exponent',ytickformat='exponent',$
        psym=10,xrange=xrange,xstyle=9,yrange=yrange,/ystyle,thick=2*linethick,color='blue',$
        xtitle=xtitle,ytitle=ytitle,charsize=charsize,/noerase      
      ; draw top axis
      cgaxis, xaxis=1, xrange = axis_range,  xtitle = axis_title, font = 1, xstyle = 1,charsize=charsize
      cgplot,/overplot,k_out,noise_out,psym=10,linestyle=0,thick=linethick,color='red'
      cgoplot,[wedge_line,wedge_line],yrange,linestyle=2,thick=linethick,color='black'
      
      cgplot,/overplot,chips_kperp,chips_pperp,psym=10,thick=linethick-2,color=cgcolor('dark gray'),$
        linestyle=0
      
      
      cgPS_Close,/PDF,/Delete_PS
    endif
  endif


  if keyword_set(final_cut_2d) or keyword_set(plot_all) then begin
    ; Run the final cut
    ; Low band

    coarse_harm_width=3
    wedge_angles=120d
    ps_foldername='ps_jan2016'
    bands=[[low],[mid],[hi]]
    pols=['xx','yy']
    
    for bandi=0,2 do begin
      for poli=0,1 do begin
        ps_wrapper,dir,obs_name,/pdf,/exact_obsnames,$
          freq_ch_range=bands[*,bandi],kperp_range_lambda_1dave=kperp_range_lambda_use,$
          kpar_range_1dave=kpar_range_1dave_use,ps_foldername=ps_foldername,$
          /plot_1to2d,coarse_harm_width=coarse_harm_width,/refresh_ps,$
          wedge_angles=wedge_angles,type_inc='res',pol_inc=pols[poli]
        file=dir+'/'+ps_foldername+'/plots/fhd_apb_EoR0_high_sem1_1_'+obs_name+$
          '_ch'+strtrim(bands[0,bandi],2)+'-'+strtrim(bands[1,bandi],2)+$
          '_bh_kperplambda'+strtrim(kperp_range_lambda_use[0],2)+'-'+strtrim(kperp_range_lambda_use[1],2)+$
          '_kpar'+string(kpar_range_1dave_use[0],format='(f4.2)')+'-'+strtrim(fix(kpar_range_1dave_use[1]),2)+$
          '_res_'+pols[poli]+'_1dcontours.pdf'
        outfile=outdir+'deep_2d_ch'+strtrim(bands[0,bandi],2)+'-'+strtrim(bands[1,bandi],2)+pols[poli]+'.pdf'
        file_copy,file,outfile,/overwrite
        if (bandi eq 1) and (poli eq 0) then begin
          file=dir+'/'+ps_foldername+'/plots/fhd_apb_EoR0_high_sem1_1_'+obs_name+$
          '_ch'+strtrim(bands[0,bandi],2)+'-'+strtrim(bands[1,bandi],2)+$
          '_bh_dencorr_2dkpower.pdf'
          file_copy,file,outdir+'deep_2d_no_contours.pdf',/overwrite
        endif
      endfor
    endfor

      
  endif
  
end

