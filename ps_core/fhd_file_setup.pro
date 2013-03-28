function fhd_file_setup, datafile, pol, type, weightfile = weightfile, variancefile = variancefile, pixelfile = pixelfile, $
                         datavar = datavar, weightvar = weightvar, variancevar = variancevar, pixelvar = pixelvar, $
                         save_path = save_path, $
                         hpx_dftsetup_savefile = hpx_dftsetup_savefile, weight_savefilebase = weight_savefilebase_in, $
                         variance_savefilebase = variance_savefilebase_in, uvf_savefilebase = uvf_savefilebase_in, $
                         savefilebase = savefilebase_in


  nfiles = n_elements(datafile)
  if nfiles gt 2 then message, 'only 1 or 2 datafiles is supported'
 
  datafile_test = file_test(datafile)
  if min(datafile_test) eq 0 then message, 'datafile not found'
  if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
     print, 'datafiles are identical'
     datafile = datafile[0]
  endif

  if n_elements(pol) ne 1 then message, 'Polarization must be specified'
  pol_enum = ['xx', 'yy']
  wh = where(pol_enum eq pol, count)
  if count eq 0 then message, 'pol ' + pol + ' not recognized.'
  pol_num = wh[0]
 
  if n_elements(type) ne 1 then message, 'Type must be specified'
  type_enum = ['dirty', 'model', 'res']
  wh = where(type_enum eq type, count)
  if count eq 0 then message, 'type ' + type + ' not recognized.'
  type_num = wh[0]
 
  type_pol_str = type + '_' + pol
  if n_elements(datavar) eq 0 then data_varname = strupcase(type_pol_str + '_cube') else data_varname = datavar 
  if n_elements(data_varname) ne nfiles then $
     if n_elements(data_varname) eq 1 then data_varname = replicate(data_varname, nfiles) $
     else message, 'datavar must be a scalar or have the same number of elements as datafile'
  if n_elements(weightvar) eq 0 then weight_varname = strupcase('weights_' + pol + '_cube') else weight_varname = weightvar
  if n_elements(weight_varname) ne nfiles then $
     if n_elements(weight_varname) eq 1 then weight_varname = replicate(weight_varname, nfiles) $
     else message, 'weightvar must be a scalar or have the same number of elements as datafile'
  if n_elements(variancevar) eq 0 then variance_varname = strupcase('variance_' + pol + '_cube') else variance_varname = variancevar
  if n_elements(variance_varname) ne nfiles then $
     if n_elements(variance_varname) eq 1 then variance_varname = replicate(variance_varname, nfiles) $
     else message, 'variancevar must be a scalar or have the same number of elements as datafile'
  wt_file_label = '_weights_' + strlowcase(pol)
  file_label = '_' + strlowcase(type_pol_str)

  if n_elements(weightfile) eq 0 then weightfile = datafile $
  else if n_elements(weightfile) ne nfiles then message, 'weightfile must have the same number of elements as datafile'
  if n_elements(variancefile) eq 0 then variancefile = datafile $
  else if n_elements(variancefile) ne nfiles then message, 'variancefile must have the same number of elements as datafile'

  if n_elements(pixelfile) gt 0 and n_elements(pixelfile) ne nfiles then $
     message, 'pixelfile must have the same number of elements as datafile'
  if n_elements(pixelvar) gt 0 and n_elements(pixelvar) ne nfiles then $
     if n_elements(pixelvar) eq 1 then pixelvar = replicate(pixelvar, nfiles) $
     else message, 'pixelvar must be a scalar or have the same number of elements as datafile'


  if n_elements(hpx_dftsetup_savefile) gt 1 then message, 'only one hpx_dftsetup_savefile allowed'
  if n_elements(savefilebase) gt 1 then message, 'only one savefilebase allowed'

  if n_elements(weight_savefilebase) gt 0 and n_elements(weight_savefilebase) ne nfiles then $
     message, 'if weight_savefilebase is specified it must have the same number of elements as data files'
  if n_elements(uvf_savefilebase) gt 0 and n_elements(uvf_savefilebase) ne nfiles then $
     message, 'if uvf_savefilebase is specified it must have the same number of elements as data files'

  if n_elements(savefilebase_in) eq 0 or n_elements(uvf_savefilebase_in) lt nfiles then begin
     if nfiles eq 1 then begin
        if n_elements(save_path) ne 0 then froot = save_path $
        else froot = file_dirname(datafile, /mark_directory)
        uvf_froot = froot
        infilebase = file_basename(datafile)
        temp2 = strpos(infilebase, '.', /reverse_search)
        general_filebase = strmid(infilebase, 0, temp2)
        if n_elements(savefilebase_in) eq 0 then savefilebase = general_filebase + file_label $
        else savefilebase = savefilebase_in
        
        ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files 
        if n_elements(uvf_savefilebase_in) eq 0 then uvf_savefilebase = general_filebase + file_label
     endif else begin
        if n_elements(save_path) ne 0 then froot = save_path $
        else froot = file_dirname(datafile[0], /mark_directory)
        infilebase = file_basename(datafile)
        temp2 = strpos(infilebase, '.', /reverse_search)
        
        if n_elements(savefilebase_in) eq 0 then begin
           fileparts_1 = strsplit(strmid(infilebase[0], 0, temp2[0]), '_', /extract)
           fileparts_2 = strsplit(strmid(infilebase[1], 0, temp2[1]), '_', /extract)
           match_test = strcmp(fileparts_1, fileparts_2)
           wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
           if count_diff eq 0 then general_filebase = strmid(infilebase[0], 0, temp2[0]) + '_joint' $
           else begin
              if count_same gt 0 then general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
                                                     + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
              else general_filebase = infilebase[0] + infilebase[1] + '_joint'
           endelse
           savefilebase = general_filebase + file_label
        endif
        
        if n_elements(uvf_savefilebase_in) eq 0 then begin
           ;; need 2 uvf files
           if n_elements(save_path) ne 0 then uvf_froot = save_path $
           else uvf_froot = file_dirname(datafile, /mark_directory)
           uvf_savefilebase = [strmid(infilebase[0], 0, temp2[0]), strmid(infilebase[1], 0, temp2[1])] + file_label
        endif
     endelse
  endif 

  if n_elements(savefilebase_in) eq 1 then begin
     if n_elements(save_path) gt 0 then froot = save_path $
     else begin
        temp = file_dirname(savefilebase_in, /mark_directory)
        if temp ne '.' then froot = temp else froot = file_dirname(datafile[0], /mark_directory)
     endelse
     savefilebase = file_basename(savefilebase_in)
     general_filebase = savefilebase
  endif

  if n_elements(uvf_savefilebase_in) eq nfiles then begin
     if n_elements(save_path) gt 0 then uvf_froot = save_path $
     else begin
        temp = file_dirname(uvf_savefilebase_in, /mark_directory)
        for i=0, nfiles-1 do if temp[i] ne '.' then uvf_froot = temp[i] else uvf_froot = file_dirname(datafile[i], /mark_directory)
     endelse
     uvf_savefilebase = file_basename(uvf_savefilebase_in)
  endif
      
  uvf_savefile = uvf_froot + uvf_savefilebase + '_uvf.idlsave'
  uf_savefile = uvf_froot + uvf_savefilebase + '_uf_plane.idlsave'
  vf_savefile = uvf_froot + uvf_savefilebase + '_vf_plane.idlsave'
  uv_savefile = uvf_froot + uvf_savefilebase + '_uv_plane.idlsave'

  kcube_savefile = froot + savefilebase + '_kcube.idlsave'
  power_savefile = froot + savefilebase + '_power.idlsave'
  fits_power_savefile = froot + savefilebase + '_power.fits'

  if n_elements(weight_savefilebase_in) eq 0 then begin
     if n_elements(save_path) gt 0 then wt_froot = save_path $
     else wt_froot = file_dirname(weightfile, /mark_directory)
     wt_infilebase = file_basename(weightfile)
     temp2 = strpos(wt_infilebase, '.', /reverse_search)
     if nfiles eq 1 then weight_savefilebase = strmid(wt_infilebase, 0, temp2) + wt_file_label $
     else begin
        weight_savefilebase = strarr(nfiles)
        for i=0, nfiles-1 do weight_savefilebase[i] = strmid(wt_infilebase[i], 0, temp2[i]) + wt_file_label
     endelse
  endif else begin
     if n_elements(save_path) gt 0 then wt_froot = save_path $
     else begin
        temp = file_dirname(weight_savefilebase_in, /mark_directory)
        if nfiles eq 1 then if temp ne '.' then wt_froot = temp else  wt_froot = file_dirname(weightfile, /mark_directory) $
        else begin
           wt_froot = strarr(nfiles)
           for i=0, nfiles-1 do if temp[i] ne '.' then wt_froot[i] = temp[i] else $
              wt_froot[i] = file_dirname(weightfile[i], /mark_directory)
        endelse
     endelse
     weight_savefilebase = file_basename(weight_savefilebase_in)
  endelse
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'

  for j=0, nfiles-1 do begin
     file_obj = obj_new('idl_savefile', datafile[j])
     varnames = file_obj->names()

     wh_nside = where(strlowcase(varnames) eq 'nside', count_nside)
     if j gt 0 then if (count_nside gt 0 and healpix eq 0) or (count_nside eq 0 and healpix eq 1) then $
        message, 'One datafile is in healpix and the other is not.'
     if count_nside gt 0 then begin
        if j eq 0 then file_obj->restore, 'nside' else begin
           nside1 = nside
           file_obj->restore, 'nside'
           if nside1 ne nside then message, 'nside parameter does not agree between datafiles'
           undefine, nside1
        endelse
        healpix = 1
        
        if n_elements(pixelfile) ne nfiles then begin
           if j eq 0 then pixelfile = datafile[j] else pixelfile = [pixelfile, datafile[j]]
        endif
        if j eq 0 then begin
           if n_elements(pixelvar) eq 0 then pixel_varname = replicate('hpx_inds', nfiles) else pixel_varname =  pixelvar
        endif
     endif else healpix = 0

     wh_obs = where(strlowcase(varnames) eq 'obs', count_obs)
     if count_obs ne 0 then file_obj->restore, 'obs' $
     else begin
        wh_obs = where(strlowcase(varnames) eq 'obs_arr', count_obs)
        if count_obs ne 0 then begin
           file_obj->restore, 'obs_arr'
           if size(obs_arr,/type) eq 10 then begin
              n_obs = n_elements(obs_arr)
              
              max_baseline_vals = dblarr(n_obs)
              obs_radec_vals = dblarr(n_obs, 2)
              zen_radec_vals = dblarr(n_obs, 2)
              for i=0, n_obs-1 do begin
                 if not healpix then if abs((*obs_arr[i]).degpix - (*obs_arr[0]).degpix) gt 0 then $
                    message, 'inconsistent degpix values in obs_arr'
                 if total(abs((*obs_arr[i]).freq - (*obs_arr[0]).freq)) gt 0 then message, 'inconsistent freq values in obs_arr'
                 if abs((*obs_arr[i]).n_freq - (*obs_arr[0]).n_freq) gt 0 then message, 'inconsistent n_freq values in obs_arr'
                 
                 max_baseline_vals[i] = (*obs_arr[i]).max_baseline
                 obs_radec_vals[i, *] = [(*obs_arr[i]).obsra, (*obs_arr[i]).obsdec]
                 zen_radec_vals[i, *] = [(*obs_arr[i]).zenra, (*obs_arr[i]).zendec]              
              endfor

              if j eq 0 then max_baseline_lambda = max(max_baseline_vals) $
              else max_baseline_lambda = max([max_baseline_lambda, max_baseline_vals])
           
              if not healpix then $
                 if j eq 0 then degpix = (*obs_arr[0]).degpix else if (*obs_arr[0]).degpix ne degpix then $
                    message, 'degpix does not agree between datafiles'
              if j eq 0 then freq = (*obs_arr[0]).freq else if total(abs(freq-(*obs_arr[0]).freq)) ne 0 then $
                 message, 'frequencies do not agree between datafiles'
              if j eq 0 then n_freq = (*obs_arr[0]).n_freq else if (*obs_arr[0]).n_freq ne n_freq then $
                 message, 'n_freq does not agree between datafiles'
stop
           endif else begin
              n_obs = n_elements(obs_arr)
              if j eq 0 then max_baseline_lambda = max(obs_arr.max_baseline) $
              else max_baseline_lambda = max([max_baseline_lambda, obs_arr.max_baseline])
              
              obs_radec_vals = [[obs_arr.obsra],[obs_arr.obsdec]]
              zen_radec_vals = [[obs_arr.zenra],[obs_arr.zendec]]
              
              if total(abs(obs_arr.n_freq - obs_arr[0].n_freq)) ne 0 then message, 'inconsistent number of frequencies in obs_arr'
              if j eq 0 then n_freq = obs_arr[0].n_freq else if obs_arr[0].n_freq ne n_freq then $
                 message, 'n_freq does not agree between datafiles'
              
              if not healpix then begin
                 if total(abs(obs_arr.degpix - obs_arr[0].degpix)) ne 0 then message, 'inconsistent degpix values in obs_arr'
                 if j eq 0 then degpix = obs_arr[0].degpix else if obs_arr[0].degpix ne degpix  then $
                    message, 'degpix does not agree between datafiles'
              endif

              obs_tags = tag_names(obs_arr)
              wh_freq = where(strlowcase(obs_tags) eq 'freq', count_freq)
              if count_freq ne 0 then freq_vals = obs_arr.freq $
              else begin
                 freq_vals = dblarr(n_freq, n_obs)
                 for i=0, n_obs-1 do freq_vals[*,i] = (*obs_arr[i].bin).freq
              endelse
              if total(abs(freq_vals - rebin(freq_vals[*,0], n_freq, n_obs))) ne 0 then message, 'inconsistent freq values in obs_arr'
              if j eq 0 then freq = freq_vals[*,0] else if total(abs(freq - freq_vals[*,0])) ne 0 then $
                 message, 'frequencies do not agree between datafiles'
              freq_resolution = freq[1]-freq[0]
              
              dt_vals = dblarr(n_obs)
              for i=0, n_obs-1 do begin
                 times = (*obs_arr[i].baseline_info).jdate
                 ;; only allow time resolutions of n*.5 sec
                 dt_vals[i] = round(((times[1]-times[0])*24*3600)*2.)/2.
              endfor
              if total(abs(dt_vals - dt_vals[0])) ne 0 then message, 'inconsistent time averaging in obs_arr'
              if j eq 0 then time_resolution = dt_vals[0] else $
                 if total(abs(time_resolution - dt_vals[0])) ne 0 then message, 'time averaging does not agree between datafiles' 

           endelse

           theta_vals = angle_difference(obs_radec_vals[*,1], obs_radec_vals[*,0], zen_radec_vals[*,1], zen_radec_vals[*,0], $
                                     /degree, /nearest)
           if j eq 0 then max_theta = max(theta_vals) else max_theta = max([max_theta, theta_vals])

        endif else message, 'no obs or obs_arr in datafile'
     endelse
     
     wh_navg = where(strlowcase(varnames) eq 'n_avg', count_obs)
     if count_obs ne 0 then begin
        if j eq 0 then file_obj->restore, 'n_avg' else begin
           n_avg1 = n_avg
           file_obj->restore, 'n_avg'
           if n_avg1 ne n_avg then message, 'n_avg parameter does not agree between datafiles'
           undefine, n_avg1
        endelse
     endif else begin
        print, 'no n_avg present, assuming n_avg=32'
        if j eq 0 then n_avg = 32 else if n_avg ne 32 then message, 'n_avg parameter does not agree between datafiles'
     endelse   

     obj_destroy, file_obj
  endfor


  n_freqbins = n_freq / n_avg
  frequencies = dblarr(n_freqbins)
  for i=0, n_freqbins-1 do begin
     frequencies[i] = mean(freq[i*n_avg:i*n_avg+(n_avg-1)]) / 1e6 ;; in MHz
  endfor
  
  ;; the max baseline in the obs structure is given in wavelengths, need to convert using the maximum frequency
  max_baseline = 3e8/max(freq)*max_baseline_lambda ;;else max_baseline = 342.497


  if healpix then if n_elements(hpx_dftsetup_savefile) eq 0 then $
     hpx_dftsetup_savefile = froot + general_filebase + '_dftsetup.idlsave'
  

  if healpix then begin
     file_struct = {datafile: datafile, weightfile: weightfile, variancefile:variancefile, pixelfile:pixelfile, $
                    datavar:data_varname, variancevar:variance_varname, weightvar:weight_varname, pixelvar:pixel_varname, $
                    frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
                    max_baseline:max_baseline, max_theta:max_theta, nside:nside, $
                    hpx_dftsetup_savefile:hpx_dftsetup_savefile, $
                    uvf_savefile:uvf_savefile, uvf_weight_savefile:uvf_weight_savefile, $
                    uf_savefile:uf_savefile, vf_savefile:vf_savefile, uv_savefile:uv_savefile, kcube_savefile:kcube_savefile, $
                    power_savefile:power_savefile, fits_power_savefile:fits_power_savefile, $
                    savefile_froot:froot, savefilebase:savefilebase, general_filebase:general_filebase, $
                    weight_savefilebase:weight_savefilebase}
  endif else begin
     file_struct = {datafile: datafile, weightfile: weightfile, variancefile:variancefile, $
                    datavar:data_varname, weightvar:weight_varname, variancevar:variance_varname, $
                    frequencies:frequencies, freq_resolution:freq_resolution, time_resolution:time_resolution, $
                    max_baseline:max_baseline, max_theta:max_theta, degpix:degpix, $
                    kcube_savefile:kcube_savefile, power_savefile:power_savefile, fits_power_savefile:fits_power_savefile, $
                    savefile_froot:froot, savefilebase:savefilebase, general_filebase:general_filebase, $
                    weight_savefilebase:weight_savefilebase}
  endelse

return, file_struct

end
