function fhd_file_setup, datafile, datavar, weightfile, weightvar, frequencies, max_baseline, $
                         degpix = degpix, nside = nside, pixelfile = pixelfile, pixelvar = pixelvar, $
                         hpx_dftsetup_savefile = hpx_dftsetup_savefile, savefilebase = savefilebase_in, $
                         weight_savefilebase = weight_savefilebase_in, uvf_savefilebase = uvf_savefilebase_in

  if n_params() ne 6 then message, 'Wrong number of parameters passed. Required parameters are: datafile (string or string array), ' $
                                   + 'datavar (string), weightfile (string or string array), weightvar (string), ' + $
                                   'freqencies (float or double array), max_baseline (float or double)'

  if n_elements(datavar) ne 1 then message, 'must have one datavar allowed'
  if n_elements(weightvar) ne 1 then message, 'must have one weightvar allowed'
  if n_elements(max_baseline) gt 1 then message, 'only one max_baseline allowed'
  if n_elements(pixelfile) gt 1 then message, 'only one pixelfile allowed'
  if n_elements(pixelvar) gt 1 then message, 'only one pixelvar allowed'
  if n_elements(degpix) gt 1 then message, 'only one degpix allowed'
  if n_elements(nside) gt 1 then message, 'only one nside allowed'
  if n_elements(hpx_dftsetup_savefile) gt 1 then message, 'only one hpx_dftsetup_savefile allowed'
  if n_elements(savefilebase) gt 1 then message, 'only one savefilebase allowed'

  freq_dims = size(frequencies, /dimension)
  if n_elements(freq_dims) gt 1 then message, 'frequencies must be a 1d array'

  nfiles = n_elements(datafile)
  if nfiles gt 2 then message, 'only 1 or 2 datafiles is supported'
  datafile_test = file_test(datafile)
  if min(datafile_test) eq 0 then message, 'datafile not found'
  if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
     print, 'datafiles are identical'
     datafile = datafile[0]
  endif

  if n_elements(weightfile) ne nfiles then message, 'same number of weight and data files must be supplied'
  weightfile_test = file_test(weightfile)
  if min(weightfile_test) eq 0 then message, 'weightfile not found'

  if n_elements(weight_savefilebase) gt 0 and n_elements(weight_savefilebase) ne nfiles then $
     message, 'if weight_savefilebase is specified it must have the same number of elements as data files'
  if n_elements(uvf_savefilebase) gt 0 and n_elements(uvf_savefilebase) ne nfiles then $
     message, 'if uvf_savefilebase is specified it must have the same number of elements as data files'

  if n_elements(savefilebase_in) eq 0 or n_elements(uvf_savefilebase_in) lt nfiles then begin
     if nfiles eq 1 then begin
        froot = file_dirname(datafile, /mark_directory)
        uvf_froot = froot
        infilebase = file_basename(datafile)
        temp2 = strpos(infilebase, '.', /reverse_search)
        if n_elements(savefilebase_in) eq 0 then savefilebase = strmid(infilebase, 0, temp2) else savefilebase = savefilebase_in
        
        ;; if we're only dealing with one file and uvf_savefilebase isn't specified then use same base for uvf files 
        if n_elements(uvf_savefilebase_in) eq 0 then uvf_savefilebase = savefilebase
     endif else begin
        froot = file_dirname(datafile[0], /mark_directory)
        infilebase = file_basename(datafile)
        temp2 = strpos(infilebase, '.', /reverse_search)
        
        if n_elements(savefilebase_in) eq 0 then begin
           fileparts_1 = strsplit(strmid(infilebase[0], 0, temp2[0]), '_', /extract)
           fileparts_2 = strsplit(strmid(infilebase[1], 0, temp2[1]), '_', /extract)
           match_test = strcmp(fileparts_1, fileparts_2)
           wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
           if count_diff eq 0 then savefilebase = strmid(infilebase[0], 0, temp2[0]) + '_joint' $
           else begin
              if count_same gt 0 then savefilebase = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) $
                                                     + '_' + strjoin(fileparts_2[wh_diff]) + '_joint' $
              else savefilebase = infilebase[0] + infilebase[1] + '_joint'
           endelse
        endif
        
        if n_elements(uvf_savefilebase_in) eq 0 then begin
           ;; need 2 uvf files
           uvf_froot = file_dirname(datafile, /mark_directory)
           uvf_savefilebase = [strmid(infilebase[0], 0, temp2[0]), strmid(infilebase[1], 0, temp2[1])]
        endif
     endelse
  endif 

  if n_elements(savefilebase_in) eq 1 then begin
     temp = file_dirname(savefilebase_in, /mark_directory)
     if temp ne '.' then froot = temp else froot = file_dirname(datafile[0], /mark_directory)
     savefilebase = file_basename(savefilebase_in)
  endif

  if n_elements(uvf_savefilebase_in) eq nfiles then begin
     temp = file_dirname(uvf_savefilebase_in, /mark_directory)
     for i=0, nfiles-1 do if temp[i] ne '.' then uvf_froot = temp[i] else uvf_froot = file_dirname(datafile[i], /mark_directory)
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
     wt_froot = file_dirname(weightfile, /mark_directory)
     wt_infilebase = file_basename(weightfile)
     temp2 = strpos(wt_infilebase, '.', /reverse_search)
     weight_savefilebase = strarr(nfiles)
     for i=0, nfiles-1 do weight_savefilebase[i] = strmid(wt_infilebase[i], 0, temp2[i]) + '_weights'
    endif else begin
     temp = file_dirname(weight_savefilebase_in, /mark_directory)
     wt_froot = strarr(nfiles)
     for i=0, nfiles-1 do if temp[i] ne '.' then  wt_froot[i] = temp[i] else wt_froot[i] = file_dirname(weightfile[i], /mark_directory)
     weight_savefilebase = file_basename(weight_savefilebase_in)
  endelse
  
  uvf_weight_savefile = wt_froot + weight_savefilebase + '_uvf.idlsave'


  if n_elements(nside) ne 0 then healpix = 1 else healpix = 0

  if healpix then begin
     if n_elements(pixelfile) eq 0 or n_elements(pixelvar) eq 0 then message, $
        'pixelfile and pixelvar keywords must be specified for healpix'

     if file_test(pixelfile) eq 0 then message, 'pixelfile not found'

     if n_elements(hpx_dftsetup_savefile) eq 0 then begin
        pix_froot = file_dirname(pixelfile, /mark_directory)
        pix_infilebase = file_basename(pixelfile)
        temp2 = strpos(pix_infilebase, '.', /reverse_search)

        hpx_dftsetup_savefile = pix_froot + strmid(pix_infilebase, 0, temp2) + '_dftsetup.idlsave'
     endif 
  endif else begin

     if n_elements(degpix) eq 0 then message, 'degpix keyword must be specified unless in healpix'
  endelse
 

  if healpix then begin
     file_struct = {datafile: datafile, weightfile: weightfile, pixelfile:pixelfile, datavar:datavar, weightvar:weightvar, $
                    pixelvar:pixelvar, frequencies:frequencies, max_baseline:max_baseline, nside:nside, $
                    hpx_dftsetup_savefile:hpx_dftsetup_savefile, $
                    uvf_savefile:uvf_savefile, uvf_weight_savefile:uvf_weight_savefile, $
                    uf_savefile:uf_savefile, vf_savefile:vf_savefile, uv_savefile:uv_savefile, kcube_savefile:kcube_savefile, $
                    power_savefile:power_savefile, fits_power_savefile:fits_power_savefile, $
                    savefile_froot:froot, savefilebase:savefilebase}
  endif else begin
     file_struct = {datafile: datafile, weightfile: weightfile, datavar:datavar, weight_var:weight_var, $
                    frequencies:frequencies, max_baseline:max_baseline, degpix:degpix, kcube_savefile:kcube_savefile, $
                    power_savefile:power_savefile, fits_power_savefile:fits_power_savefile, $
                    savefile_froot:froot, savefilebase:savefilebase}
  endelse

return, file_struct

end
