function fhd_file_setup, filename, weightfile = weightfile, $
    variancefile = variancefile, beamfile = beamfile, pixelfile = pixelfile, $
    dirtyvar = dirtyvar, modelvar = modelvar, weightvar = weightvar, $
    variancevar = variancevar, beamvar = beamvar, pixelvar = pixelvar, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, $
    save_path = save_path, savefilebase = savefilebase_in, $
    uvf_input = uvf_input, sim = sim, refresh_info = refresh_info, $
    uvf_options = uvf_options, ps_options = ps_options

  ;; check to see if filename is an info file
  wh_info = strpos(filename, '_info.idlsave')
  if max(wh_info) gt -1 then begin
    if n_elements(filename) gt 1 then message, 'only 1 info file can be passed'
    info_file = filename

    restore, info_file

    if n_elements(save_path) ne 0 then begin
      froot = save_path
    endif else begin
      info_filebase = cgRootName(info_file, directory=froot)
    endelse

    if not tag_exist(metadata_struct, 'nfiles') then begin
      print, 'Info file is very old, creating a new one.'
      refresh_info = 1
    endif

    if keyword_set(refresh_info) then begin
      if n_elements(metadata_struct) gt 0 then begin
        datafile = metadata_struct.datafile
      endif else begin
        datafile = file_struct_arr[0].datafile
      endelse

      datafile_test = file_test(datafile)
      if min(datafile_test) eq 0 then begin
        ;; can't find datafile. try some other potential folders
        if n_elements(save_path) ne 0 then begin
          datafile_test = file_test(save_path + file_basename(datafile))
          if min(datafile_test) gt 0 then datafile = save_path + file_basename(datafile)
        endif

        if min(datafile_test) eq 0 then begin
          datafile_test = file_test(file_dirname(info_file, /mark_directory) + $
            file_basename(datafile))
          if min(datafile_test) gt 0 then begin
            datafile = file_dirname(info_file, /mark_directory) + file_basename(datafile)
          endif
        endif

        if min(datafile_test) eq 0 then begin
          message, 'refresh_info is set but datafile cannot be found'
        endif
      endif

    endif
  endif

  if keyword_set(refresh_info) or max(wh_info) eq -1 then begin
    if n_elements(datafile) eq 0 then datafile = filename

    ;; check whether datafiles have polarization identifiers in filename
    pol_exist = stregex(datafile, '[xy][xy]', /boolean, /fold_case)
    if max(pol_exist) gt 0 then begin
      wh_pol_exist = where(pol_exist gt 0, count_pol_exist, ncomplement = count_no_pol)
      if count_no_pol gt 0 then begin
        message, 'some datafiles have pol identifiers and some do not'
      endif
      pols = strlowcase(stregex(datafile, '[xy][xy]', /extract, /fold_case))
      if n_elements(weightfile) gt 0 then begin
        wt_pols = strlowcase(stregex(weightfile, '[xy][xy]', /extract, /fold_case))
        match, pols, wt_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then begin
          message, 'weightfile polarizations must match datafile polarizations'
        endif
      endif
      if n_elements(variancefile) gt 0 then begin
        var_pols = strlowcase(stregex(variancefile, '[xy][xy]', /extract, /fold_case))
        match, pols, var_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then begin
          message, 'variancefile polarizations must match datafile polarizations'
        endif
      endif
      if n_elements(beamfile) gt 0 then begin
        bm_pols = strlowcase(stregex(beamfile, '[xy][xy]', /extract, /fold_case))
        match, pols, bm_pols, suba, subb, count = count_pol_match
        if count_pol_match ne n_elements(pols) then begin
          message, 'beamfile polarizations must match datafile polarizations'
        endif
      endif

      pol_inc = pols[0]
      pol_num = intarr(n_elements(pols))
      for i=0, n_elements(pols)-1 do begin
        wh_pol = where(pol_inc eq pols[i], count_pol)
        if count_pol eq 1 then pol_num[i] = wh_pol[0] else begin
          pol_inc = [pol_inc, pols[i]]
          pol_num[i] = n_elements(pol_inc)-1
        endelse
      endfor

      npol = n_elements(pol_inc)
      nfiles = ceil(n_elements(datafile)/float(npol))
      if nfiles gt 2 then begin
        message, 'Only 1 or 2 datafiles supported per pol'
      endif
      if n_elements(datafile) ne npol*nfiles then begin
        message, 'The same number of files must be included per pol'
      endif
      temp = strarr(npol, nfiles)
      temp_wt = strarr(npol, nfiles)
      temp_var = strarr(npol, nfiles)
      temp_bm = strarr(npol, nfiles)
      for i=0, npol-1 do begin
        wh_pol =  where(pol_num eq i, count_pol)
        if count_pol ne nfiles then begin
          message, 'The same number of files must be included per pol'
        endif
        temp[i,*] = datafile[wh_pol]
        if n_elements(weightfile) gt 0 then begin
          temp_wt[i,*] = weightfile[wh_pol]
        endif else begin
          temp_wt = temp
        endelse
        if n_elements(variancefile) gt 0 then begin
          temp_var[i,*] = variancefile[wh_pol]
        endif else begin
          temp_var = temp
        endelse
        if n_elements(beamfile) gt 0 then begin
          temp_bm[i,*] = beamfile[wh_pol]
        endif else begin
          temp_bm = temp
        endelse
        if n_elements(pixelfile) gt 0 then begin
          temp_pix[i,*] = pixelfile[wh_pol]
        endif else begin
          temp_pix = temp
        endelse
      endfor
      datafile = temporary(temp)
      weightfile = temporary(temp_wt)
      variancefile = temporary(temp_var)
      beamfile = temporary(temp_bm)
      pixelfile = temporary(temp_pix)
    endif else begin
      ;; no pol identifiers
      if n_elements(pol_inc) eq 0 then pol_inc = ['xx', 'yy']
      pol_enum = ['xx', 'yy']
      npol = n_elements(pol_inc)
      pol_num = intarr(npol)
      for i=0, npol-1 do begin
        wh = where(pol_enum eq pol_inc[i], count)
        if count eq 0 then begin
          message, 'pol ' + pol_inc[i] + ' not recognized.'
        endif
        pol_num[i] = wh[0]
      endfor
      pol_inc = pol_enum[pol_num[uniq(pol_num, sort(pol_num))]]

      nfile_dims = size(datafile, /dimension)
      if n_elements(nfile_dims) gt 1 then datafile = datafile[0,*]
      nfiles = n_elements(datafile)
      if nfiles gt 2 then message, 'Only 1 or 2 datafiles supported'
      if nfiles eq 2 then if datafile[0] eq datafile[1] then begin
        print, 'datafiles are identical'
        datafile = datafile[0]
        nfiles = 1
      endif


      if n_elements(weightfile) gt 0 and n_elements(weightfile) ne nfiles then begin
        message, 'weightfile must have the same number of elements as datafile'
      endif
      if n_elements(variancefile) gt 0 and n_elements(variancefile) ne nfiles then begin
        message, 'variancefile must have the same number of elements as datafile'
      endif
      if n_elements(beamfile) gt 0 and n_elements(beamfile) ne nfiles then begin
        message, 'beamfile must have the same number of elements as datafile'
      endif

      temp = strarr(npol, nfiles)
      temp_wt = strarr(npol, nfiles)
      temp_var = strarr(npol, nfiles)
      temp_bm = strarr(npol, nfiles)
      for i=0, nfiles-1 do begin
        temp[*,i] = datafile[i]
        if n_elements(weightfile) gt 0 then begin
          temp_wt[*,i] = weightfile[i]
        endif else begin
          temp_wt = temp
        endelse
        if n_elements(variancefile) gt 0 then begin
          temp_var[*,i] = variancefile[i]
        endif else begin
          temp_var = temp
        endelse
        if n_elements(beamfile) gt 0 then begin
          temp_bm[*,i] = beamfile[i]
        endif else begin
          temp_bm = temp
        endelse
        if n_elements(pixelfile) gt 0 then begin
          temp_pix[*,i] = pixelfile[i]
        endif else begin
          temp_pix = temp
        endelse
      endfor
      datafile = temp
      weightfile = temporary(temp_wt)
      variancefile = temporary(temp_var)
      beamfile = temporary(temp_bm)
      pixelfile = temporary(temp_pix)
    endelse

    if n_elements(savefilebase_in) gt 1 then message, 'only one savefilebase allowed'

    if nfiles eq 1 then infile_label = '' else begin
      data_filebase = [cgRootName(datafile[0,0]), cgRootName(datafile[0,1])]

      ;; make this general to all pols
      if max(pol_exist) gt 0 then begin
        for i=0, 1 do begin
          data_filebase[i] = strjoin(strsplit(data_filebase[i], '_?[xy][xy]_?', $
            /regex, /extract), '_')
        endfor
      endif

      fileparts_1 = strsplit(data_filebase[0], '_', /extract)
      fileparts_2 = strsplit(data_filebase[1], '_', /extract)
      match_test = strcmp(fileparts_1, fileparts_2)
      wh_diff = where(match_test eq 0, count_diff, complement = wh_same, $
        ncomplement = count_same)

      if count_diff gt 0 then begin
        infile_label = [strjoin(fileparts_1[wh_diff]), strjoin(fileparts_2[wh_diff])]
      endif else begin
        infile_label = strarr(2)
      endelse
    endelse

    if n_elements(savefilebase_in) eq 0 then begin
      if nfiles eq 1 then begin
        data_filebase = cgRootName(datafile[0], directory=datafile_dir)

        ;; make this general to all pols
        if max(pol_exist) gt 0 then begin
          data_filebase= strjoin(strsplit(data_filebase, '_?[xy][xy]_?', $
            /regex, /extract), '_')
        endif

        if n_elements(save_path) ne 0 then begin
          froot = save_path
        endif else begin
          froot = datafile_dir
        endelse

        general_filebase = data_filebase

      endif else begin
        if n_elements(save_path) ne 0 then begin
          froot = save_path
        endif else begin
          temp = cgRootName(datafile[0], directory=froot)
        endelse

        if count_diff eq 0 then general_filebase = data_filebase[0] + '_joint' else begin
          if count_same gt 0 then begin
            general_filebase = strjoin(fileparts_1[wh_same], '_') + '__' + $
              strjoin(fileparts_1[wh_diff]) + '_' + strjoin(fileparts_2[wh_diff]) + '_joint'
            endif else begin
              general_filebase = data_filebase[0] + data_filebase[1] + '_joint'
            endelse
        endelse

      endelse
    endif else begin
      if n_elements(save_path) gt 0 then froot = save_path else begin
        savefilebase_in_base = cgRootName(savefilebase_in[0], directory=froot)
        if froot eq '.' then temp = cgRootName(datafile[0], directory=froot)
      endelse

      general_filebase = savefilebase_in_base
    endelse

    if not keyword_set(refresh_info) then begin
      ;; check for pre-saved info file. If it exists, restore it and return structure.
      ;; Otherwise construct structure.
      if keyword_set(uvf_input) then uvf_info_tag = '_uvf' else uvf_info_tag = ''
      info_file = froot + general_filebase + uvf_info_tag + '_info.idlsave'
      info_filetest = file_test(info_file)
      if info_filetest eq 1 then begin
        ;; check if info file has the right polarizations
        info_metadata = getvar_savefile(info_file, 'metadata_struct')
        match2, pol_inc, info_metadata.pol_inc, suba, subb
        if min([suba, subb]) ge 0 then begin
          ;; looks good, restore & check for directory structure changes
          file_struct_arr = fhd_file_setup(info_file, weightfile = weightfile, $
              variancefile = variancefile, beamfile = beamfile, pixelfile = pixelfile, $
              dirtyvar = dirtyvar, modelvar = modelvar, weightvar = weightvar, $
              variancevar = variancevar, beamvar = beamvar, pixelvar = pixelvar, $
              freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
              freq_flag_name = freq_flag_name, $
              save_path = save_path, savefilebase = savefilebase_in, $
              uvf_input = uvf_input, sim = sim, refresh_info = refresh_info, $
              uvf_options = uvf_options, ps_options = ps_options)

          return, file_struct_arr
        endif
      endif
    endif

    datafile_test = file_test(datafile)
    if min(datafile_test) eq 0 then message, 'datafile not found'

    void = getvar_savefile(datafile[0], names = varnames)

    if keyword_set(sim) then begin
      if n_elements(modelvar) eq 0 then begin
        if max(strmatch(varnames, 'model*',/fold_case)) then begin
          if keyword_set(uvf_input) then begin
            model_varname = strupcase('model_uv_arr')
          endif else begin
            if max(pol_exist) gt 0 then begin
              model_varname = strupcase('model_cube')
            endif else begin
              model_varname = strupcase('model_' + pol_inc + '_cube')
            endelse
          endelse
        endif else begin
          print, 'No model cubes found, checking for dirty or residual cubes'
          if max(strmatch(varnames, 'dirty*',/fold_case)) then begin
            if keyword_set(uvf_input) then begin
              model_varname = strupcase('dirty_uv_arr')
            endif else begin
              if max(pol_exist) gt 0 then begin
                model_varname = strupcase('dirty_cube')
              endif else begin
                model_varname = strupcase('dirty_' + pol_inc + '_cube')
              endelse
            endelse
          endif else if max(strmatch(varnames, 'res*',/fold_case)) then begin
            if keyword_set(uvf_input) then begin
              model_varname = strupcase('res_uv_arr')
            endif else begin
              if max(pol_exist) gt 0 then begin
                model_varname = strupcase('res_cube')
              endif else begin
                model_varname = strupcase('res_' + pol_inc + '_cube')
              endelse
            endelse
          endif
        endelse
      endif else model_varname = modelvar

      if n_elements(model_varname) ne npol then begin
        if n_elements(model_varname) eq 1 then begin
          model_varname = replicate(model_varname, npol)
        endif
      endif else begin
        message, 'modelvar must be a scalar or have the same number of elements as pol_inc'
      endelse

      type_inc = ['model']
      if npol gt 1 then begin
        cube_varname = transpose(model_varname)
      endif else begin
        cube_varname = model_varname
      endelse
      ntypes = n_elements(type_inc)
      type_pol_str = strarr(npol, ntypes)
      for i=0, npol-1 do type_pol_str[i, *] = type_inc + '_' + pol_inc[i]
    endif else begin
      ;; check for the existence of dirty, model, residual cubes
      if max(strmatch(varnames, 'dirty*',/fold_case)) then begin
        if n_elements(dirtyvar) eq 0 then begin
          if keyword_set(uvf_input) then begin
            dirty_varname = strupcase('dirty_uv_arr')
          endif else begin
            if max(pol_exist) gt 0 then begin
              dirty_varname = strupcase('dirty_cube')
            endif else begin
              dirty_varname = strupcase('dirty_' + pol_inc + '_cube')
            endelse
          endelse
        endif else dirty_varname = dirtyvar
        if n_elements(dirty_varname) ne npol then begin
          if n_elements(dirty_varname) eq 1 then begin
            dirty_varname = replicate(dirty_varname, npol)
          endif else begin
            message, 'dirtyvar must be a scalar or have the same number of elements as pol_inc'
          endelse
        endif

        type_inc = ['dirty']
        if npol gt 1 then begin
          cube_varname = transpose(dirty_varname)
        endif else begin
          cube_varname = dirty_varname
        endelse
      endif

      if max(strmatch(varnames, 'model*',/fold_case)) then begin
        if n_elements(modelvar) eq 0 then begin
          if keyword_set(uvf_input) then begin
            model_varname = strupcase('model_uv_arr')
          endif else begin
            if max(pol_exist) gt 0 then begin
              model_varname = strupcase('model_cube')
            endif else begin
              model_varname = strupcase('model_' + pol_inc + '_cube')
            endelse
          endelse
        endif else model_varname = modelvar
        if n_elements(model_varname) ne npol then begin
          if n_elements(model_varname) eq 1 then begin
            model_varname = replicate(model_varname, npol)
          endif else begin
            message, 'modelvar must be a scalar or have the same number of elements as pol_inc'
          endelse
        endif
        if n_elements(type_inc) eq 0 then begin
          type_inc = ['model']
          if npol gt 1 then begin
            cube_varname = transpose(model_varname)
          endif else begin
            cube_varname = model_varname
          endelse
        endif else begin
          type_inc = [type_inc, 'model']
          if npol gt 1 then begin
            cube_varname = [cube_varname,transpose(model_varname)]
          endif else begin
            cube_varname = [cube_varname, model_varname]
          endelse
        endelse
      endif

      if max(strmatch(varnames, 'res*',/fold_case)) then begin
        if n_elements(residualvar) eq 0 then begin
          if keyword_set(uvf_input) then begin
            residual_varname = strupcase('res_uv_arr')
          endif else begin
            if max(pol_exist) gt 0 then begin
              residual_varname = strupcase('res_cube')
            endif else begin
              residual_varname = strupcase('res_' + pol_inc + '_cube')
            endelse
          endelse
        endif else dirty_varname = residualvar
        if n_elements(residual_varname) ne npol then begin
          if n_elements(residual_varname) eq 1 then begin
            residual_varname = replicate(residual_varname, npol)
          endif
        endif else begin
          message, 'residualvar must be a scalar or have the same number of elements as pol_inc'
        endelse

        if n_elements(type_inc) eq 0 then begin
          type_inc = ['res']
          if npol gt 1 then begin
            cube_varname = transpose(residual_varname)
          endif else begin
            cube_varname = residual_varname
          endelse
        endif else begin
          type_inc = [type_inc, 'res']
          if npol gt 1 then begin
            cube_varname = [cube_varname,transpose(residual_varname)]
          endif else begin
            cube_varname = [cube_varname, residual_varname]
          endelse
        endelse
      endif else if (Max(type_inc EQ 'dirty')<Max(type_inc EQ 'model')) then begin
        ;; residual can be constructed from dirty-model
        type_inc = [type_inc, 'res']
        if npol gt 1 then begin
          cube_varname = [cube_varname,transpose(strarr(npol))]
        endif else begin
          cube_varname = [cube_varname, '']
        endelse
      endif
      ntypes = n_elements(type_inc)
      type_pol_str = strarr(npol, ntypes)
      for i=0, npol-1 do type_pol_str[i, *] = type_inc + '_' + pol_inc[i]
    endelse

    if n_elements(weightvar) eq 0 then begin
      if keyword_set(uvf_input) then begin
        weight_varname = strupcase('weights_uv_arr')
      endif else begin
        if max(pol_exist) gt 0 then begin
          weight_varname = strupcase('weights_cube')
        endif else begin
          weight_varname = strupcase('weights_' + pol_inc + '_cube')
        endelse
      endelse
    endif else weight_varname = weightvar
    if n_elements(weight_varname) ne npol then begin
      if n_elements(weight_varname) eq 1 then begin
        weight_varname = replicate(weight_varname, npol)
      endif else begin
        message, 'weightvar must be a scalar or have the same number of elements as pol_inc'
      endelse
    endif
    if n_elements(variancevar) eq 0 then begin
      if keyword_set(uvf_input) then begin
        variance_varname = strupcase('variance_uv_arr')
      endif else begin
        if max(pol_exist) gt 0 then begin
          variance_varname = strupcase('variance_cube')
        endif else begin
          variance_varname = strupcase('variance_' + pol_inc + '_cube')
        endelse
      endelse
    endif else variance_varname = variancevar
    if n_elements(variance_varname) ne npol then begin
      if n_elements(variance_varname) eq 1 then begin
        variance_varname = replicate(variance_varname, npol)
      endif else begin
        message, 'variancevar must be a scalar or have the same number of elements as pol_inc'
      endelse
    endif
    if n_elements(beamvar) eq 0 then begin
      if keyword_set(uvf_input) then begin
        if max(pol_exist) gt 0 then begin
          beam_varname = strupcase('beam2_image')
        endif else begin
          beam_varname = strupcase('beam2_' + pol_inc + '_image')
        endelse
      endif else begin
        if max(pol_exist) gt 0 then begin
          beam_varname = strupcase('beam_squared_cube')
        endif else begin
          beam_varname = strupcase('beam_' + pol_inc + '_cube')
        endelse
      endelse
    endif else beam_varname = beamvar
    if n_elements(beam_varname) ne npol then begin
      if n_elements(beam_varname) eq 1 then begin
        beam_varname = replicate(beam_varname, npol)
      endif else begin
        message, 'beamvar must be a scalar or have the same number of elements as pol_inc'
      endelse
    endif
    if n_elements(pixelvar) eq 0 then begin
      pixel_varname = strupcase('hpx_inds')
    endif else pixel_varname = pixelvar
    if n_elements(pixel_varname) ne npol then begin
      if n_elements(pixel_varname) eq 1 then begin
        pixel_varname = replicate(pixel_varname, npol)
      endif else begin
        message, 'pixelvar must be a scalar or have the same number of elements as pol_inc'
      endelse
    endif
    n_obs = lonarr(npol, nfiles)
    for j=0, nfiles*npol-1 do begin
      pol_i = j mod npol
      file_i = j / npol

      void = getvar_savefile(datafile[pol_i, file_i], names = varnames)

      wh_nside = where(strlowcase(varnames) eq 'nside', count_nside)
      if count_nside gt 0 then this_healpix = 1 else this_healpix = 0

      for type_i=0, ntypes-1 do begin
        if cube_varname[type_i, pol_i] ne '' then begin
          wh = where(strlowcase(varnames) eq strlowcase(cube_varname[type_i, pol_i]), count)
          if count eq 0 then begin
            message, cube_varname[type_i, pol_i] + ' is not present in datafile ' + $
              '(datafile=' + datafile[pol_i, file_i] + ')'
          endif

          data_size = getvar_savefile(datafile[pol_i, file_i], $
            cube_varname[type_i, pol_i], /return_size)

          if data_size[0] eq 0 then message, 'data cube has no size'
          if data_size[n_elements(data_size)-2] eq 10 then begin
            ;; data is a pointer
            if keyword_set(uvf_input) then begin
              if data_size[0] ne 2 then begin
                message, 'Data are in a pointer array, format unknown'
              endif
              if data_size[1] ne npol then begin
                message, 'Data are in a pointer array, format unknown'
              endif
              data = getvar_savefile(datafile[pol_i, file_i], cube_varname[type_i, pol_i])
              dims2 = size(*data[0], /dimension)
              iter=1
              while min(dims2) eq 0 and iter lt data_size[2] do begin
                print, 'warning: some frequency slices have null pointers'
                dims2 = size(*data[iter], /dimension)
                iter = iter+1
              endwhile
              if min(dims2) eq 0 then begin
                message, 'data cube is empty (file: ' + datafile[pol_i, file_i] + $
                  ', cube name: ' + cube_varname[type_i, pol_i] + ')'
              endif
              this_data_dims = [dims2, data_size[2], data_size[1]]
              undefine_fhd, data
            endif else begin
              message, 'Data is in a pointer array, format unknown'
            endelse

          endif else this_data_dims = data_size[1:data_size[0]]

          if type_i eq 0 and j eq 0 then begin
            data_dims = this_data_dims
          endif else if total(abs(this_data_dims - data_dims)) ne 0 then begin
            message, 'data dimensions in files do not match'
          endif
        endif
      endfor


      void = getvar_savefile(weightfile[pol_i, file_i], names = wt_varnames)
      void = getvar_savefile(variancefile[pol_i, file_i], names = var_varnames)
      if n_elements(beamfile) gt 0 then begin
        void = getvar_savefile(beamfile[pol_i, file_i], names = bm_varnames)
      endif
      wh = where(strlowcase(wt_varnames) eq strlowcase(weight_varname[pol_i]), count)

      if count eq 0 then begin
        message, weight_varname[pol_i] + ' is not present in weightfile (weightfile=' + $
          weightfile[pol_i, file_i] + ')'
      endif

      wt_size = getvar_savefile(weightfile[pol_i, file_i], weight_varname[pol_i], $
        /return_size)
      if wt_size[0] eq 0 then message, 'weight cube has no size'
      if wt_size[n_elements(wt_size)-2] eq 10 then begin
        ;; weights cube is a pointer
        if keyword_set(uvf_input) then begin
          if wt_size[0] ne 2 then begin
            message, 'Weights are in a pointer array, format unknown'
          endif
          if wt_size[1] ne npol then begin
            message, 'Weights are in a pointer array, format unknown'
          endif
          weights = getvar_savefile(weightfile[pol_i, file_i], weight_varname[pol_i])
          dims2 = size(*weights[0], /dimension)
          iter=1
          while min(dims2) eq 0 and iter lt data_size[2] do begin
            print, 'warning: some frequency slices have null pointers'
            dims2 = size(*weights[iter], /dimension)
            iter = iter+1
          endwhile
          wt_dims = [dims2, wt_size[2], wt_size[1]]
          undefine_fhd, weights
        endif else begin
          message, 'Weights are in a pointer array, format unknown'
        endelse

      endif else wt_dims = wt_size[1:wt_size[0]]

      if total(abs(wt_dims - data_dims)) ne 0 then begin
        message, 'weight and data dimensions in files do not match'
      endif

      wh = where(strlowcase(var_varnames) eq strlowcase(variance_varname[pol_i]), count)

      if count eq 0 then begin
        print, variance_varname[pol_i] + ' is not present in variancefile ' + $
          '(variancefile=' + variancefile[pol_i, file_i] + '). ' + $
          'Using weights^2 instead of variances'
        no_var = 1
      endif else begin
        no_var = 0
        var_size = getvar_savefile(variancefile[pol_i, file_i], variance_varname[pol_i], /return_size)
        if var_size[0] eq 0 then message, 'Variance cube has no size'
        if var_size[n_elements(var_size)-2] eq 10 then begin
          ;; Variance cube is a pointer
          if keyword_set(uvf_input) then begin
            if var_size[0] ne 2 then begin
              message, 'Variance are in a pointer array, format unknown'
            endif
            if var_size[1] ne npol then begin
              message, 'Variance are in a pointer array, format unknown'
            endif
            variance = getvar_savefile(variancefile[pol_i, file_i], variance_varname[pol_i])
            dims2 = size(*variance[0], /dimension)
            iter=1
            while min(dims2) eq 0 and iter lt data_size[2] do begin
              print, 'warning: some frequency slices have null pointers'
              dims2 = size(*variance[iter], /dimension)
              iter = iter+1
            endwhile
            var_dims = [dims2, var_size[2], var_size[1]]
            undefine_fhd, variance
          endif else begin
            message, 'Variance is in a pointer array, format unknown'
          endelse

        endif else var_dims = var_size[1:var_size[0]]

        if total(abs(var_dims - data_dims)) ne 0 then begin
          message, 'variance and data dimensions in files do not match'
        endif
      endelse

      if n_elements(beamfile) gt 0 then begin
        wh = where(strlowcase(bm_varnames) eq strlowcase(beam_varname[pol_i]), count)
        if count eq 0 then begin
          if beamfile[pol_i, file_i] ne datafile[pol_i, file_i] then begin
            message, beam_varname[pol_i] + ' is not present in beamfile (beamfile=' + $
              beamfile[pol_i, file_i] + ')'
          endif else begin
            undefine, beamfile
          endelse
        endif
      endif

      if n_elements(beamfile) gt 0 then begin
        bm_size = getvar_savefile(beamfile[pol_i, file_i], beam_varname[pol_i], /return_size)
        if bm_size[0] eq 0 then message, 'beam cube has no size'
        if bm_size[n_elements(bm_size)-2] eq 10 then begin
          ;; beam cube is a pointer
          message, 'Beam is in a pointer array, format unknown'
        endif else bm_dims = bm_size[1:bm_size[0]]

        if total(abs(bm_dims - data_dims)) ne 0 then begin
          print, 'Warning: Beam and data dimensions in files do not match'
        endif
      endif

      if j gt 0 then if (this_healpix eq 1 and healpix eq 0) or (this_healpix eq 0 and healpix eq 1) then begin
        message, 'One datafile is in healpix and the other is not.'
      endif
      if this_healpix eq 1 then begin
        if j eq 0 then nside = getvar_savefile(datafile[pol_i, file_i], 'nside') else begin
          nside1 = nside
          nside = getvar_savefile(datafile[pol_i, file_i], 'nside')
          if nside1 ne nside then begin
            message, 'nside parameter does not agree between datafiles'
          endif
          undefine, nside1
        endelse
        healpix = 1

        void = getvar_savefile(pixelfile[pol_i, file_i], names = pix_varnames)
        wh = where(strlowcase(pix_varnames) eq strlowcase(pixel_varname[pol_i]), count)

        if count eq 0 then begin
          message, pixel_varname[pol_i] + ' is not present in pixelfile (pixelfile=' + $
            pixelfile[pol_i, file_i] + ')'
        endif

        pix_size = getvar_savefile(pixelfile[pol_i, file_i], pixel_varname[pol_i], /return_size)
        if pix_size[0] gt 0 then pix_dims = pix_size[1:pix_size[0]]
        if total(abs(pix_dims - data_dims[0])) ne 0 then begin
          message, 'Number of Healpix pixels does not match data dimension'
        endif

        this_pix_vals = getvar_savefile(pixelfile[pol_i, file_i], pixel_varname[pol_i])
        if pol_i eq 0 and file_i eq 0 then begin
          pix_vals = this_pix_vals
        endif else if total(abs(this_pix_vals - pix_vals)) ne 0 then begin
          message, 'pixel values in files do not match'
        endif
      endif else healpix = 0

      wh_obs = where(strmatch(varnames, 'obs*',/fold_case) gt 0, count_obs)
      if count_obs eq 1 then obs_arr = getvar_savefile(datafile[pol_i, file_i], varnames[wh_obs[0]])
      if count_obs gt 1 then message, 'more than one obs structure in datafile'

      if count_obs ne 0 then begin
        n_obs[pol_i, file_i] = n_elements(obs_arr)

        if j eq 0 then begin
          max_baseline_lambda = max(obs_arr.max_baseline)
        endif else begin
          max_baseline_lambda = max([max_baseline_lambda, obs_arr.max_baseline])
        endelse

        obs_radec_vals = [[obs_arr.obsra],[obs_arr.obsdec]]
        zen_radec_vals = [[obs_arr.zenra],[obs_arr.zendec]]

        if total(abs(obs_arr.n_freq - obs_arr[0].n_freq)) ne 0 then begin
          message, 'inconsistent number of frequencies in obs_arr'
        endif
        if j eq 0 then n_freq = obs_arr[0].n_freq else if obs_arr[0].n_freq ne n_freq then begin
          message, 'n_freq does not agree between datafiles'
        endif

        if total(abs(obs_arr.degpix - obs_arr[0].degpix)) ne 0 then begin
          message, 'inconsistent degpix values in obs_arr'
        endif
        if j eq 0 then degpix = obs_arr[0].degpix else if obs_arr[0].degpix ne degpix  then begin
          message, 'degpix does not agree between datafiles'
        endif

        if total(abs(obs_arr.kpix - obs_arr[0].kpix)) ne 0 then begin
          message, 'inconsistent kpix values in obs_arr'
        endif
        if j eq 0 then kpix = obs_arr[0].kpix else if obs_arr[0].kpix ne kpix then begin
          message, 'kpix does not agree between datafiles'
        endif

        if total(abs(obs_arr.dimension - obs_arr[0].dimension)) ne 0 then begin
          message, 'inconsistent dimension values in obs_arr'
        endif
        if j eq 0 then fhd_dim = obs_arr[0].dimension else if obs_arr[0].dimension ne fhd_dim then begin
          message, 'dimension does not agree between datafiles'
        endif

        if total(abs(obs_arr.elements - obs_arr[0].elements)) ne 0 then begin
          message, 'inconsistent elements values in obs_arr'
        endif
        if j eq 0 then fhd_elem = obs_arr[0].elements else if obs_arr[0].elements ne fhd_elem then begin
          message, 'elements does not agree between datafiles'
        endif

        if fhd_dim ne fhd_elem then message, 'fhd image is not square in x & y'
        kspan = kpix * fhd_dim

        obs_tags = tag_names(obs_arr)
        wh_freq = where(strlowcase(obs_tags) eq 'freq', count_freq)
        if count_freq ne 0 then freq_vals = obs_arr.freq else begin
          wh_bin = where(strlowcase(obs_tags) eq 'bin', count_bin)
          wh_base_info = where(strlowcase(obs_tags) eq 'baseline_info', count_base_info)
          if count_bin gt 0 or count_base_info gt 0 then begin
            freq_vals = dblarr(n_freq, n_obs[pol_i, file_i])
            if count_bin ne 0 then begin
              for i=0, n_obs[pol_i, file_i]-1 do freq_vals[*,i] = (*obs_arr[i].bin).freq
            endif else begin
              for i=0, n_obs[pol_i, file_i]-1 do freq_vals[*,i] = (*obs_arr[i].baseline_info).freq
            endelse
          endif else stop
        endelse
        if total(abs(freq_vals - rebin(freq_vals[*,0], n_freq, n_obs[pol_i, file_i]))) ne 0 then begin
          message, 'inconsistent freq values in obs_arr'
        endif
        if j eq 0 then freq = freq_vals[*,0] else if total(abs(freq - freq_vals[*,0])) ne 0 then begin
          message, 'frequencies do not agree between datafiles'
        endif
        freq_resolution = freq[1]-freq[0]

        dt_vals = dblarr(n_obs[pol_i, file_i])
        wh_timeres = where(strlowcase(obs_tags) eq 'time_res', count_timeres)
        if count_timeres ne 0 then dt_vals = obs_arr.time_res else begin
          for i=0, n_obs[pol_i, file_i]-1 do begin
            times = (*obs_arr[i].baseline_info).jdate
            ;; only allow time resolutions of n*.5 sec
            dt_vals[i] = round(((times[1]-times[0])*24*3600)*2.)/2.
          endfor
          if total(abs(dt_vals - dt_vals[0])) ne 0 then begin
            message, 'inconsistent time averaging in obs_arr'
          endif
        endelse
        if j eq 0 then time_resolution = dt_vals[0] else begin
          if total(abs(time_resolution - dt_vals[0])) ne 0 then begin
            message, 'time averaging does not agree between datafiles'
          endif
        endelse
        theta_vals = angle_difference(obs_radec_vals[*,1], obs_radec_vals[*,0], $
          zen_radec_vals[*,1], zen_radec_vals[*,0], /degree, /nearest)
        if j eq 0 then begin
          max_theta = max(theta_vals)
        endif else begin
          max_theta = max([max_theta, theta_vals])
        endelse

        if j eq 0 then n_vis = dblarr(npol, nfiles)
        n_vis[pol_i, file_i] = total(obs_arr.n_vis,/double)

        if j eq 0 then n_vis_freq = dblarr(npol, nfiles, n_freq)
        n_vis_freq_arr = dblarr([n_obs[pol_i, file_i], n_freq])
        for i=0, n_obs[pol_i, file_i]-1 do n_vis_freq_arr[i, *] = obs_arr[i].nf_vis

        if tag_exist(obs_arr, 'beam_integral') or tag_exist(obs_arr, 'primary_beam_sq_area') then begin
          if j eq 0 then beam_int = fltarr(npol, nfiles, n_freq)
          beam_int_arr = fltarr([n_obs[pol_i, file_i], n_freq])
          if tag_exist(obs_arr,'primary_beam_sq_area') then begin
            for i=0, n_obs[pol_i, file_i]-1 do begin
              if obs_arr[i].primary_beam_sq_area(pol_i) eq !Null then begin
                beam_int_arr[i, *] = 0
              endif else begin
                beam_int_arr[i, *] = *(obs_arr[i].primary_beam_sq_area(pol_i))
              endelse
            endfor
          endif else begin
            for i=0, n_obs[pol_i, file_i]-1 do begin
              if obs_arr[i].beam_integral(pol_i) eq !Null then begin
                beam_int_arr[i, *] = 0
              endif else begin
                beam_int_arr[i, *] = *(obs_arr[i].beam_integral(pol_i))
              endelse
            endfor
          endelse
          if max(beam_int_arr) eq 0 then begin
            ;; all the pointers were null -- treat as if it doesn't exist in obs structure
            undefine, beam_int_arr
          endif else begin
            if min(beam_int_arr) eq 0 then begin
              message, 'some beam_integrals are null, others are not.'
            endif

            wh_freq_all_0 = where(total(n_vis_freq_arr, 1) eq 0, count_freq_all_0)
            beam_int[pol_i, file_i, *] = total(beam_int_arr * n_vis_freq_arr, 1) / $
              total(n_vis_freq_arr, 1)
            if count_freq_all_0 gt 0 then beam_int[pol_i, file_i, wh_freq_all_0]=0
          endelse
        endif

        if tag_exist(obs_arr[0], 'vis_noise') then begin
          vis_noise_arr = fltarr([n_obs[pol_i, file_i], n_freq])

          noise_dims = size(*obs_arr[0].vis_noise, /dimension)

          if noise_dims[0] lt npol or noise_dims[1] ne n_freq then begin
            if min(pol_exist) eq 1 and n_elements(*obs_arr[0].vis_noise) eq n_freq then begin
              for i=0, n_obs[pol_i, file_i]-1 do begin
                vis_noise_arr[i, *] = (*obs_arr[i].vis_noise)
              endfor
            endif else begin
              message, 'vis_noise dimensions do not match npol, n_freq'
            endelse
          endif else begin
            wh_pol = where(strmatch(obs_arr[0].pol_names, pol_inc[pol_i], /fold_case), count)
            if count ne 1 then begin
              message, 'pol can not be uniquely identified in obs_arr.pol_names'
            endif
            for i=0, n_obs[pol_i, file_i]-1 do vis_noise_arr[i, *] = (*obs_arr[i].vis_noise)[wh_pol,*]
          endelse

          if j eq 0 then vis_noise = fltarr(npol, nfiles, n_freq)

          vis_noise[pol_i, file_i, *] = sqrt(total(vis_noise_arr^2.*n_vis_freq_arr, 1)/total(n_vis_freq_arr, 1))
          wh_zero = where(total(n_vis_freq_arr, 1) eq 0, count_zero)
          if count_zero gt 0 then vis_noise[*, *, wh_zero] = 0
          undefine, wh_zero, count_zero
        endif

        n_vis_freq[pol_i, file_i, *] = total(n_vis_freq_arr, 1)

        undefine_fhd, obs_arr
      endif else begin
        message, 'no obs or obs_arr in datafile'
      endelse


      if healpix then data_nfreq = data_dims[1] else data_nfreq = data_dims[2]
      wh_navg = where(strlowcase(varnames) eq 'n_avg', count_navg)
      if count_navg ne 0 then begin
        if j eq 0 then n_avg = getvar_savefile(datafile[pol_i, file_i], 'n_avg') else begin
          n_avg1 = n_avg
          n_avg = getvar_savefile(datafile[pol_i, file_i], 'n_avg')
          if n_avg1 ne n_avg then begin
            message, 'n_avg parameter does not agree between datafiles'
          endif
          undefine, n_avg1
        endelse
      endif else begin
        print, 'no n_avg present, calculate from data frequency length; n_avg = ', $
          n_freq/data_nfreq
        if j eq 0 then begin
          n_avg = n_freq/data_nfreq
        endif else if n_avg ne n_freq/data_nfreq then begin
          message, 'calculated n_avg does not agree between datafiles'
        endif
      endelse
      if n_freq/n_avg ne data_nfreq then begin
        message, 'number of frequencies does not match number of data slices'
      endif

    endfor

    ;; fix uvf savefiles for gridded uv
    if healpix eq 0 and keyword_set(uvf_input) then begin
      uvf_savefile = datafile
      uvf_weight_savefile = weightfile
    endif

    n_freqbins = n_freq / n_avg
    frequencies = dblarr(n_freqbins)

    if n_elements(vis_noise) gt 0 then vis_noise_avg = fltarr(npol, nfiles, n_freqbins)
    if n_elements(beam_int) gt 0 then beam_int_avg = fltarr(npol, nfiles, n_freqbins)
    n_vis_freq_avg = dblarr(npol, nfiles, n_freqbins)
    if n_avg gt 1 then begin
      for i=0, n_freqbins-1 do begin
        frequencies[i] = mean(freq[i*n_avg:i*n_avg+(n_avg-1)]) / 1e6 ;; in MHz

        inds_arr = indgen(n_freq)
        if n_elements(vis_noise) gt 0 then begin
          inds_use = inds_arr[i*n_avg:i*n_avg+(n_avg-1)]
          wh_zero = where(total(total(n_vis_freq[*,*,inds_use],2),1) eq 0, $
            count_zero, complement = wh_gt0, ncomplement = count_gt0)

          if count_gt0 eq 0 then continue

          if count_zero gt 0 then inds_use = inds_use[wh_gt0]

          if count_gt0 eq 1 then begin
            vis_noise_avg[*,*,i] = vis_noise[*,*,inds_use]
            if n_elements(beam_int) gt 0 then beam_int_avg[*,*,i] = beam_int[*,*,inds_use]
            n_vis_freq_avg[*,*,i] = n_vis_freq[*,*,inds_use]
          endif else begin
            vis_noise_avg[*,*,i] = sqrt(total(vis_noise[*,*,inds_use]^2.*n_vis_freq[*,*,inds_use], 3) / $
              total(n_vis_freq[*,*,inds_use],3))
            if n_elements(beam_int) gt 0 then begin
              beam_int_avg[*,*,i] = total(beam_int[*,*,inds_use]*n_vis_freq[*,*,inds_use], 3) / $
                total(n_vis_freq[*,*,inds_use],3)
            endif
            n_vis_freq_avg[*,*,i] = total(n_vis_freq[*,*,inds_use], 3)
          endelse

        endif
      endfor
      if n_elements(vis_noise) gt 0 then vis_noise = vis_noise_avg
      if n_elements(beam_int) gt 0 then beam_int = beam_int_avg

      n_vis_freq = n_vis_freq_avg
    endif else begin
      frequencies = freq / 1e6 ;; in MHz
    endelse

    ;; check frequency spacing to determine if a frequency DFT is required.
    if not ps_options.freq_dft then begin
      freq_diff = frequencies - shift(frequencies, 1)
      freq_diff = freq_diff[1:*]
      if max(abs(freq_diff-freq_diff[0])) gt 1e-12 then begin
        ps_options = create_ps_options(ps_options = ps_options, $
          freq_dft = 1)
      endif
    endif

    metadata_struct = {datafile: datafile, weightfile: weightfile, $
      variancefile:variancefile, cube_varname:cube_varname, $
      weight_varname:weight_varname, variance_varname:variance_varname, $
      frequencies:frequencies, freq_resolution:freq_resolution, $
      time_resolution:time_resolution, $
      n_vis:n_vis, n_vis_freq:n_vis_freq, max_baseline_lambda:max_baseline_lambda, $
      max_theta:max_theta, degpix:degpix, kpix:kpix, kspan:kspan, $
      general_filebase:general_filebase, infile_label:infile_label, $
      type_pol_str:type_pol_str, type_inc:type_inc, n_obs:n_obs, pol_inc:pol_inc, $
      nfiles:nfiles}

    if n_elements(beamfile) gt 0 then begin
      metadata_struct = create_struct(metadata_struct, 'beamfile', beamfile, $
        'beam_varname', beam_varname)
    endif

    if healpix then begin
      metadata_struct = create_struct(metadata_struct, 'pixelfile', pixelfile, $
        'pixel_varname', pixel_varname, 'nside', nside)
    endif

    if no_var then metadata_struct = create_struct(metadata_struct, 'no_var', 1)

    if n_elements(vis_noise) gt 0 then begin
      metadata_struct = create_struct(metadata_struct, 'vis_noise', vis_noise)
    endif

    if n_elements(beam_int) gt 0 then begin
      metadata_struct = create_struct(metadata_struct, 'beam_int', beam_int)
    endif

    save, filename = info_file, metadata_struct
  endif

  npol = n_elements(metadata_struct.pol_inc)
  nfiles = metadata_struct.nfiles
  ntypes = n_elements(metadata_struct.type_inc)
  if tag_exist(metadata_struct, 'nside') then healpix = 1 else healpix = 0
  if tag_exist(metadata_struct, 'no_var') then no_var = 1 else no_var = 0

  pol_exist = stregex(metadata_struct.datafile, '[xy][xy]', /boolean, /fold_case)

  if n_elements(freq_flags) ne 0 then begin
    freq_mask = intarr(n_elements(metadata_struct.frequencies)) + 1
    freq_mask[freq_flags] = 0

    if n_elements(freq_ch_range) gt 0 then begin
      freq_mask = freq_mask[min(freq_ch_range):max(freq_ch_range)]
    endif
    wh_freq_use = where(freq_mask gt 0, count_freq_use)
    if count_freq_use lt 3 then message, 'Too many frequencies flagged'

  endif

  if n_elements(freq_flags) ne 0 then begin
    if min(freq_flags) lt 0 or max(freq_flags) gt n_elements(metadata_struct.frequencies) then begin
      message, 'invalid freq_flags'
    endif
  endif

  file_tags = create_file_tags(freq_ch_range = freq_ch_range, freq_flags = freq_flags, $
    freq_flag_name = freq_flag_name, uvf_options = uvf_options, ps_options = ps_options)

  wt_file_label = '_weights_' + strlowcase(metadata_struct.pol_inc)
  file_label = strarr(npol, ntypes)
  for i=0, npol-1 do begin
    file_label[i,*] = '_' + strlowcase(metadata_struct.type_inc) + '_' + $
      strlowcase(metadata_struct.pol_inc[i])
  endfor
  savefilebase = metadata_struct.general_filebase + file_tags.uvf_tag + file_label

  uvf_savefilebase = strarr(npol, nfiles, ntypes)
  uvf_label = strarr(npol, nfiles, ntypes)
  for pol_i=0, npol-1 do begin
    for file_i=0, nfiles-1 do begin
      if max(pol_exist) eq 1 then begin
        uvf_savefilebase[pol_i, file_i,*] = cgRootName(metadata_struct.datafile[pol_i, file_i]) + $
          file_tags.uvf_tag + '_' + metadata_struct.type_inc
        endif else begin
          uvf_savefilebase[pol_i, file_i,*] = cgRootName(metadata_struct.datafile[pol_i, file_i]) + $
            file_tags.uvf_tag + file_label[pol_i, *]
        endelse
      uvf_label[pol_i, file_i,*] = metadata_struct.infile_label[file_i] + '_' + $
        file_label[pol_i, *]
    endfor
  endfor
  if n_elements(size(uvf_savefilebase, /dim)) ne 3 then begin
    uvf_savefilebase = reform(uvf_savefilebase, [npol, nfiles, ntypes])
  endif

  ;; add sw tag to general_filebase so that plotfiles havefile_tags.uvf_tag in them
  general_filebase = metadata_struct.general_filebase + file_tags.uvf_tag

  subfolders = {data: 'data' + path_sep(), plots: 'plots' + path_sep(), $
    uvf: 'uvf_cubes' + path_sep(), beams: 'beam_cubes' + path_sep(), $
    kspace: 'kspace_cubes' + path_sep(), slices: 'slices' + path_sep(), $
    bin_2d: '2d_binning' + path_sep(), bin_1d: '1d_binning' + path_sep()}

  uvf_path = froot + subfolders.data + subfolders.uvf
  uvf_slice_path = uvf_path + subfolders.slices

  ;; file for saving ra/decs of pixels & their DFT values
  radec_file = uvf_path + metadata_struct.general_filebase + $
    '_radec.idlsave'

  uvf_savefile = uvf_path + uvf_savefilebase + '_uvf.idlsave'
  uf_savefile = uvf_slice_path + uvf_savefilebase + '_uf_plane.idlsave'
  vf_savefile = uvf_slice_path + uvf_savefilebase + '_vf_plane.idlsave'
  uv_savefile = uvf_slice_path + uvf_savefilebase + '_uv_plane.idlsave'

  uf_raw_savefile = uvf_slice_path + uvf_savefilebase + '_uf_plane_raw.idlsave'
  vf_raw_savefile = uvf_slice_path + uvf_savefilebase + '_vf_plane_raw.idlsave'
  uv_raw_savefile = uvf_slice_path + uvf_savefilebase + '_uv_plane_raw.idlsave'

  uf_sum_savefile = uvf_slice_path + savefilebase + '_sum_uf_plane.idlsave'
  vf_sum_savefile = uvf_slice_path + savefilebase + '_sum_vf_plane.idlsave'
  uv_sum_savefile = uvf_slice_path + savefilebase + '_sum_uv_plane.idlsave'
  uf_diff_savefile = uvf_slice_path + savefilebase + '_diff_uf_plane.idlsave'
  vf_diff_savefile = uvf_slice_path + savefilebase + '_diff_vf_plane.idlsave'
  uv_diff_savefile = uvf_slice_path + savefilebase + '_diff_uv_plane.idlsave'

  if tag_exist(metadata_struct, 'beamfile') then begin
    beam_savefile = froot + subfolders.data + subfolders.beams + uvf_savefilebase + '_beam2.idlsave'
  endif

  kspace_path = froot + subfolders.data + subfolders.kspace
  kspace_slice_path = kspace_path + subfolders.slices

  xz_savefile = kspace_slice_path + savefilebase + file_tags.power_tag + '_xz_plane.idlsave'
  yz_savefile = kspace_slice_path + savefilebase + file_tags.power_tag + '_yz_plane.idlsave'
  xy_savefile = kspace_slice_path + savefilebase + file_tags.power_tag + '_xy_plane.idlsave'

  kcube_savefile = kspace_path + savefilebase + file_tags.kcube_tag + '_kcube.idlsave'
  power_savefile = kspace_path + savefilebase + file_tags.power_tag + '_power.idlsave'
  fits_power_savefile = kspace_path + savefilebase + file_tags.power_tag + '_power.fits'

  weight_savefilebase = strarr(npol, nfiles)
  for pol_i=0, npol-1 do begin
    for file_i=0, nfiles-1 do begin
      if max(pol_exist) eq 1 then begin
        weight_savefilebase[pol_i, file_i] = cgRootName(metadata_struct.weightfile[pol_i, file_i]) + $
          file_tags.uvf_tag + '_weights'
        endif else begin
          weight_savefilebase[pol_i, file_i] = cgRootName(metadata_struct.weightfile[pol_i, file_i]) + $
            file_tags.uvf_tag + wt_file_label[pol_i]
        endelse
    endfor
  endfor

  var_savefilebase = strarr(npol, nfiles)
  for i=0, n_elements(weight_savefilebase)-1 do begin
    var_savefilebase[i] = strjoin(strsplit(weight_savefilebase[i], 'weights', /regex, /extract, $
      /preserve_null), 'variance')
  endfor

  uvf_weight_savefile = uvf_path + weight_savefilebase + '_uvf.idlsave'
  uf_weight_savefile = uvf_slice_path + weight_savefilebase + '_uf_plane.idlsave'
  vf_weight_savefile = uvf_slice_path + weight_savefilebase + '_vf_plane.idlsave'
  uv_weight_savefile = uvf_slice_path + weight_savefilebase + '_uv_plane.idlsave'

  uf_var_savefile = uvf_slice_path + var_savefilebase + '_uf_plane.idlsave'
  vf_var_savefile = uvf_slice_path + var_savefilebase + '_vf_plane.idlsave'
  uv_var_savefile = uvf_slice_path + var_savefilebase + '_uv_plane.idlsave'

  for pol_i=0, npol-1 do begin
    for type_i=0, ntypes-1 do begin

      if ntypes gt 1 then begin
        data_varname = metadata_struct.cube_varname[type_i, pol_i]
      endif else begin
        data_varname = metadata_struct.cube_varname[pol_i]
      endelse
      if data_varname ne '' then begin
        derived_uvf_inputfiles = strarr(nfiles,2)
        derived_uvf_varname = strarr(nfiles,2)
      endif else begin
        if healpix or not keyword_set(uvf_input) then begin
          derived_uvf_inputfiles = uvf_savefile[pol_i, *, 0:1]
          derived_uvf_varname = strarr(nfiles, 2) + 'data_cube'
        endif else begin
          derived_uvf_inputfiles = strarr(nfiles,2)
          derived_uvf_varname = strarr(nfiles,2)
          for j=0, nfiles-1 do begin
            derived_uvf_inputfiles[j,*] = metadata_struct.datafile[pol_i, j]
            derived_uvf_varname[j,*] = [metadata_struct.cube_varname[0, pol_i], $
              metadata_struct.cube_varname[1, pol_i]]
          endfor
        endelse
      endelse

      file_struct = {datafile:reform(metadata_struct.datafile[pol_i,*]), $
        weightfile:reform(metadata_struct.weightfile[pol_i,*]), $
        variancefile:reform(metadata_struct.variancefile[pol_i,*]), $
        datavar:data_varname, weightvar:metadata_struct.weight_varname[pol_i], $
        variancevar:metadata_struct.variance_varname[pol_i], $
        frequencies:metadata_struct.frequencies, $
        freq_resolution:metadata_struct.freq_resolution, $
        time_resolution:metadata_struct.time_resolution, $
        n_obs:reform(metadata_struct.n_obs[pol_i,*]), $
        n_vis:reform(metadata_struct.n_vis[pol_i,*]), $
        n_vis_freq:reform(metadata_struct.n_vis_freq[pol_i,*,*]), $
        max_baseline_lambda:metadata_struct.max_baseline_lambda, $
        max_theta:metadata_struct.max_theta, degpix:metadata_struct.degpix, $
        kpix:metadata_struct.kpix, kspan:metadata_struct.kspan, $
        uf_savefile:reform(uf_savefile[pol_i,*,type_i]), $
        vf_savefile:reform(vf_savefile[pol_i,*,type_i]), $
        uv_savefile:reform(uv_savefile[pol_i,*,type_i]), $
        uf_raw_savefile:reform(uf_raw_savefile[pol_i,*,type_i]), $
        vf_raw_savefile:reform(vf_raw_savefile[pol_i,*,type_i]), $
        uv_raw_savefile:reform(uv_raw_savefile[pol_i,*,type_i]), $
        uf_sum_savefile:uf_sum_savefile[pol_i, type_i], $
        vf_sum_savefile:vf_sum_savefile[pol_i,type_i], $
        uv_sum_savefile:uv_sum_savefile[pol_i,type_i], $
        uf_diff_savefile:uf_diff_savefile[pol_i,type_i], $
        vf_diff_savefile:vf_diff_savefile[pol_i,type_i], $
        uv_diff_savefile:uv_diff_savefile[pol_i,type_i], $
        uf_weight_savefile:reform(uf_weight_savefile[pol_i, *]), $
        vf_weight_savefile:reform(vf_weight_savefile[pol_i, *]), $
        uv_weight_savefile:reform(uv_weight_savefile[pol_i, *]), $
        uf_var_savefile:reform(uf_var_savefile[pol_i, *]), $
        vf_var_savefile:reform(vf_var_savefile[pol_i, *]), $
        uv_var_savefile:reform(uv_var_savefile[pol_i, *]), $
        xz_savefile:xz_savefile[pol_i,type_i], $
        yz_savefile:yz_savefile[pol_i,type_i], $
        xy_savefile:xy_savefile[pol_i,type_i], $
        kcube_savefile:kcube_savefile[pol_i,type_i], $
        power_savefile:power_savefile[pol_i,type_i], $
        fits_power_savefile:fits_power_savefile[pol_i,type_i],$
        savefile_froot:froot, subfolders:subfolders, $
        savefilebase:savefilebase[pol_i,type_i], $
        general_filebase:general_filebase, $
        weight_savefilebase:reform(weight_savefilebase[pol_i, *]), $
        derived_uvf_inputfiles:derived_uvf_inputfiles, $
        derived_uvf_varname:derived_uvf_varname, $
        file_label:file_label[pol_i,type_i], $
        uvf_label:reform(uvf_label[pol_i,*,type_i]), $
        wt_file_label:wt_file_label[pol_i], $
        uvf_tag:file_tags.uvf_tag, kcube_tag:file_tags.kcube_tag, $
        power_tag:file_tags.power_tag, $
        type_pol_str:metadata_struct.type_pol_str[pol_i,type_i], $
        pol_index:pol_i, type_index:type_i, pol:metadata_struct.pol_inc[pol_i], $
        type:metadata_struct.type_inc[type_i], nfiles:nfiles}

      if tag_exist(metadata_struct, 'beamfile') then begin
        file_struct = create_struct(file_struct, 'beamfile', reform(metadata_struct.beamfile[pol_i,*]), $
          'beamvar', metadata_struct.beam_varname[pol_i], 'beam_savefile', reform(beam_savefile[pol_i, *]))
      endif

      if healpix or not keyword_set(uvf_input) then begin
        file_struct = create_struct(file_struct, 'radec_file', radec_file)

        file_struct = create_struct(file_struct, 'uvf_savefile', reform(uvf_savefile[pol_i,*,type_i]), $
          'uvf_weight_savefile', uvf_weight_savefile[pol_i,*])
      endif

      if healpix then file_struct = create_struct(file_struct, 'pixelfile', $
        metadata_struct.pixelfile[pol_i,*], 'pixelvar', $
        metadata_struct.pixel_varname[pol_i,*], 'nside', metadata_struct.nside)

      if no_var then file_struct = create_struct(file_struct, 'no_var', 1)

      if n_elements(freq_mask) gt 0 then begin
        file_struct = create_struct(file_struct, 'freq_mask', freq_mask)
      endif

      if tag_exist(metadata_struct, 'vis_noise') gt 0 then begin
        file_struct = create_struct(file_struct, 'vis_noise', reform(metadata_struct.vis_noise[pol_i,*,*]))
      endif

      if tag_exist(metadata_struct, 'beam_int') then begin
        file_struct = create_struct(file_struct, 'beam_int', reform(metadata_struct.beam_int[pol_i,*,*]))
      endif

      if pol_i eq 0 and type_i eq 0 then begin
        file_struct_arr = replicate(file_struct, ntypes, npol)
      endif else begin
        file_struct_arr[type_i, pol_i] = file_struct
      endelse
    endfor
  endfor

  return, reform(file_struct_arr, ntypes*npol)

end
