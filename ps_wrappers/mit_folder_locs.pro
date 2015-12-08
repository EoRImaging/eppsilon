function mit_folder_locs, folder_names_in, rts = rts

  folder_names = folder_names_in
  if keyword_set(rts) then begin
    if n_elements(folder_name) eq 0 then folder_name = '/data3/MWA/bpindor/RTS/dec_11/'
    
  endif else begin
  
    for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'
      start_path = '/nfs/mwa-09/r1/djc/'
      folder_test = file_test(folder_names[i], /directory)
      if folder_test eq 0 then begin
        pos_eor2013 = strpos(folder_names[i], 'EoR2013')
        if pos_eor2013 gt -1 then begin
          test_name = start_path + strmid(folder_names[i], pos_eor2013)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        pos_aug23 = strpos(folder_names[i], 'Aug23')
        if pos_aug23 gt -1 then begin
          test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_aug23)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        pos_aug26 = strpos(folder_names[i], 'Aug26')
        if pos_aug26 gt -1 then begin
          test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_aug26)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        pos_week1 = strpos(folder_names[i], 'week1')
        if pos_week1 gt -1 then begin
          test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_week1)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        test_name = start_path + 'EoR2013/Aug23/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      if folder_test eq 0 then begin
        test_name = start_path + 'EoR2013/Aug26/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      if folder_test eq 0 then begin
        test_name = start_path + 'EoR2013/week1/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      
      
      start_path2 = '/nfs/mwa-03/r1/'
      if folder_test eq 0 then begin
        pos_eor2013 = strpos(folder_names[i], 'EoR2013')
        if pos_eor2013 gt -1 then begin
          test_name = start_path2 + strmid(folder_names[i], pos_eor2013)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        test_name = start_path2 + 'EoR2013/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      
      
      if folder_test eq 0 then message, 'folder not found'
    endfor
  endelse
  
  return, folder_names
  
end