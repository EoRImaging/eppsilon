function enterprise_folder_locs, folder_names_in, rts = rts

  folder_names = folder_names_in
  
  if keyword_set(rts) then begin
  
    ;; check for folder existence, otherwise look for common folder names to figure out full path.
    ;; If none found, try '/data3/MWA/bpindor/RTS/'
    start_path = '/data3/MWA/bpindor/RTS/'
    folder_test = file_test(folder_names, /directory)
    if folder_test eq 0 then begin
      pos_RTS = strpos(folder_names, 'RTS')
      if pos_RTS gt -1 then begin
        test_name = start_path + strmid(folder_names, pos_RTS)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names = test_name
      endif
    endif
    if folder_test eq 0 then begin
      test_name = start_path + folder_names
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_names = test_name
    endif
    
    if folder_test eq 0 then message, 'folder not found'
    
    
  endif else begin
  
    for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path.
      start_path = '/data4/MWA/'
      folder_test = file_test(folder_names[i], /directory)
      if folder_test eq 0 then begin
        pos_aug23 = strpos(folder_names, 'FHD_Aug23')
        if pos_aug23 gt -1 then begin
          test_name = start_path + strmid(folder_names[i], pos_aug23)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        test_name = start_path + 'FHD_Aug23/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      
      if folder_test eq 0 then message, 'folder not found'
    endfor
  endelse
  
  return, folder_names
end