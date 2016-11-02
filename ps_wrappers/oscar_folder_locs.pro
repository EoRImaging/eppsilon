function oscar_folder_locs, folder_names_in

  folder_names = folder_names_in
  for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path.'
;      start_path = '/gpfs/scratch/alanman'
      start_path = '/gpfs/data/jpober/alanman'
      folder_test = file_test(folder_names[i], /directory)
      if folder_test eq 0 then begin
        pos_fhdout = strpos(folder_names[i], 'FHD_out')
        if pos_fhdout gt -1 then begin
          test_name = start_path + strmid(folder_names[i], pos_fhdout)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then message, 'folder not found'
  endfor
  
  return, folder_names
  
end
