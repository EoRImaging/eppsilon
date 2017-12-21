function get_folder, folder_name_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder

  if n_elements(loc_name) eq 0 then begin
    spawn, 'hostname', hostname
    if stregex(hostname, 'mit.edu', /boolean) eq 1 then loc_name = 'mit'
    if stregex(hostname, 'enterprise', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'constellation', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'defiant', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'intrepid', /boolean) eq 1 then loc_name = 'enterprise'
    if n_elements(loc_name) eq 0 then loc_name = hostname
  endif
  case loc_name of
    'mit':  folder_name = mit_folder_locs(folder_name_in, rts = rts)
    'enterprise': folder_name = enterprise_folder_locs(folder_name_in, rts = rts)
    else: begin
      folder_test = file_test(folder_name_in, /directory)
      if min(folder_test) eq 0 then message, 'machine not recognized and folder not found'
      folder_name = folder_name_in
    endelse
  endcase

  if keyword_set(rts) and n_elements(dirty_folder) gt 0 then begin
    case loc_name of
      'mit':  dirty_folder = mit_folder_locs(dirty_folder, rts = rts)
      'enterprise': dirty_folder = enterprise_folder_locs(dirty_folder, rts = rts)
      else: begin
        folder_test = file_test(dirty_folder, /directory)
        if min(folder_test) eq 0 then message, 'machine not recognized and dirty_folder not found'
      endelse
    endcase
  endif

return, folder_name
end
