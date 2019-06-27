function get_folder, folder_name_in, loc_name = loc_name,  rts = rts, $
    dirty_folder = dirty_folder

  if n_elements(loc_name) eq 0 then begin
    spawn, 'hostname', hostname
    if stregex(hostname, 'enterprise', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'constellation', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'defiant', /boolean) eq 1 then loc_name = 'enterprise'
    if stregex(hostname, 'intrepid', /boolean) eq 1 then loc_name = 'enterprise'
    if n_elements(loc_name) eq 0 then loc_name = hostname
  endif
  case loc_name of
    'enterprise': folder_name = enterprise_folder_locs(folder_name_in, rts = rts)
    else: begin
      folder_test = file_test(folder_name_in, /directory)
      if min(folder_test) eq 0 then message, 'machine not recognized and folder not found'
      folder_name = folder_name_in
    endelse
  endcase

  if keyword_set(rts) and n_elements(dirty_folder) gt 0 then begin
    case loc_name of
      'enterprise': dirty_folder = enterprise_folder_locs(dirty_folder, rts = rts)
      else: begin
        folder_test = file_test(dirty_folder, /directory)
        if min(folder_test) eq 0 then message, 'machine not recognized and dirty_folder not found'
      endelse
    endcase
  endif

;; make sure there are no path separators at the end
if n_elements(folder_name) gt 1 then begin
  for i=0, n_elements(folder_name) do begin
    folder_name[i] = expand_path(folder_name[i])
  endfor
endif else begin
  folder_name = expand_path(folder_name)
endelse

return, folder_name
end
