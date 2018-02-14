function check_old_path, file_struct, tag_name, index = index

  tag_loc = where(strlowcase(tag_names(file_struct)) eq strlowcase(tag_name), count_loc)
  if count_loc eq 0 then begin
    message, 'tag_name ' + tag_name + ' not found in file_struct'
  endif
  tag_loc = tag_loc[0]

  if n_elements(index) gt 0 then begin
    old_file = file_struct.savefile_froot + file_basename(file_struct.(tag_loc)[index])
  endif else begin
    old_file = file_struct.savefile_froot + file_basename(file_struct.(tag_loc))
  endelse

  test_uvf_old = file_valid(old_file)

  if test_uvf_old then begin
    if n_elements(index) gt 0 then begin
      file_struct.(tag_loc)[index] = old_file
    endif else begin
      file_struct.(tag_loc) = old_file
    endelse
  endif

return, test_uvf_old
end
