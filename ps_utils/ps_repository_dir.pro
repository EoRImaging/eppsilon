function ps_repository_dir

  case strlowcase(!version.os_family) OF
    'windows': split_delim = ';'
    'unix':    split_delim = ':'
  endcase
  
  path_dirs = strsplit(!path, split_delim, /extract)
  
  ps_cor_loc = strpos(path_dirs, 'ps_core')
  wh_pscore = where(ps_cor_loc gt 0, count_pscore)
  
  if count_pscore eq 0 then message, 'ps_core cannot be located in !path variable.'
  
  ps_dir = file_dirname(path_dirs[wh_pscore[0]])
  
  return, ps_dir
  
end