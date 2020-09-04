function get_cube_uvf_input, savefile, varname, n_freq, pol_index

  cube = getvar_savefile(savefile, varname)
  if size(cube, /n_dim) eq 0 then begin
    message, string(savefile) + ' not found. Verify path on filesystem.'
  endif

  cube_size = size(cube)
  if cube_size[n_elements(cube_size)-2] eq 10 then begin
    ;; cube is a pointer
    dims2 = size(*cube[0], /dimension)
    iter=1
    while min(dims2) eq 0 and iter lt cube_size[2] do begin
      print, 'warning: some frequency slices have null pointers'
      dims2 = size(*cube[iter], /dimension)
      iter = iter+1
    endwhile

    temp = complex(fltarr([dims2, n_freq]))
    for i = 0, n_freq-1 do begin
      foo = *cube[pol_index, i]
      if foo ne !null then begin
        temp[*,*,i] = temporary(foo)
      endif
    endfor
    undefine_fhd, cube

    cube = temporary(temp)

  endif

  return, cube
end
