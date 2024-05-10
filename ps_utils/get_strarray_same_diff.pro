function get_strarray_same_diff, strarray, delimiter, same_str=same_str, $
  diff_count=diff_count, match_arr=match_arr

  compile_opt idl2, logical_predicate

  n_str = n_elements(strarray)
  parts_list = strsplit(strarray, delimiter, /extract, $
    count=parts_count)

  match_arr = intarr(n_str-1, max(parts_count))
  for fi=1, n_str-1 do begin
    this_match_test = strcmp(parts_list[0], parts_list[fi])
    this_same = where(this_match_test eq 1)
    match_arr[fi-1, this_same] = 1
  endfor

  same_inds = where(total(match_arr, 1) eq (n_str-1), $
    count_same, complement=diff_inds, ncomplement=diff_count)

  if count_same gt 0 then begin
    same_str = strjoin((parts_list[0])[same_inds], delimiter)
  endif else begin
    same_str = ''
  endelse

  diff_parts = strarr(n_str)
  if diff_count gt 0 then begin
    for fi=0, n_str-1 do begin
      diff_parts[fi] = strjoin((parts_list[fi])[diff_inds], delimiter)
    endfor
  endif

  return, diff_parts
end
