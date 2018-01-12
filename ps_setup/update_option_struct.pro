function update_option_struct, struct, tags, values, remove_tags = remove_tags, $
    return_new = return_new

  if n_elements(tags) ne n_elements(values) then begin
    message, 'number of tags must equal number of values'
  endif

  existing_tags = strlowcase(tag_names(struct))

    if n_elements(remove_tags) gt 0 then begin
      for i=0, n_elements(tags) - 1 do begin
        if max(strcmp(remove_tags, tags[i])) gt 0 then begin
          message, 'tag name "' + tags[i] + '" is in tags to modify and in ' + $
            'remove_tags list, which is not allowed'
        endif
      endfor

      for i=0, n_elements(remove_tags) - 1 do begin
        if not tag_exist(struct, remove_tags[i]) then begin
          message, 'tag to remove "' + remove_tags[i] + '" is not in struct'
        endif
      endfor
    endif

  if n_elements(remove_tags) gt 0 or keyword_set(return_new) then begin

    tags_to_keep = list()
    values_to_keep = list()
    for i=0, n_elements(existing_tags)-1 do begin
      if n_elements(remove_tags) gt 0 then begin
        if max(strcmp(remove_tags, existing_tags[i])) eq 0 then keep = 1 else keep = 0
      endif else keep = 1

      if keep then begin
        tags_to_keep.add, existing_tags[i]
        if n_elements(tags) gt 0 then begin
          if max(strcmp(tags, existing_tags[i])) eq 1 then begin
            tindex = (where(strcmp(tags, existing_tags[i]) EQ 1))[0]
            values_to_keep.add, values[tindex]
          endif else begin
            values_to_keep.add, struct.(i)
          endelse
        endif else begin
          values_to_keep.add, struct.(i)
        endelse
      endif
    endfor

    if n_elements(tags_to_keep) lt 1 then message, 'all tags are marked for deletion'
    tags_to_keep = tags_to_keep.toarray()
    new_struct = create_struct(tags_to_keep[0], values_to_keep[0])
    for i=1, n_elements(tags_to_keep) - 1 do begin
      new_struct = create_struct(new_struct, tags_to_keep[i], values_to_keep[i])
    endfor
    return, new_struct
  endif else begin
    for i=0, n_elements(tags) - 1 do begin
      if tag_exist(struct, tags[i]) then begin
        tindex = (where(strcmp(existing_tags, tags[i]) EQ 1))[0]
        struct.(tindex) = values[i]
      endif else begin
        struct = create_struct(struct, tags[i], values[i])
      endelse
    endfor

    return, struct
  endelse
end
