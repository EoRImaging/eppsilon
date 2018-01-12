function avg_array, array, factor = factor, axis = axis

  dims = size(array, /dimension)
  new_dims = dims
  new_dims[axis-1] = floor(new_dims[axis-1] / factor)

  temp = complex(fltarr(new_dims))

end
