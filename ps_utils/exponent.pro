function exponent, axis, index, number

  ;; A special case.
  if number eq 0 then return, '0' 

  return, number_formatter(number, format='(e10.0)',/print_exp)

end
