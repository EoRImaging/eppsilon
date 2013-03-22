function exponent, axis, index, number

  ;; A special case.
  if number eq 0 then return, '0' 

  return, number_formatter(number, format='(e13.1)',/print_exp)

end
