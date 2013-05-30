;; set the periodic keyword for cases in which n_samples is even but
;; there is a single maximum (ie k=0 exists) -- this implies an asymmetric window
;; function and is correct for FFTs with even numbers of samples

function spectral_window, n_samples, type = type, periodic = periodic

  type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']

  if n_elements(type) eq 0 then type = 'Blackman-Harris'
  wh_type = where(strlowcase(type_list) eq strlowcase(type), count_type)
  if count_type eq 0 then message, 'Spectral window type not recognized.' else type = type_list[wh_type[0]]

  if n_samples lt 2 then message, 'n_samples must be greater than 2'

  if n_samples mod 2 eq 0 and keyword_set(periodic) then n_use = n_samples+1 else n_use = n_samples
  
  cos_term1 = cos(2.*!pi*findgen(n_use)/(n_use-1))
  cos_term2 = cos(4.*!pi*findgen(n_use)/(n_use-1))
  cos_term3 = cos(6.*!pi*findgen(n_use)/(n_use-1))
  
  case type of
     'Hann': window = 0.5 * (1-cos_term1)
     'Hamming': window = 0.54 - 0.46 * cos_term1
     'Blackman': window = (1-0.16)/2. - 0.5 * cos_term1 + (0.16/2.) * cos_term2
     'Nutall': window = 0.355768 - 0.487396 * cos_term1 + 0.144232 * cos_term2 - 0.012604 * cos_term3
     'Blackman-Nutall': window = 0.3635819 - 0.4891775 * cos_term1 + 0.1365995 * cos_term2 - 0.0106411 * cos_term3
     'Blackman-Harris': window = 0.35875 - 0.48829 * cos_term1 + 0.14128 * cos_term2 - 0.01168 * cos_term3
  endcase

  return, window[0:n_samples-1]

end
