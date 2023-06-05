pro plot_psf, frequencies, plot_path, freq_name=freq_name, dft_z_use=dft_z_use

    ; fhd high band frequencies at 40 kHz cadence
    ; freqs = findgen(768,increment=0.04, start=167.04)
    
    ; create file name
    file_name = "dftz_" + dft_z_use
    if n_elements(freq_name) gt 0 then begin
        file_name += "_freq_" + freq_name
    endif else begin
        file_name += "_freq_standard"
    endelse
   
    save_path = plot_path + file_name
    save_path += ".png"
    print, file_name 
    print, save_path
    ; compute kz values
    n_freq = n_elements(frequencies) ;Selected number of frequencies

    z_mpc_mean = z_mpc(frequencies, hubble_param = hubble_param, f_delta = f_delta, $
    redshifts = redshifts, comov_dist_los = comov_dist_los, even_freq = even_freq, $
    z_mpc_delta = z_mpc_delta, z_mpc_length = z_mpc_length)

    kperp_lambda_conv = z_mpc_mean / (2.*!dpi)
    delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
    delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
    delay_params = [delay_delta, delay_max]

    ; z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
    kz_mpc_range =  (2.*!dpi) / (z_mpc_delta)
    kz_mpc_delta = (2.*!dpi) / z_mpc_length

    ;; need to figure out where the zero bin is (for later computations)
    ;; The location of the zero bin also dictates what the most negative kz bin is (min_kz)
    ;; which we use to construct the kz_mpc_orig array
    ;; The calculation depends a bit on whether we have an odd or even number of frequencies.
    ;; The following calculation is taken directly from the IDL FFT documentation
    int_arr = findgen((n_freq - 1)/2) + 1
    is_n_freq_even = (n_freq mod 2) eq 0
    fft_shift_val = -(n_freq/2 + 1)
    if (is_n_freq_even) then begin
        kz_integers = [0.0, int_arr, n_freq/2, -n_freq/2 + int_arr]
    endif else begin
        kz_integers = [0.0, int_arr, fft_shift_val + int_arr]
    endelse
    kz_integers = shift(kz_integers, fft_shift_val)
    if where(kz_integers eq min(kz_integers)) ne 0 then begin
        message, 'something went very wrong with shifting!'
    endif

    min_kz = min(kz_integers)*kz_mpc_delta
    kz_mpc_orig = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta + min_kz

    n_kz = n_elements(kz_mpc_orig)

    n_val = kz_integers
    wh_pos = where(n_val gt 0, n_pos)
    wh_neg = where(n_val lt 0, n_neg)

    ;; this should be zero by construction, but add a test just in case
    if kz_mpc_orig[where(n_val eq 0)] ne 0 then begin
        message, 'something has gone terribly wrong with calculating the kz_mpc_orig values.'
    endif

    kz_mpc = kz_mpc_orig[where(n_val ge 0)]

    if is_n_freq_even then begin
        ;; if n_freq is even, there's an extra positive mode (because there's always a zero mode.)
        if max(abs(n_val(wh_pos))) le max(abs(n_val(wh_neg))) then begin
        message, 'something has gone terribly wrong with calculating the kz_integers.'
        endif
        ;; remove unmatched max k positive mode:
        trim_max_pos = 1
        n_val = n_val[0:-2]
        kz_mpc = kz_mpc[0:-2]
        kz_mpc_orig_trim = kz_mpc_orig[0:-2]
        wh_pos = wh_pos[0:-2]
        n_pos = n_elements(wh_pos)
    endif else begin
        kz_mpc_orig_trim = kz_mpc_orig
        trim_max_pos = 0
    endelse

    ;; set up a hyperfine kz array
    ;; typically we would just use kz_mpc_orig_trim
    kz_fine_len = ((n_elements(kz_mpc_orig_trim)-1) * 100) + 1
    kz_fine_delta = kz_mpc_delta / 100
    kz_fine = findgen(kz_fine_len,increment=kz_fine_delta,start=kz_mpc_orig_trim[0])

    z_reg = (frequencies-frequencies[0])*z_mpc_delta/f_delta
    z_exp_zreg = exp(-1.*dcomplex(0,1)*matrix_multiply(z_reg, kz_mpc_orig_trim, /btranspose))
    z_exp_zregf = exp(-1.*dcomplex(0,1)*matrix_multiply(z_reg, kz_fine, /btranspose))


    ;; This one uses the true comov_dist_los values, so it is a little more different compared to the fft
    ;; (which run large to small)
    ;; in some simple testing, this had a difference ratio vs the FFT (abs(diff)/abs(fft)) of ~10^1
    z_true = reverse(comov_dist_los-min(comov_dist_los))
    z_exp_ztrue =  exp(-1.*dcomplex(0,1)*matrix_multiply(reverse(comov_dist_los-min(comov_dist_los)), kz_mpc_orig_trim, /btranspose))
    z_exp_ztruef =  exp(-1.*dcomplex(0,1)*matrix_multiply(reverse(comov_dist_los-min(comov_dist_los)), kz_fine, /btranspose))

    case dft_z_use of
        'true': z_exp = z_exp_ztrue
        'regular': z_exp = z_exp_zreg
    endcase

    case dft_z_use of
        'true': z_expf = z_exp_ztruef
        'regular': z_expf = z_exp_zregf
    endcase

    case dft_z_use of
        'true': z_use = z_true
        'regular': z_use = z_reg
    endcase

    n_kz_trim = n_elements(kz_mpc_orig_trim)
    n_kz_fine = n_elements(kz_fine)
    ;; data_ones has shape n_Freq)
    data_ones = replicate(1, n_freq)
    ;; z_exp has shape (n_freq, n_kz)
    data_ones_ft_fine = z_mpc_delta * matrix_multiply(data_ones, z_expf)
    data_ones_ft = z_mpc_delta * matrix_multiply(data_ones, z_exp)

    ;; generate plot
    if n_elements(save_path) ne 0 then begin
        cgps_open, save_path, /encapsulated, charsize=0.7
    endif
    
    ; plot final freqs and z array as scatter
    ; get buffers for scatter plots xrange
    z_buff = (max(z_use) - min(z_use)) / 30
    f_buff = (max(frequencies) - min(frequencies)) / 30

    inds = indgen(24, start=4, increment=n_freq/24)
    cb_data = replicate(0, n_freq)
    cb_data[inds] = 1
    cb_fft = fft(cb_data) * n_freq * z_mpc_delta
    cb_inds = where(abs(cb_fft)>0)
    print, kz_mpc_orig_trim
    cgplot, frequencies, data_ones - 1.2 , psym=6, linestyle=6, position=[.01, .8, .99, .88], symsize=.1, ytickformat="(A1)", $
        color='royal blue', yrange=[-.4, .4], xstyle=9, xrange=[min(frequencies) - f_buff, max(frequencies) + f_buff], xtitle='frequency (MHz)'
    cgaxis, xaxis=1, xstyle=9, xrange = [min(z_use)-z_buff, max(z_use)+z_buff],xtitle='dft z', /save
    cgplot, z_use, data_ones - .8, psym=6, linestyle=6, color='sky blue', symsize=.1, /over

    ; plot psf of fine transform

    cgplot, kz_fine, data_ones_ft_fine, ytitle='fine kz psf', xrange=[-.1,1], yrange=[-100,100],xtickformat="(A1)", $
        position=[.07, .33, .99, .71],/noerase
    ; cgplot, kz_mpc_orig_trim, cb_fft, /over, color='powder blue'
    cgplot, [0.30785424, 0.30785424], [-100, 100], /over, color='grey'
    cgplot, [0.61570847, 0.61570847], [-100, 100], /over, color='grey'
    cgplot, [0.92356277, 0.92356277], [-100, 100], /over, color='grey'

    cgplot, kz_mpc_orig_trim, data_ones_ft, linestyle=6, color='red', psym=6, symsize=.1, /over
    ; plot psf of coarse transform
    data_ft_scaled = abs(data_ones_ft) / max(abs(data_ones_ft))
    cgplot, kz_mpc_orig_trim, data_ft_scaled, xtitle='kz (mpc)', ytitle='abs scaled coarse psf', ylog=1, xrange=[-.1,1], yrange=[1e-4,1], $
        position=[.07, .07, .99, .31],/noerase
    
    cgplot, [0.30785424, 0.30785424], [1e-4,1], /over, color='grey'
    cgplot, [0.61570847, 0.61570847], [1e-4,1], /over, color='grey'
    cgplot, [0.92356277, 0.92356277], [1e-4,1], /over, color='grey'
    ; plot 
    ; for i = 0, n_elements(kz_mpc_orig_trim)-1 do begin
       ;  cgplot, [kz_mpc_orig_trim[i], kz_mpc_orig_trim[i]], !Y.CRange, /over, LineStyle=2, color='grey'
    ; endfor
    cgText, 0.5, .95, /normal, file_name, $
      charsize=cgdefcharsize()*1.25, alignment=0.5
    
    if n_elements(save_path) ne 0 then begin
        cgps_close, png=1, density=600, delete_ps=1
    endif 

end