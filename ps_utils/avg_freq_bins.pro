function avg_freq_bins, frequencies, num_coarse_bins, freq_bin_avg, freq_bin_mask=freq_bin_mask, avg_freqs=avg_freqs

    ;; compute number of frequency channels in each coarse bin
    if n_elements(frequencies) mod num_coarse_bins ne 0 then begin
        message, 'number of coarse bins must divide evenly into frequencies'
    endif else begin
        coarse_bin_width = n_elements(frequencies) / num_coarse_bins
    endelse

    ;; apply mask if applicable
    if n_elements(freq_bin_mask) ne 0 then begin
        if n_elements(freq_bin_mask) ne coarse_bin_width then begin
            message, 'freq_bin_mask must have length of number of coarse bin channels'
        endif
    endif

    avg_n_freqs = num_coarse_bins * n_elements(freq_bin_avg)
    temp_freqs = fltarr(avg_n_freqs)

    ;; perform binning over each coarse channel
    for i = 0, num_coarse_bins - 1 do begin
        bin_step = 0
        for j = 0, n_elements(freq_bin_avg) - 1 do begin
            freq_inds = findgen(freq_bin_avg[j]) + bin_step + (i * coarse_bin_width)
            this_freqs = frequencies(freq_inds)
            mask_inds = findgen(freq_bin_avg[j]) + bin_step
            this_mask = freq_bin_mask[mask_inds]
            wh_unflagged = where(this_mask eq 1, count_unflagged)
            if count_unflagged gt 0 then begin
                temp_freqs[j + i * n_elements(freq_bin_avg)] = mean(this_freqs[wh_unflagged])
            endif
            bin_step += freq_bin_avg[j]
        endfor
    endfor

    unflagged_freqs = where(temp_freqs ne 0, freq_use)
    if freq_use gt 0 then begin
        avg_freqs = temp_freqs[unflagged_freqs]
    endif else begin
        message, 'all input frequencies were flagged'
    endelse

    return, avg_freqs

end