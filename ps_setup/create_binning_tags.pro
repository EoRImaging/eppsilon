function create_binning_tags, file_struct_arr=file_struct_arr, $
    binning_2d_options = binning_2d_options, binning_1d_options = binning_1d_options, $
    plot_2d_options = plot_2d_options, plot_options = plot_options

    plot_fadd_2d = ''
    plot_fadd_kslice = ''
    if plot_2d_options.kperp_linear_axis and plot_2d_options.kpar_linear_axis then begin
        plot_fadd_2d = plot_fadd_2d + '_linaxes'
    endif else begin
        if plot_2d_options.kperp_linear_axis then begin
        plot_fadd_2d = plot_fadd_2d + '_kperplinaxis'
        ;; Right now, if kperp_linear_axis is set, the power slices will all have linear axes
        ;; basically kpar_linear_axis is ignored for slice plots.
        ;; uvf slice plots always have linear axes
        plot_fadd_kslice = '_linaxes'
        endif
        if plot_2d_options.kpar_linear_axis then plot_fadd_2d = plot_fadd_2d + '_kparlinaxis'
    endelse

    fadd_2dbin = ''
    if n_elements(binning_2d_options) gt 0 then begin
        if binning_2d_options.no_kzero then fadd_2dbin = fadd_2dbin + '_nok0'
        if binning_2d_options.log_kpar then fadd_2dbin = fadd_2dbin + '_logkpar'
        if binning_2d_options.log_kperp then fadd_2dbin = fadd_2dbin + '_logkperp'
    endif

    fadd_1dbin = ''
    if binning_1d_options.log_k then fadd_1dbin = fadd_1dbin + '_logk'

    fadd_kpar_1d = fadd_1dbin
    fadd_kperp_1d = fadd_1dbin

    fadd_1dmask = ''
    if tag_exist(binning_1d_options, 'kperp_range_1dave') then begin
        fadd_1dmask = fadd_1dmask + '_kperp' + $
            number_formatter(binning_1d_options.kperp_range_1dave[0]) + '-' + $
            number_formatter(binning_1d_options.kperp_range_1dave[1])
        note_1d = 'kperp: [' + $
            number_formatter(binning_1d_options.kperp_range_1dave[0]) + ',' + $
            number_formatter(binning_1d_options.kperp_range_1dave[1]) + ']'

        ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
        ;; Convert to 1/Mpc for internal code usage
        if plot_options.hinv then begin
            kperp_range_1d_use = binning_1d_options.kperp_range_1dave * hubble_param
        endif else begin
            kperp_range_1d_use = binning_1d_options.kperp_range_1dave
        endelse
    endif else begin
        if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') then begin
            fadd_1dmask = fadd_1dmask + '_kperplambda' + $
                number_formatter(binning_1d_options.kperp_range_lambda_1dave[0]) + '-' + $
                number_formatter(binning_1d_options.kperp_range_lambda_1dave[1])
            note_1d = 'kperp: [' + $
                number_formatter(binning_1d_options.kperp_range_lambda_1dave[0]) + ',' + $
                number_formatter(binning_1d_options.kperp_range_lambda_1dave[1]) + ']'
        endif else begin
            ;; if no range set default to same range as is used in 2D plots
            kperp_range_lambda_1dave = [5., min([file_struct_arr.kspan/2., $
                file_struct_arr.max_baseline_lambda])]
            binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
                kperp_range_lambda_1dave = kperp_range_lambda_1dave)
        endelse
    endelse

    if tag_exist(binning_1d_options, 'kx_range_1dave') then begin
        fadd_1dmask = fadd_1dmask + '_kx' + $
            number_formatter(binning_1d_options.kx_range_1dave[0]) + '-' + $
            number_formatter(binning_1d_options.kx_range_1dave[1])
        note_1d = 'kx: [' + $
            number_formatter(binning_1d_options.kx_range_1dave[0]) + ',' + $
            number_formatter(binning_1d_options.kx_range_1dave[1]) + ']'

        ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
        ;; Convert to 1/Mpc for internal code usage
        if plot_options.hinv then begin
            kx_range_1d_use = binning_1d_options.kx_range_1dave * hubble_param
        endif else begin
            kx_range_1d_use = binning_1d_options.kx_range_1dave
        endelse
    endif else begin
        if tag_exist(binning_1d_options, 'kx_range_lambda_1dave') then begin
            fadd_1dmask = fadd_1dmask + '_kxlambda' + $
                number_formatter(binning_1d_options.kx_range_lambda_1dave[0]) + '-' + $
                number_formatter(binning_1d_options.kx_range_lambda_1dave[1])
            note_1d = 'kx: [' + $
                number_formatter(binning_1d_options.kx_range_lambda_1dave[0]) + ',' + $
                number_formatter(binning_1d_options.kx_range_lambda_1dave[1]) + ']'
        endif
    endelse

    if tag_exist(binning_1d_options, 'ky_range_1dave') then begin
        fadd_1dmask = fadd_1dmask + '_ky' + $
            number_formatter(binning_1d_options.ky_range_1dave[0]) + '-' + $
            number_formatter(binning_1d_options.ky_range_1dave[1])
        note_1d = 'ky: [' + $
            number_formatter(binning_1d_options.ky_range_1dave[0]) + ',' + $
            number_formatter(binning_1d_options.ky_range_1dave[1]) + ']'

        ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
        ;; Convert to 1/Mpc for internal code usage
        if plot_options.hinv then begin
            ky_range_1d_use = binning_1d_options.ky_range_1dave * hubble_param
        endif else begin
            ky_range_1d_use = binning_1d_options.ky_range_1dave
        endelse
    endif else begin
        if tag_exist(binning_1d_options, 'ky_range_lambda_1dave') then begin
            fadd_1dmask = fadd_1dmask + '_kylambda' + $
                number_formatter(binning_1d_options.ky_range_lambda_1dave[0]) + '-' + $
                number_formatter(binning_1d_options.ky_range_lambda_1dave[1])
            note_1d = 'ky: [' + $
                number_formatter(binning_1d_options.ky_range_lambda_1dave[0]) + ',' + $
                number_formatter(binning_1d_options.ky_range_lambda_1dave[1]) + ']'
        endif
    endelse

    if tag_exist(binning_1d_options, 'kperp_range_lambda_kparpower') then begin
        fadd_kpar_1d = fadd_kpar_1d + '_kperplambda' + $
            number_formatter(binning_1d_options.kperp_range_lambda_kparpower[0]) + '-' + $
            number_formatter(binning_1d_options.kperp_range_lambda_kparpower[1])
        note_kpar_1d = 'kperp: [' + $
            number_formatter(binning_1d_options.kperp_range_lambda_kparpower[0]) + ',' + $
            number_formatter(binning_1d_options.kperp_range_lambda_kparpower[1]) + ']'
    endif else begin
        fadd_kpar_1d = fadd_1dbin + fadd_1dmask
        if n_elements(note_1d) gt 0 then note_kpar_1d = note_1d

        if tag_exist(binning_1d_options, 'kperp_range_lambda_1dave') then begin
            binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
            kperp_range_lambda_kparpower = binning_1d_options.kperp_range_lambda_1dave)
        endif
    endelse

    if tag_exist(binning_1d_options, 'kpar_range_1dave') then begin
        fadd_1dmask = fadd_1dmask + '_kpar' + $
            number_formatter(binning_1d_options.kpar_range_1dave[0]) + '-' + $
            number_formatter(binning_1d_options.kpar_range_1dave[1])
        if n_elements(note_1d) eq 0 then begin
            note_1d=''
        endif else begin
            note_1d = note_1d + '; '
        endelse
        note_1d = note_1d + 'kpar: [' + $
            number_formatter(binning_1d_options.kpar_range_1dave[0]) + ',' + $
            number_formatter(binning_1d_options.kpar_range_1dave[1]) + ']'

        ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
        ;; Convert to 1/Mpc for internal code usage
        if plot_options.hinv then begin
            kpar_range_1d_use = binning_1d_options.kpar_range_1dave * hubble_param
        endif else begin
            kpar_range_1d_use = binning_1d_options.kpar_range_1dave
        endelse
    endif

    if not tag_exist(binning_1d_options, 'kpar_range_kperppower') and $
        tag_exist(binning_1d_options, 'kpar_range_1dave') then begin

        binning_1d_options = create_binning_1d_options(binning_1d_options = binning_1d_options, $
        kpar_range_kperppower = kpar_range_1dave)
    endif

    if tag_exist(binning_1d_options, 'kpar_range_kperppower') then begin
        fadd_kperp_1d = fadd_kperp_1d + '_kpar' + $
            number_formatter(binning_1d_options.kpar_range_kperppower[0]) + '-' + $
            number_formatter(binning_1d_options.kpar_range_kperppower[1])
        note_kperp_1d = 'kpar: [' + $
            number_formatter(binning_1d_options.kpar_range_kperppower[0]) + ',' + $
            number_formatter(binning_1d_options.kpar_range_kperppower[1]) + ']'

        ;; if we're plotting in [k]=h/Mpc then assume these numbers are in h/Mpc.
        ;; Convert to 1/Mpc for internal code usage
        if plot_options.hinv then begin
            kpar_range_kperppower_use = binning_1d_options.kpar_range_kperppower * hubble_param
        endif else begin
            kpar_range_kperppower_use = binning_1d_options.kpar_range_kperppower
        endelse
    endif
    fadd_1dbin = fadd_1dbin + fadd_1dmask

    binning_tags={fadd_2dbin: fadd_2dbin, fadd_1dbin: fadd_1dbin, fadd_1dmask: fadd_1dmask, $
        fadd_kpar_1d: fadd_kpar_1d, fadd_kperp_1d: fadd_kperp_1d, $
        plot_fadd_2d: plot_fadd_2d, plot_fadd_kslice: plot_fadd_kslice}

    if n_elements(kperp_range_1d_use) then binning_tags = create_struct(binning_tags, $
        'kperp_range_1d_use', kperp_range_1d_use)
    if n_elements(kpar_range_1d_use) then binning_tags = create_struct(binning_tags, $
        'kpar_range_1d_use', kpar_range_1d_use)
    if n_elements(kx_range_1d_use) then binning_tags = create_struct(binning_tags, $
        'kx_range_1d_use', kx_range_1d_use)
    if n_elements(ky_range_1d_use) then binning_tags = create_struct(binning_tags, $
        'ky_range_1d_use', ky_range_1d_use)
    if n_elements(kpar_range_kperppower_use) then binning_tags = create_struct(binning_tags, $
        'kpar_range_kperppower_use', kpar_range_kperppower_use)

    if n_elements(note_1d) then binning_tags = create_struct(binning_tags, $
        'note_1d', note_1d)
    if n_elements(note_kpar_1d) then binning_tags = create_struct(binning_tags, $
        'note_kpar_1d', note_kpar_1d)
    if n_elements(note_kperp_1d) then binning_tags = create_struct(binning_tags, $
        'note_kperp_1d', note_kperp_1d)

    return, binning_tags
end