function create_plot_options, plot_options = plot_options, $
    hinv = hinv, plot_path = plot_path, plot_filebase = plot_filebase, $
    note = note, individual_plots = individual_plots, $
    png = png, eps = eps, pdf = pdf, return_new = return_new

  update_tags = list()
  update_values = list()
  if n_elements(plot_options) eq 0 then begin
    ;; default to hinv
    if n_elements(hinv) eq 0 then hinv = 1

    ;; default to combined plots:
    if n_elements(individual_plots) eq 0 then individual_plots=0

    ;; default to plotting to the screen
    if n_elements(png) eq 0 then png=0
    if n_elements(eps) eq 0 then eps=0
    if n_elements(pdf) eq 0 then pdf=0

    if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then begin
      pub = 1

      if (keyword_set(png) + keyword_set(eps) + keyword_set(pdf)) gt 1 then begin
        if keyword_set(png) then begin
          print, 'only one of eps, pdf and png can be set, using png'
          eps = 0
          pdf = 0
        endif else if keyword_set(pdf) then begin
          print, 'only one of eps, pdf and png can be set, using pdf'
          eps = 0
        endif
      endif

      if keyword_set(png) then begin
        plot_exten = '.png'
        delete_ps = 1
      endif else if keyword_set(pdf) then begin
        plot_exten = '.pdf'
        delete_ps = 1
      endif else if keyword_set(eps) then begin
        plot_exten = '.eps'
        delete_ps = 0
      endif
    endif else begin
      pub = 0
      plot_exten = ''
      delete_ps = 0
    endelse

    plot_options = {hinv: hinv, individual_plots: individual_plots, $
      png: png, eps: eps, pdf:pdf, pub: pub, plot_exten: plot_exten, $
      delete_ps: delete_ps}

  endif else begin

    if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then begin
      pub = 1

      if (keyword_set(png) + keyword_set(eps) + keyword_set(pdf)) gt 1 then begin
        if keyword_set(png) then begin
          print, 'only one of eps, pdf and png can be set, using png'
          eps = 0
          pdf = 0
        endif else if keyword_set(pdf) then begin
          print, 'only one of eps, pdf and png can be set, using pdf'
          eps = 0
        endif
      endif

      if keyword_set(png) then begin
        eps = 0
        pdf = 0
        plot_exten = '.png'
        delete_ps = 1
      endif else if keyword_set(pdf) then begin
        eps = 0
        png = 0
        plot_exten = '.pdf'
        delete_ps = 1
      endif else if keyword_set(eps) then begin
        png = 0
        pdf = 0
        plot_exten = '.eps'
        delete_ps = 0
      endif
    endif else begin
      ;; if the set one is turned off, turn off saving plots
      if (n_elements(png) gt 0 and plot_options.png gt 0) or $
          (n_elements(pdf) gt 0 and plot_options.pdf gt 0) or $
          (n_elements(eps) gt 0 and plot_options.eps gt 0) then begin
        pub = 0
        plot_exten = ''
      endif
    endelse

    if n_elements(hinv) gt 0 then begin
      update_tags.add, 'hinv'
      update_values.add, hinv
    endif
    if n_elements(individual_plots) gt 0 then begin
      update_tags.add, 'individual_plots'
      update_values.add, individual_plots
    endif
    if n_elements(png) gt 0 then begin
      update_tags.add, 'png'
      update_values.add, png
    endif
    if n_elements(pdf) gt 0 then begin
      update_tags.add, 'pdf'
      update_values.add, pdf
    endif
    if n_elements(eps) gt 0 then begin
      update_tags.add, 'eps'
      update_values.add, eps
    endif
    if n_elements(pub) gt 0 then begin
      update_tags.add, 'pub'
      update_values.add, pub
    endif
    if n_elements(plot_exten) gt 0 then begin
      update_tags.add, 'plot_exten'
      update_values.add, plot_exten
    endif
    if n_elements(delete_ps) gt 0 then begin
      update_tags.add, 'delete_ps'
      update_values.add, delete_ps
    endif
  endelse

  if n_elements(plot_path) gt 0 then begin
    update_tags.add, 'plot_path'
    update_values.add, plot_path
  endif
  if n_elements(plot_filebase) gt 0 then begin
    update_tags.add, 'plot_filebase'
    update_values.add, plot_filebase
  endif
  if n_elements(note) gt 0 then begin
    update_tags.add, 'note'
    update_values.add, note
  endif

  update_tags = update_tags.toarray()

  if n_elements(update_tags) gt 0 or keyword_set(return_new) then begin
    plot_options = update_option_struct(plot_options, update_tags, $
      update_values, return_new = return_new)
  endif

  return, plot_options
end
