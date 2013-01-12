pro pson, filename = filename, paper = paper, margin = margin, $
  page_size = page_size, inches = inches, aspect = aspect, $
  landscape = landscape, quiet = quiet, eps = eps

;- Check arguments
if (n_elements(filename) eq 0) then filename='idl.ps'
if (n_elements(paper) eq 0) then paper = 'LETTER'
if (n_elements(margin) eq 0) then begin
  margin = 2.5
endif else begin
  if keyword_set(inches) then margin = margin * 2.54
endelse

;- Check if Postscript mode is active
if !d.name eq 'PS' then begin
  message, 'POSTSCRIPT output is already active', /continue
  return
endif

;- Get ratio of character width/height to
;- screen width/height
xratio = float(!d.x_ch_size) / float(!d.x_vsize)
yratio = float(!d.y_ch_size) / float(!d.y_vsize)

;- Save current device information in common block
common pson_information, info
info = {device:!d.name, window:!d.window, font:!p.font, $
        background:!p.background, color:!p.color, filename:filename, $
        xratio:xratio, yratio:yratio}

;- Get size of page (centimeters)
widths  = [[ 8.5,  8.5, 11.0,  7.25] * 2.54, 21.0, 29.7]
heights = [[11.0, 14.0, 17.0, 10.50] * 2.54, 29.7, 42.0]
names   = ['LETTER', 'LEGAL', 'TABLOID', 'EXECUTIVE', $
  'A4', 'A3']
index = where(strupcase(paper) eq names, count)
if (count ne 1) then begin
  message, 'PAPER selection not supported', /continue
  return
endif
page_width  = widths[index[0]]
page_height = heights[index[0]]

;- If page size was supplied, use it
if (n_elements(page_size) eq 2) then begin
  page_width  = page_size[0]
  page_height = page_size[1]
  if keyword_set(inches) then begin
    page_width  = page_width * 2.54
    page_height = page_height * 2.54
  endif
endif

;- Compute aspect ratio of page when margins are subtracted
page_aspect = float(page_height - 2.0 * margin) / $
              float(page_width  - 2.0 * margin)

;- Get aspect ratio of current graphics window
if (!d.window ge 0) then begin
  win_aspect = float(!d.y_vsize) / float(!d.x_vsize)
endif else begin
  win_aspect = 512.0 / 640.0
endelse

;- If aspect ratio was supplied, use it
if (n_elements(aspect) eq 1) then $
  win_aspect = float(aspect)

;- Compute size of drawable area
;- (method used here is the same as the Printer method)
case keyword_set(landscape) of
  0 : begin
        if (win_aspect ge page_aspect) then begin
          ysize = page_height - 2.0 * margin
          xsize = ysize / win_aspect
        endif else begin
          xsize = page_width - 2.0 * margin
          ysize = xsize * win_aspect
        endelse
      end
  1 : begin
        if (win_aspect ge (1.0 / page_aspect)) then begin
          ysize = page_width - 2.0 * margin
          xsize = ysize / win_aspect
        endif else begin
          xsize = page_height - 2.0 * margin
          ysize = xsize * win_aspect
        endelse
      end
endcase

;- Compute offset of drawable area from page edges
;- (landscape method here is different than
;-  the Printer method)
if (keyword_set(landscape) eq 0) then begin
  xoffset = (page_width  - xsize) * 0.5
  yoffset = (page_height - ysize) * 0.5
endif else begin
  xoffset = (page_width  - ysize) * 0.5
  yoffset = (page_height - xsize) * 0.5 + xsize
endelse

;- Switch to Postscript device
;- Note (1): Default units are centimeters
set_plot, 'PS'
device, landscape=keyword_set(landscape), scale_factor=1.0
device, xsize=xsize, ysize=ysize, $
  xoffset=xoffset, yoffset=yoffset
device, filename=filename, /color, bits_per_pixel=8
if keyword_set(eps) then device, /encapsulated, preview = 2

;- Set character size
xcharsize = round(info.xratio * !d.x_vsize)
ycharsize = round(info.yratio * !d.y_vsize)
device, set_character_size=[xcharsize, ycharsize]

;- Report to user
if (keyword_set(quiet) eq 0) then $
  print, filename, $
    format='("Started POSTSCRIPT output to ", a)'

END
