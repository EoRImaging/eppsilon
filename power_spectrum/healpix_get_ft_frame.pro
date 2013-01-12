pro healpix_get_ft_frame, info_file

 info = read_info_file(info_file)

 nfiles = n_elements(info.files)
 
 ;; fix header problem from RTS files. Only needs to be done once, but
 ;; will just skip over if already done.
 if info.format eq 'fits' then rts_fits_cleanup, info.files, /quiet


 file_frequencies = dblarr(nfiles)
 file_npix = dblarr(nfiles)
 for i=0, nfiles-1 do begin
    map_file = info.files[i]

    if info.format eq 'fits' then begin
       fits_npix = getsize_fits(map_file, nmaps=nmaps, obs_npix = obs_npix, type = ftype, ext=0, nside=nside)
       read_fits_map, map_file, map_out, hdr, ehdr, nside=nside, ordering=ordering, coordsys=coordsys
       
       tbinfo, ehdr, ehdr_info    
       temp = strpos(strlowcase(ehdr_info.ttype), 'pixnum')
       wh = where(temp ne -1, count)
       if count eq 1 then begin
          pix_ind = wh[0]
          pix_nums = long(map_out[*,pix_ind])
       endif else message, 'Could not identify healpix pixel numbers'
       ;; check whether weights are included.
       dims = size(map_out, /dimensions)
       case dims[1] of
          6: begin
             wh = where(strlowcase(ehdr_info.ttype) eq 'stokes i' or strlowcase(ehdr_info.ttype) eq 'i', count)
             if count eq 1 then begin
                data_ind = wh[0]
                weight_ind = -1
             endif else message, 'Could not identify Stokes I column'
          end
          9: begin
             wh = where(strlowcase(ehdr_info.ttype) eq 'stokes i' or strlowcase(ehdr_info.ttype) eq 'i', count)
             wh2 = where(strlowcase(ehdr_info.ttype) eq 'iweight', count2)
             if count eq 1 and count2 eq 1 then begin
                data_ind = wh[0]
                weight_ind = wh2[0]
             endif else message, 'Could not identify Stokes I and I weight columns'
          end
          else: message, 'Unknown fits format'
       endcase
       
       ;; get frequency for each file
       fits_read, map_file, subint_data, subint_header, exten_no = 2, /no_pdu
       tbinfo, subint_header, subint_info
       
       freq_ind = where(subint_info.ttype eq 'FREQ')
       freq_list = tbget(subint_info, subint_data, freq_ind+1)
       ;; if info.freq_res gt 0 then $
       ;;    frequencies[i] = round((mean(freq_list) - info.freq_res/2.) / info.freq_res) * info.freq_res + info.freq_res/2. $
       ;; else frequencies[i] = mean(freq_list)
       file_frequencies[i] = mean(freq_list)
       ;; print, n_elements(freq_list), minmax(freq_list)
       
       file_npix[i] = obs_npix

    endif else begin
       ;; data is in text files not fits files
       pix_ind = info.text_cols.pix_ind
       data_ind = info.text_cols.data_ind
       weight_ind = info.text_cols.weight_ind

       ;; only need to read the pixnum column for now
       ncol = pix_ind+1
       fmt = ''
       for j=0, ncol-1 do begin
          if pix_ind eq j then fmt = fmt + 'L' else fmt = fmt + 'X'
          if j lt ncol-1 then fmt = fmt + ','
       endfor

       readcol, map_file, format = fmt, v1, count = nread, nlines = nlines, /silent
       if nread ne nlines then stop

       pix_nums = v1  
       obs_npix = n_elements(pix_nums)
       nside = 1024
       file_frequencies[i] = info.frequencies[i]
    endelse
 
    if i eq 0 then begin
       ref_nside = nside
       ref_obs_npix = obs_npix
       ref_pix_nums = reform(pix_nums, obs_npix, 1)
       ref_data_ind = data_ind
       ref_weight_ind = weight_ind
       nsets = 1
       min_pix = min(pix_nums)
    endif else begin
       if nside ne ref_nside then message, 'file ' + map_file + ' does not have the same nside as the first listed file.'
       if data_ind ne ref_data_ind or weight_ind ne ref_weight_ind then $
          message, 'file ' + map_file + ' does not have the same header format as the first listed file.'
       
       min_pix = min([min_pix, pix_nums])
       if obs_npix ne ref_obs_npix then begin
          if obs_npix lt ref_obs_npix then pix_nums = [pix_nums, lonarr(ref_obs_npix - obs_npix) -1] $
          else begin
             new_ref = lonarr(obs_npix, nsets) -1
             for j = 0, nsets-1 do new_ref[0:ref_obs_npix-1,j] = ref_pix_nums[*,j]
             
             ref_obs_npix = obs_npix
             ref_pix_nums = new_ref
          endelse
       endif
       
       match = 0
       for j=0, nsets -1 do if total(abs(pix_nums[sort(pix_nums)] - ref_pix_nums[sort(ref_pix_nums[*,j]),j])) eq 0 then match = 1
       
       if match eq 0 then begin
          print, 'This file (' + map_file + $
                 ') does not have the same pixels as any previous file. Only pixels that appear in all files will be used.'
          nsets = nsets + 1
          
          new_ref = lonarr(ref_obs_npix, nsets) -1
          for j=0, nsets -2 do new_ref[*,j] = ref_pix_nums[*,j]
          new_ref[*,nsets-1] = pix_nums
          
          ref_pix_nums = new_ref
       endif
       
    endelse
 endfor
 pix_nums=0

 n_freq = n_elements(info.frequencies)
 
 if nfiles gt 1 then begin
    file_freq_diffs = file_frequencies - shift(file_frequencies, 1)
    file_freq_diffs = file_freq_diffs[1:*]
    temp = floor(alog10(min(file_freq_diffs)))
    binsize = 10d ^ (temp - 1)
    file_freq_hist = histogram(file_freq_diffs, min = binsize/2d, binsize = binsize, locations = locs, reverse_indices = ri)
    wh = where(file_freq_hist gt 10*mean(file_freq_hist))
    file_freq_res = round(mean(file_freq_diffs[ri[ri[wh[0]]:ri[wh[0]+1]-1]])/binsize)*binsize
 endif else file_freq_res = info.frequencies

 if n_freq ne nfiles or info.frequencies[0] eq -1 then begin
    print, 'Number of frequencies in info file does not match number of files. Using frequencies from files.'
    frequencies = file_frequencies

    freq_resolution = file_freq_res
 endif else begin
    info_freq_diffs = info.frequencies - shift(info.frequencies, 1)
    info_freq_diffs = info_freq_diffs[1:*]
    
    temp = floor(alog10(min(info_freq_diffs)))
    binsize = 10d ^ (temp - 1)
    info_freq_hist = histogram(info_freq_diffs, min = binsize/2d, binsize = binsize, locations = locs, reverse_indices = ri)
    wh = where(info_freq_hist gt 10*mean(file_freq_hist))
    info_freq_res = round(mean(info_freq_diffs[ri[ri[wh[0]]:ri[wh[0]+1]-1]])/binsize)*binsize
    
    frequencies = info.frequencies
    freq_resolution = info_freq_res
 endelse

 ;; Now have to find common pixels among all sets
 if nsets gt 1 then begin
    ;; histogram all the pixel values and keep all the ones that appear nsets times
    pix_hist = histogram(reform(ref_pix_nums, nsets*ref_obs_npix), min=min_pix, locations =locs)
    wh = where(pix_hist eq nsets, count)
    if count gt 0 then pixels = locs[wh]
    pix_hist = 0
    wh = 0

    if min(pixels) lt 0 then pixels = pixels[where(pixels ge 0)]
 endif else pixels = reform(ref_pix_nums)

 ;; get pixel center theta/phi and limit to given ranges
 pix2ang_ring, nside, pixels, thetas, phis

 theta_range = info.theta_range *!pi/180d
 phi_range = info.phi_range *!pi/180d

 wh_pix_out = where(thetas lt theta_range[0] or thetas gt theta_range[1] or phis lt phi_range[0] or phis gt phi_range[1], $
                    n_pix_out, complement = wh_pix_in, ncomplement = n_pix_in)
 if n_pix_in eq 0 then message, 'no pixels in theta/phi range'
 if n_pix_out gt 0 then pixels = pixels[wh_pix_in]

 ;; replace code below with new routine. Haven't run this yet,
 ;; there will be bugs to fix.
 healpix_setup_ft,  pixels, nside, new_pix_vec, limits, degpix, kx_rad_vals, ky_rad_vals

 ;; pix2vec_ring, nside, pixels, pix_center_vec
 ;; ;; find mid point (work in x/y because of possible jumps in phi)
 ;; vec_mid = [mean(pix_center_vec[*,0]), mean(pix_center_vec[*,1]), mean(pix_center_vec[*,2])]
 ;; theta0 = acos(vec_mid[2])
 ;; phi0 = atan(vec_mid[1], vec_mid[0])

 ;; dists = sqrt((pix_center_vec[*,0]-vec_mid[0])^2d + (pix_center_vec[*,1]-vec_mid[1])^2d + (pix_center_vec[*,2]-vec_mid[2])^2d)
 ;; radius = max(dists)

 ;; disc_covers = 0
 ;; nloops = 0
 ;; while disc_covers lt 1 do begin
 ;;    query_disc, nside, vec_mid, radius, listpix, nlist, /inc
 ;;    min_pix = min([pixels, listpix])
 ;;    wh2 = where(histogram(listpix, min=min_pix) eq 0 and histogram(pixels, min=min_pix) gt 0, count2)
 ;;    if count2 gt 0 then radius = radius * (1 + (nloops+1)*0.1d) else disc_covers = 1
 ;;    nloops = nloops+1
 ;; endwhile
 
 ;; ;; remove pixels from listpix that are in my image -- only want nearby pixels not in my image
 ;; min_pix = min([pixels, listpix])
 ;; max_pix = max([pixels, listpix])
 ;; wh = where(histogram(listpix, min = min_pix, max = max_pix) gt 0 and histogram(pixels, min = min_pix, max = max_pix) eq 0, count)
 ;; if count gt 0 then outside_pix = wh + min_pix else stop
 ;; listpix=0
 ;; wh=0

 ;; pix2vec_ring, nside, outside_pix, out_center_vec

 ;; ;; define new coordinate system
 ;; ;; To get to current location, need to first rotate around z by
 ;; ;; phi, then around y by -theta, then around z by -phi
 ;; ;; use inverse to rotate back to vertical
 ;; rot_matrix = get_rot_matrix(theta0, phi0, /inverse)

 ;; new_pix_vec = rot_matrix ## pix_center_vec
 ;; new_out_vec = rot_matrix ## out_center_vec

 ;; if windowavailable(1) then wset, 1 else window, 1
 ;; surface, dist(5), /nodata, /save, xrange = [-1, 1], yrange = [-1, 1], zrange = [-1, 1], xtitle = 'x', ytitle = 'y'
 ;; plots, out_center_vec[*,0], out_center_vec[*,1], out_center_vec[*,2], psym = 4, color = 200, /T3D
 ;; plots, pix_center_vec[*,0], pix_center_vec[*,1], pix_center_vec[*,2], psym = 4, color = 254, /T3D

 ;; plots, new_out_vec[*,0], new_out_vec[*,1], new_out_vec[*,2], psym = 4, color = 100, /T3D
 ;; plots, new_pix_vec[*,0], new_pix_vec[*,1], new_pix_vec[*,2], psym = 4, color = 75, /T3D

 ;; pred_angle = healpix_rot(new_pix_vec[*,0], new_pix_vec[*,1])

 ;; x_rot = new_pix_vec[*,0] * cos(pred_angle) - new_pix_vec[*,1] * sin(pred_angle)
 ;; y_rot = new_pix_vec[*,0] * sin(pred_angle) + new_pix_vec[*,1] * cos(pred_angle)
 ;; x_out_rot = new_out_vec[*,0] * cos(pred_angle) - new_out_vec[*,1] * sin(pred_angle)
 ;; y_out_rot = new_out_vec[*,0] * sin(pred_angle) + new_out_vec[*,1] * cos(pred_angle)

 ;; lims = healpix_limits(x_rot, y_rot, x_out_rot, y_out_rot)

 ;; wh_inside = where(x_rot ge lims[0] and x_rot le lims[2] and y_rot ge lims[1] and y_rot le lims[3], count_in)
 ;; wh_out_inside = where(x_out_rot ge lims[0] and x_out_rot le lims[2] and y_out_rot ge lims[1] and y_out_rot le lims[3], $
 ;;                       count_out)
 ;; if count_out ne 0 then begin
 ;;    print, 'Amoeba failed to find limits that exclude all outside pixels. Set limits by hand.'
 ;;    allgood = 0
 ;; endif else begin
 ;;    if windowavailable(2) then wset, 2 else window, 2
 ;;    plot, x_rot, y_rot, /nodata, xrange = minmax(x_out_rot), yrange = minmax(y_out_rot), color=0
 ;;    oplot, x_out_rot, y_out_rot, psym = 3, color = 100
 ;;    oplot, x_rot, y_rot, psym = 3, color = 75
            
 ;;    x_range_plot = [lims[0], replicate(lims[2], 2), replicate(lims[0], 2)]
 ;;    y_range_plot = [replicate(lims[1], 2), replicate(lims[3], 2), lims[1]]
 ;;    oplot, x_range_plot, y_range_plot, psym = -3, color = 0

 ;;    print, n_elements(pixels)
 ;;    print, minmax(pixels)
 ;;    print, theta0, phi0
 ;;    print, vec_mid
 ;;    print, pred_angle*180/!pi
 ;;    print, lims
 ;;    print, count_in

 ;;    answer=''
 ;;    read, answer, prompt = 'Rotation and limits were calculated programatically, are they all right? (y/n)'
 ;;    case answer of
 ;;       'y': allgood = 1
 ;;       'yes': allgood = 1
 ;;       else: allgood = 0
 ;;    endcase
 ;; endelse

 ;; if allgood eq 1 then begin
 ;;    rot_angle = pred_angle
 ;;    limits = lims
 ;; endif else begin

 ;;    test_file = file_test(info.metadata_file) *  (1 - file_test(info.metadata_file, /zero_length))  
 ;;    if test_file eq 1 then begin
 ;;       ;; want to restore ONLY the rotation angle and limits set previously
 ;;       sObj = obj_new('IDL_Savefile', info.metadata_file) 
 ;;       sObj -> restore, ['rot_angle', 'limits']
 ;;       nloop = 1
 ;;    endif else begin
 ;;       rot_angle = pred_angle
 ;;       limits = [min(x_new), min(y_new), max(x_new), max(y_new)]
 ;;       nloop = 0
 ;;    endelse
    
 ;;    while allgood eq 0 do begin

 ;;       rot_good = 0
 ;;       while rot_good eq 0 do begin
 ;;          x_rot = new_pix_vec[*,0] * cos(rot_angle) - new_pix_vec[*,1] * sin(rot_angle)
 ;;          y_rot = new_pix_vec[*,0] * sin(rot_angle) + new_pix_vec[*,1] * cos(rot_angle)
 ;;          x_out_rot = new_out_vec[*,0] * cos(rot_angle) - new_out_vec[*,1] * sin(rot_angle)
 ;;          y_out_rot = new_out_vec[*,0] * sin(rot_angle) + new_out_vec[*,1] * cos(rot_angle)
 ;;         if nloop eq 0 then limits = [min(x_rot), min(y_rot), max(x_rot), max(y_rot)]
                   
 ;;          if windowavailable(1) then wset, 1 else window, 1
 ;;          plot, x_rot, y_rot, /nodata, xrange = minmax(x_out_rot), yrange = minmax(x_out_rot), psym = 3, color = 0
 ;;          oplot, x_out_rot, y_out_rot, psym = 3, color = 254
 ;;          oplot, x_rot, y_rot, psym = 3, color = 0
         
 ;;          x_range_plot = [limits[0], replicate(limits[2], 2), replicate(limits[0], 2)]
 ;;          y_range_plot = [replicate(limits[1], 2), replicate(limits[3], 2), limits[1]]
 ;;          oplot, x_range_plot, y_range_plot, psym = -3, color = 0
         
 ;;          rot_str = strsplit(string(rot_angle*180d/!pi), /extract)
 ;;          answer=''
 ;;          read, answer, prompt = 'Current rotation is ' + rot_str + '. Is this optimal? (y/n) '
 ;;          case answer of
 ;;             'y': rot_good = 1
 ;;             'yes': rot_good = 1
 ;;             else: begin
 ;;                new_rot = 0d           
 ;;                valid = 0           
 ;;                while valid eq 0 do begin  
 ;;                   on_ioerror, bad_rot  
 ;;                   read, new_rot, prompt = 'Enter new rotation angle in degrees: '
 ;;                   ;; If we get here, new_rot is a good number.
 ;;                   valid = 1  
 ;;                   bad_rot: if ~ valid then print, 'You entered an invalid number.'  
 ;;                endwhile
 ;;                rot_angle = new_rot *!pi/180d
 ;;             end
 ;;          endcase
 ;;       endwhile
       
 ;;       if nloop eq 0 then limits = [min(x_rot), min(y_rot), max(x_rot), max(y_rot)]
 ;;       limitsgood = intarr(4)
       
 ;;       lim_str = ['x min', 'y min','x max',  'y max']
 ;;       for i=0, 3 do begin
          
 ;;          while limitsgood[i] lt 1 do begin        
 ;;             plot, x_rot, y_rot, xrange = xrange, yrange = yrange, psym = 3, color = 0
 ;;             oplot, x_out_rot, y_out_rot, psym = 3, color = 254
            
 ;;             x_range_plot = [limits[0], replicate(limits[2], 2), replicate(limits[0], 2)]
 ;;             y_range_plot = [replicate(limits[1], 2), replicate(limits[3], 2), limits[1]]
 ;;             oplot, x_range_plot, y_range_plot, psym = -3, color = 0
             
 ;;             answer=''
 ;;             read, answer, prompt = 'Current ' + lim_str[i] + ' is ' + strsplit(string(limits[i]), /extract) + $
 ;;                   '. Is this optimal? (y/n) '
 ;;             case answer of
 ;;                'y': limitsgood[i] = 1
 ;;                'yes': limitsgood[i] = 1
 ;;                else: begin
 ;;                   new_lim = 0d           
 ;;                   valid = 0           
 ;;                   while valid eq 0 do begin  
 ;;                      on_ioerror, bad_lim  
 ;;                      read, new_lim, prompt = 'Enter new ' + lim_str[i] + ': '
 ;;                      ;; If we get here, new_lim is a good number.
 ;;                      valid = 1  
 ;;                      bad_lim: if ~ valid then print, 'You entered an invalid number.'  
 ;;                   endwhile
 ;;                   limits[i] = new_lim
 ;;                end
 ;;             endcase
 ;;          endwhile
 ;;       endfor
       
 ;;       answer=''
 ;;       read, answer, prompt = 'Are the rotation and all limits right? (y/n)'
 ;;       case answer of
 ;;          'y': allgood = 1
 ;;          'yes': allgood = 1
 ;;          else: ;; do nothing
 ;;       endcase
       
 ;;       nloop = nloop+1
 ;;    endwhile
 ;; endelse

 ;; x_new = x_rot
 ;; y_new = y_rot

 ;; limit pixels to those inside limits
 wh_inside = where(x_new ge limits[0] and x_new le limits[2] and y_new ge limits[1] and y_new le limits[3], count)
 pix_inside = pixels[wh_inside]
 x_rad_vals = x_new[wh_inside]
 y_rad_vals = y_new[wh_inside]
 npix_rad = count

 ;; get current colortable so it can be restored later
 tvlct, r, g, b, /get
 loadct,39

 wset, 1
 plot, x_new, y_new, psym = 3, color = 0
 oplot, x_rad_vals,  y_rad_vals, psym = 3, color = 75

 ;; Calculate k step size and range, which are given by spatial resolution & size of field
 ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
 ang_resolution = sqrt(3d/!pi) * 3600d/nside * (1d/60d) * (!pi/180d)

 x_rad_length = limits[2] - limits[0] + ang_resolution
 y_rad_length = limits[3] - limits[1] + ang_resolution
 
 kxy_rad_range = (2d*!pi) / ang_resolution
 kx_rad_delta = (2d*!pi) / x_rad_length
 ky_rad_delta = (2d*!pi) / y_rad_length

 ;; define locations (in k) to take FT
 n_kx = round(kxy_rad_range/kx_rad_delta) + 1
 n_ky = round(kxy_rad_range/ky_rad_delta) + 1
 if (ceil(n_kx/2d)-floor(n_kx/2d)) gt 0 then kx_rad_vals = (dindgen(n_kx)-floor(n_kx/2d)) * kx_rad_delta $
 else kx_rad_vals = (dindgen(n_kx)-n_kx/2+1) * kx_rad_delta
 if (ceil(n_ky/2d)-floor(n_ky/2d)) gt 0 then ky_rad_vals = (dindgen(n_ky)-floor(n_ky/2d)) * ky_rad_delta $
 else ky_rad_vals = (dindgen(n_ky)-n_ky/2+1) * ky_rad_delta
 
 ;; convert to mK (from Jy/beam)
 ;; pixel size in strad. in HEALPix is given by Omega = pi/(3 Nside^2)
 pix_size_str = !pi / (3d * nside^2)
 
 ;; calculate beam diameter in radians = c/freq*max_baseline
 max_baseline = 342.497
 ;; beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
 ;; beam_area_str = !pi * beam_diameter_rad^2d /4d

 ;; powers of ten are from: Jy, c^2, mK, MHz, kB
 ;; mk_conv_factors = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2 * 1.38)
 mk_conv_factors = dblarr(n_freq) + 2d * max_baseline^2d / (!dpi * 1.38)

 pixels = pix_inside
 save, file = info.metadata_file, pix_ind, data_ind, weight_ind, rot_angle, limits, pixels, x_rad_vals, y_rad_vals, ang_resolution, $
       frequencies, freq_resolution, mk_conv_factors, kx_rad_vals, ky_rad_vals

 tvlct, r, g, b

end
