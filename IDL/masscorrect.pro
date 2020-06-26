PRO masscorrect, filename
	openr, lun, filename, /get_lun
	files = []
	attributes = []
	attr = ' '
	line = ' '
	while not eof(lun) do begin & $
		readf, lun, $
		format = '(A89, A)', line, attr & $
		line = strtrim(line, 2) & attr = strtrim(attr, 2) & $
			files = [files, line] & $
			attributes = [attributes, attr]
	endwhile
	free_lun, lun

	len = n_elements(files)
	idx = 0
	avgbrightness_setpoint = 0
	while (idx lt len) do begin
		jupCorr = getJupCorr(files[idx])
		if idx eq 0 then begin
			avgbrightness_setpoint = avg(jupCorr)
		endif else begin
			scalefactor = avg(jupCorr) / avgbrightness_setpoint
			jupCorr = jupCorr / scalefactor
		endelse
		file_orig_arr = strsplit(files[idx], '/', /EXTRACT)
		file_orig = file_orig_arr[n_elements(file_orig_arr)-1]
		file_new_arr = strsplit(file_orig, '.', /EXTRACT)
		file_new = file_new_arr[0] + '_corrected.cmap.fits'
		print, file_new
		writefits, file_new, jupCorr
		idx = idx + 1
	endwhile
END
