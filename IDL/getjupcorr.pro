FUNCTION getJupCorr, filename
	mumin = 0.4
	k = 0.2

	;; Reads selected file and associated angle backplanes
	jup = readfits(filename, header, /SILENT)

	mufilename = repstr(filename, 'cmap', 'mu')
	mu = readfits(mufilename, /SILENT)
	mu0filename = repstr(filename, 'cmap', 'mu0')
	mu0 = readfits(mu0filename, /SILENT)

	;; In case of 360x180 pixel map, converts to 720x360
	if n_elements(jup[*, 0]) lt 700 then begin $
		jup = rebin(jup, 720, 360)
		mu = rebin(mu, 720, 360)
		mu0 = rebin(mu0, 720, 360)
	endif
;______________________________________________________________________________

	;; Performs angle corrections
	nonz = where(mu0*mu gt mumin and jup gt 0.0)
	jupcorr = jup
	jupcorr(nonz) = jupcorr(nonz) / mu(nonz)^k / mu0(nonz)^(1.0 - k)
;______________________________________________________________________________

	;; Creates NEB average value array

	yneb_0 = 204
	yneb_1 = 222

	yrange_length = yneb_1 - yneb_0 + 1

	x_length = n_elements(jupcorr[*, yneb_0])
	jupcorr_avg = dblarr(x_length)
	for y = yneb_0, yneb_1 do begin
		jupcorr_avg = jupcorr_avg + jupcorr[*, y]
	endfor
	jupcorr_avg = jupcorr_avg / yrange_length
	
;______________________________________________________________________________

	;; Shifts for convenience of display
	if (jupcorr_avg[0] gt mean(jupcorr_avg)) then nshift = 360 else nshift = 0

	jupcorr = shift(jupcorr, -nshift)
	jupcorr = normalize(jupcorr)
	return, jupcorr
END
