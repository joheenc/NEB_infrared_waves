; Tyler Hackett
; Apr 27 2019
;
; Compares the intensity of NEB waves to their latitudinal coverage.
;

inputFile  = '~/datafiles/cmaps_waves_classified.txt'
;inputFile  = 'out_test_cleaned.txt'

;
; ==== Get image data ====
;

openr, lun, inputFile, /get_lun
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

	label = ''
	imgFile = ''

numImages = n_elements(files)-1

;print, numImages, FORMAT='%d images'

imgIdx = 0 ;files[0] is null; files[1] is true beginning of array.

;; SET DESIRED REFERENCE DATE
date0 = 365.25 * 1995

a1995 = 365.25 * 1995 - date0
a1996 = 365.25 * 1996 - date0
a1997 = 365.25 * 1997 - date0
a1998 = 365.25 * 1998 - date0
a1999 = 365.25 * 1999 - date0
a2000 = 365.25 * 2000 - date0
a2001 = 365.25 * 2001 - date0
a2002 = 365.25 * 2002 - date0
a2003 = 365.25 * 2003 - date0
a2004 = 365.25 * 2004 - date0
a2005 = 365.25 * 2005 - date0
a2006 = 365.25 * 2006 - date0
a2007 = 365.25 * 2007 - date0
a2008 = 365.25 * 2008 - date0
a2009 = 365.25 * 2009 - date0
a2010 = 365.25 * 2010 - date0
a2011 = 365.25 * 2011 - date0
a2012 = 365.25 * 2012 - date0
a2013 = 365.25 * 2013 - date0
a2014 = 365.25 * 2014 - date0
a2015 = 365.25 * 2015 - date0
a2016 = 365.25 * 2016 - date0
a2017 = 365.25 * 2017 - date0
a2018 = 365.25 * 2018 - date0
a2019 = 365.25 * 2019 - date0
a2020 = 365.25 * 2020 - date0

;; NEB expansions:
startdates = [a1996 + 3 * 30.42, a1999, a2004 + 3 * 30.42, a2009 + 4 * 30.42, $
 a2012 + 2 * 30.42, a2015 + 5 * 30.42, a2017]
enddates = [a1997 + 10 * 30.42, a2002 + 4 * 30.42, a2007 + 3 * 30.42, a2011, $
 a2013, a2016 + 4 * 30.42, a2019]

numExpansions = n_elements(startdates)

timestamps  = fltarr(numImages)
expansions  = fltarr(numImages)
nebPower_16 = fltarr(numImages)
nebPower_17 = fltarr(numImages)

while (imgIdx le numImages) do begin

	imgIdx = imgIdx + 1 ;array begins at files[1]
	if imgIdx eq numImages then break
	filename = files[imgIdx]

	header = headfits(filename, /silent)
	datestr = sxpar(header, 'DATE_OBS')
	typestr = strlen(datestr)
	if typestr eq 10 then begin &
		year = strmid(datestr, 0, 4)
		month = strmid(datestr, 5, 2)
		day = strmid(datestr, 8, 2)
	endif
	if typestr eq 8 then begin &
		year = strjoin([strtrim(19, 2), strtrim(strmid(datestr, 6, 2), 2)])
		month = strmid(datestr, 3, 2)
		day = strmid(datestr, 0, 2)
	endif
	timestamps[imgIdx] = year * 365.25 + (month - 1) * 30.42 + (day - 1) - date0

	for i = 0, numExpansions-1 do begin

		if (timestamps[imgIdx] gt startdates[i]) and $
			(timestamps[imgIdx] lt enddates[i]) then begin

			expansions[imgIdx] = 1
			break

		endif

	endfor
	;; Reads selected file and associated angle backplanes
	jup = readfits(filename, header, /SILENT)

	mufilename = repstr(filename, 'cmap', 'mu')
	mu = readfits(mufilename, /SILENT)
	mu0filename = repstr(filename, 'cmap', 'mu0')
	mu0 = readfits(mu0filename, /SILENT)

	;CATCH, errorStatus

	;if errorStatus ne 0 then continue

	;; In case of 360x180 pixel map, converts to 720x360
	if n_elements(jup[*, 0]) ne 720 then begin $
		jup = rebin(jup, 720, 360)
		mu = rebin(mu, 720, 360)
		mu0 = rebin(mu0, 720, 360)
	endif
;______________________________________________________________________________

	;; Performs angle corrections
	mumin = 0.4
	nonz = where(mu0*mu gt mumin and jup gt 0.0)
	jupcorr = jup
	k = 0.2
	jupcorr(nonz) = jupcorr(nonz) / mu(nonz)^k / mu0(nonz)^(1.0 - k)
;______________________________________________________________________________

	;; Creates NEB average value array

	yneb_0 = 204
	yneb_1 = 226

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

	jupcorr_ave = shift(jupcorr_avg, -nshift)
	jupcorr     = shift(jupcorr, -nshift)

	midneb = (yneb_1 + yneb_0) / 2
	mu0e = shift(mu0[*, midneb], -nshift)
	mue = shift(mu[*, midneb], -nshift)

	fnshift = float(nshift)
	longitude = 360.0 + 0.5d0 * fnshift - 0.5d0 * findgen(720) ;; !!!Backwards: left-to-right degrees
;______________________________________________________________________________

	;; Restricts intensity to observed longitudes
	good = where((mu0e ^ k * mue ^ (1 - k)) gt 0.7) ;; OBSERVED

	;; If you uncomment the tvscl linei below, you will 
	;; see the angle-corrected region of the cmap:
	
	;tvscl, jupcorr

	;; Extract the NEB region from jupcorr
	neb = jupcorr[good, 208:222]

	neb = normalize(neb)

	; Compute the mean of each column of neb.
	; The dimensions after computing the mean
	; should match the number of columns in
	; the original array. If it does not, use
	; DIMENSIONS=1 instead of DIMENSIONS=2
	print, size(neb, /DIMENSIONS)
	neb = mean(neb, DIMENSION=2)
	print, size(neb, /DIMENSIONS)

	;; NOTE: if neb is not 720 elements wide, consider
	;; padding the array with zeros. The following code
	;; will assume neb is 720 elements wide.

	; Get the power spectrum
	; See also: fft_powerspectrum()
	nebPowerSpectrum = fft(neb)^2 / n_elements(neb)

	; Now you must extract the power of the relevant wavenumbers.
	; Check my notes. I may have left some info on the wavenumbers.
	; For now, I will assume it was 16 and 17:
	nebPower_16[imgIdx-1] = 2.0 * nebPowerSpectrum[16]
	nebPower_17[imgIdx-1] = 2.0 * nebPowerSpectrum[17]

	;; Note: the factor of 2.0 is just convention, since the power
	;; is split in half across the positive and negative frequencies.

endwhile

;; You now have arrays containing the power of the waves in the NEB.
;; You also have a corresponding array of timestamps for each of these.
;; For example: nebPower_16[100] will correspond to the power of the
;; waves on the date timestamp[100]. Using these values, we can create
;; a time-series plot of the wave power over time.

;; I made a function for creating the kinds of plots you may have
;; seen in my or Sandra's previous reports. It's just a matter of
;; tracking it down hopefully. It may be this one:
;; /home/tylerjh/Fall_2018/waves/neb_scatter.pro

end
