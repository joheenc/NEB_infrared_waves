FUNCTION getNEB_nocorr, filename

	;; Reads selected file
	jup = readfits(filename, header, /SILENT)
	
	yneb_0 = 204
	yneb_1 = 222

	yrange_length = yneb_1 - yneb_0 + 1

	x_length = n_elements(jup[*, yneb_0])

	img = jup[*, yneb_0:yneb_1]
	img = normalize(img)

	RETURN, img
END


