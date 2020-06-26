
FUNCTION normalize, data
	
	max = max(data)
	min = min(data)

	RETURN, (data - min)/(max-min)

END


