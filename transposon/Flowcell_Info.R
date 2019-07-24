###
fancy_scientific <- function(l) { 
	# turn in to character string in scientific notation 
	l <- format(l, scientific = TRUE) 
	# quote the part before the exponent to keep all the digits 
	l <- gsub("^(.*)e", "'\\1'e", l) 
	# turn the 'e+' into plotmath format 
	l <- gsub("e", "%*%10^", l) 
	l <- gsub("\\+",'',l)
	# return this as an expression 
	parse(text=l) 
} 
