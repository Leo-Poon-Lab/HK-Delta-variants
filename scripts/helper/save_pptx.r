library(officer)

save_pptx <- function(plot=p, width=8, height=8/sqrt(2), file="tmp.pptx"){
	plot <- rvg::dml(ggobj = plot)
	read_pptx() %>%
	# add slide ----
	add_slide() %>%
	# specify object and location of object ----
	ph_with(plot, ph_location(width=width, height=height)) %>%
	# export slide -----
	base::print(
		target = file
	)
}