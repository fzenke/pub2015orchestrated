BEGIN {
	stimno = 0
}
{ 
	t = $1

	if ( t > timepoint ) {
		print stimno
		exit
	}

	active = $2
	stimno = $3
}

