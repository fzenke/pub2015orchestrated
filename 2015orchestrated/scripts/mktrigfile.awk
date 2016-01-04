#!/usr/bin/awk -f
BEGIN {
	active = -1
}
{
	if ( $2 == 1 ) {
		if ( active == -1 ) {
			active = $3
			tstart = $1
		} 
	} else { # wait until not active any more
	    if ( active > -1 ) {
			print tstart " " active " " $1-tstart
			active = -1
		}
	}
}
