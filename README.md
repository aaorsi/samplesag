# samplesag
Extracts galaxy properties from a SAG run, including observables such as galaxy luminosities and magnitudes

Checklist for version 1 release:
- Make readsag_hdf5.c to read all relevant variables in the hdf5 files
- Compute zcold in get_zcold.c 
- Compute the line luminosities, including IGM extinction
- Write running script in a clever way. 
	- Using a python script, perhaps?
