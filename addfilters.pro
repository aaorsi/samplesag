; Script to add filters into the filters files

pro addfilters

	FilterData = '/fdg/aaorsi/LightCones/Stand/data/filters.dat'
	FilterInfo = '/fdg/aaorsi/LightCones/Stand/data/filters.info'

	Lastid = 223

	AddDir 	   = '/home/aaorsi/LSSTFilters/throughputs-1.2/baseline/rescaled_filters/'



	NewFilesArr = AddDir + ['rescaled_u.dat','rescaled_g.dat','rescaled_r.dat','rescaled_i.dat','rescaled_z.dat','rescaled_y4.dat']
	NewHead		= 'LSST March2012 total '+['u','g','r','i','z','y']
	NFilters = n_elements(NewFilesArr)
	
	openw,1,FilterData,/append
	openw,2,FilterInfo,/append

	_id = 0
	for i = 0,NFilters-1 do begin
		
		Filter = NewFilesArr[i]
		readcol,Filter,nn,tt,/sile
		aa = nn*10.0
		nlines =n_elements(aa)
				
		printf,1,'# '+NewHead[i]
		for j = 0,nlines-1 do $
			printf,1,aa[j],tt[j]
		printf,2,++Lastid,nlines,'		#'+NewHead[i]
	endfor

	close,1
	close,2
end


