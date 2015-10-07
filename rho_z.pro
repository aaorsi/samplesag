;Script to create a lookup table with rho and z

pro rho_z

	File = '../out/rho_z.txt'

	OmegaM = 0.28
	H0 = 100.0
	
	zMax = 20.
	zMin = 0.

	nZ   = 1000.
	binZ = (zMax - zMin)/(nZ-1)
	
	ZArr = findgen(nZ)*binZ + zMin

	openw,1,File
	for i = 0l,nZ-1 do begin
		z = ZArr[i]
		printf,1,z,dp(z,OM=OmegaM,H0 = h0)
	endfor
	close,1

end
