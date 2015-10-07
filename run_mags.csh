#!/bin/tcsh -f

#$ -S /bin/tcsh
#$ -cwd



set RootDir = '/data2/aaorsi/SAG/Output/sag4_emlines/Cones_Stand'
set RunName	= 'Cones_Stand'
#set LightConeFile = '/fdg/ammunoz1/SUBMM2013/Conitos/Lightcone_aaorsi_130709_Catalog_130715_9.61sqdeg/lightcone_cube_coords_ord_1.dat'
set LightConeFile	= '/fdg/ammunoz1/SUBMM2013/Conitos/Lightcone_aaorsi_130709_Catalog_130725_9.61sqdeg/lightcone_cube_coords_ord_1.dat'
set rho_zFile	 = '../out/rho_z.txt'
#set OutFile = '/fdg/aaorsi/LightCones/Stand/out/mags_lightcone_v2'
set OutFile = '/fdg/aaorsi/LightCones/Stand/out/lightcone_igm'

set FilterInfo = '../data/filters.info'
set FilterData = '../data/filters.dat'
set VegaFile   = '/fdg/aaorsi/SAG/Version/sag4_emlines/data/A0V_KURUCZ_92.SED'

set NBinSed = 200

set NMags   = 6
set idMags  = (224 225 226 227 228 229)
set NameMags= (u g r i z y)
set NBoxes  = 64
set minStellarMass = 1.0e8
set IGMAtt	= 1

set OutputFormat = 'BOTH' # 'FITS' , 'HDF5' or 'BOTH'


echo $RootDir $RunName $rho_zFile $OutFile $FilterInfo $FilterData $VegaFile $LightConeFile $NBinSed $NMags $idMags $NBoxes $minStellarMass $IGMAtt $OutputFormat $NameMags

./get_mags $RootDir $RunName $rho_zFile $OutFile $FilterInfo $FilterData $VegaFile $LightConeFile \
$NBinSed $NMags $idMags $NBoxes $minStellarMass $IGMAtt $OutputFormat $NameMags
