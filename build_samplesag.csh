#!/bin/tcsh -f
#$ -S /bin/tcsh
#$ -cwd

set OutDir      = './'
set sag_hdf5dir = './'
set ranName     = 'test'

set snaplist    = (99 61)
set iv0         = 0
set ivf         = 63

set OutputFormat = 'ASCII'

#File locations
set FilterInfo = '../data/filters.info'
set FilterData = '../data/filters.dat'
set VegaFile   = '/fdg/aaorsi/SAG/Version/sag4_emlines/data/A0V_KURUCZ_92.SED'

set idMags  = (224 225 226 227 228 229)
set NameMags= (u g r i z y)
set AB      = (1 1 1 1 1 1)

set nsnap       = `echo $#snaplist` 
set nmags       = `echo $#idMags`
@ nvol          = $ivf - $iv0 + 1

./get_mags $sag_hdf5dir $ranName $nsnap $snaplist $nvol $iv0 $ivf \
           $OutputFormat $FilterInfo $FilterData $VegaFile \
           $nmags $idMags $NameMags $AB \
           props redshift mstellar gt 1e8 mag_uox mag_gox mag_zox \
           l_lyalphax l_oii3727x mcold lt 1e12


