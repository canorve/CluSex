#! /usr/bin/env python

import sys
import os.path


import configparser


def readcon(params,confile):


    convalues = read_config(confile)



    if convalues.get('run1'):
        params.run1 = convalues['run1']


    if convalues.get('deblend_nthresh1'):
        params.dn1 = convalues['deblend_nthresh1']

    if convalues.get('deblend_mincont1'):
        params.dm1 = convalues['deblend_mincont1']

    if convalues.get('analysis_thresh1'):
        params.at1 = convalues['analysis_thresh1']

    if convalues.get('detect_thresh1'):
        params.dt1 = convalues['detect_thresh1']

    if convalues.get('detect_minarea1'):
        params.da1 = convalues['detect_minarea1']

    if convalues.get('back_size1'):
        params.bs1 = convalues['back_size1']

    if convalues.get('run2'):
        params.run2 = convalues['run2']

    if convalues.get('deblend_nthresh2'):
        params.dn2 = convalues['deblend_nthresh2']

    if convalues.get('deblend_mincont2'):
        params.dm2 = convalues['deblend_mincont2']

    if convalues.get('analysis_thresh2'):
        params.at2 = convalues['analysis_thresh2']

    if convalues.get('detect_thresh2'):
        params.dt2 = convalues['detect_thresh2']

    if convalues.get('detect_minarea2'):
        params.da2 = convalues['detect_minarea2']

    if convalues.get('back_size2'):
        params.bs2 = convalues['back_size2']

    if convalues.get('back_filtersize2'):
        params.bf2 = convalues['back_filtersize2']


    if convalues.get('SatDs9'):
        params.satfileout = convalues['SatDs9']


    if convalues.get('SatScale'):
        params.satscale = convalues['SatScale']

    if convalues.get('SatOffset'):
        params.satoffset = convalues['SatOffset']

    if convalues.get('satlevel'):
        params.satlevel = convalues['satlevel']


    if convalues.get('makemask'):
        params.create = convalues['makemask']

    if convalues.get('scale'):
        params.scale = convalues['scale']


    if convalues.get('Offset'):
        params.offset = convalues['Offset']


    if convalues.get('OutCatalog'):
        params.output2 = convalues['OutCatalog']


    if convalues.get('MinSatSize'):
        params.minsatsize = convalues['MinSatSize']


    if convalues.get('SatQ'):
        params.satq = convalues['SatQ']


    if convalues.get('image'):
        params.image = convalues['image']


    if convalues.get('FracTol'):
        params.tol = convalues['FracTol']

    if convalues.get('ReduCoef'):
        params.red = convalues['ReduCoef']

    if convalues.get('SatMethod'):
        params.satmethod = convalues['SatMethod']

    if convalues.get('magzpt'):
        params.zpt = convalues['magzpt']


    if convalues.get('gain'):
        params.gain = convalues['gain']


    if convalues.get('pixscale'):
        params.plate = convalues['pixscale']


    if convalues.get('seeing'):
        params.seeing = convalues['seeing']


    if convalues.get('JoinScale'):
        params.joinscale = convalues['JoinScale']


    if convalues.get('ScaleCor'):
        params.scalecor = convalues['ScaleCor']


    if convalues.get('RegDs9'):
        params.regoutfile = convalues['RegDs9']


    




def read_config(confile):

	# Create a ConfigParser object
	config = configparser.ConfigParser()

	# Read the configuration file
	config.read(confile)

    # Access values from the configuration file
	image = config.get('General', 'image')
	magzpt = config.getfloat('General', 'MAG_ZEROPOINT')
	gain = config.getfloat('General', 'GAIN')
	pixscale= config.getfloat('General', 'PIXEL_SCALE')
	seeing = config.getfloat('General', 'SEEING_FWHM')
	makemask = config.getboolean('General', 'MakeMask')
	OutCatalog = config.get('General', 'OutCatalog')
	RegDs9 = config.get('General', 'RegDs9')
	run1 = config.getboolean('General', 'run1')
	run2 = config.getboolean('General', 'run2')


	deblend_nthresh1 = config.getfloat('Run1', 'DEBLEND_NTHRESH1')
	deblend_mincont1 = config.getfloat('Run1', 'DEBLEND_MINCONT1')
	analysis_thresh1 = config.getfloat('Run1', 'ANALYSIS_THRESH1')
	detect_thresh1 = config.getfloat('Run1', 'DETECT_THRESH1')
	detect_minarea1 = config.getfloat('Run1', 'DETECT_MINAREA1')
	back_size1 = config.getfloat('Run1', 'BACK_SIZE1')
	back_filtersize1 = config.getfloat('Run1', 'BACK_FILTERSIZE1')


	deblend_nthresh2 = config.getfloat('Run2', 'DEBLEND_NTHRESH2')
	deblend_mincont2 = config.getfloat('Run2', 'DEBLEND_MINCONT2')
	analysis_thresh2 = config.getfloat('Run2', 'ANALYSIS_THRESH2')
	detect_thresh2 = config.getfloat('Run2', 'DETECT_THRESH2')
	detect_minarea2 = config.getfloat('Run2', 'DETECT_MINAREA2')
	back_size2 = config.getfloat('Run2', 'BACK_SIZE2')
	back_filtersize2 = config.getfloat('Run2', 'BACK_FILTERSIZE2')


	scale = config.getfloat('Sizes', 'Scale', fallback=1)
	Offset = config.getfloat('Sizes', 'Offset', fallback=0)
	ReduCoef = config.getfloat('Sizes', 'ReduCoef')
	FracTol = config.getfloat('Sizes', 'FracTol')
	JoinScale = config.getfloat('Sizes', 'JoinScale')
	ScaleCor = config.getfloat('Sizes', 'ScaleCor')


	satlevel = config.getfloat('Saturation', 'SATUR_LEVEL')
	SatDs9 = config.get('Saturation', 'SatDs9')
	SatScale = config.getfloat('Saturation', 'SatScale')
	SatOffset = config.getfloat('Saturation', 'SatOffset')
	MinSatSize = config.getfloat('Saturation', 'MinSatSize')
	SatQ = config.getfloat('Saturation', 'SatQ')
	SatMethod = config.getint('Saturation', 'SatMethod')


	# Return a dictionary with the retrieved values
	config_values = {
		'image': image,
		'magzpt': magzpt,
		'gain': gain,
		'pixscale': pixscale,
		'seeing': seeing,
		'makemask': makemask,
		'OutCatalog': OutCatalog,
		'RegDs9': RegDs9,

		'run1': run1,
		'deblend_nthresh1': deblend_nthresh1,
		'deblend_mincont1': deblend_mincont1,
		'analysis_thresh1': analysis_thresh1,
		'detect_thresh1': detect_thresh1,
		'detect_minarea1': detect_minarea1,
		'back_size1': back_size1,
		'back_filtersize1': back_filtersize1,

		'run2': run2,
		'deblend_nthresh2': deblend_nthresh2,
		'deblend_mincont2': deblend_mincont2,
		'analysis_thresh2': analysis_thresh2,
		'detect_thresh2': detect_thresh2,
		'detect_minarea2': detect_minarea2,
		'back_size2': back_size2,
		'back_filtersize2': back_filtersize2,

		'scale': scale,
		'Offset': Offset,
		'ReduCoef': ReduCoef,
		'FracTol': FracTol,
		'JoinScale': JoinScale,
		'ScaleCor': ScaleCor,

		'satlevel': satlevel, 
		'SatDs9': SatDs9, 
		'SatScale': SatScale, 
		'SatOffset': SatOffset, 
		'MinSatSize': MinSatSize, 
		'SatQ': SatQ, 
		'SatMethod': SatMethod, 



	}

	return config_values


