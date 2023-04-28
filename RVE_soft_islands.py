# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #


# / Code description \
# This python script is ready to use in Abaqus CAE finite elements software and generates a fully parameterized 
# representative volume element of a soft islands-patterned microstructure. Geometry, material and discretization parameters 
# can be changed by the user in the Variables section. The code automatically generates RVEs for all input parameter combinations 
# and handles meshing and application of the boundary conditions as specified.

# Contact: wilmers@uni-wuppertal.de


# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// IMPORT \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

from abaqus import *
from abaqusConstants import *
import __main__
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from textRepr import *
from math import *
import sys
import regionToolset
import displayGroupMdbToolset as dgm
import displayGroupOdbToolset as dgo
import connectorBehavior
import mesh
import os
from time import time, sleep
from datetime import datetime
import itertools
from random import *
import numpy as np
from numpy.linalg import norm



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// PRE \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

# make all commands readable
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# work directory
wdir = str(os.getcwd())+'\\'

# log function for both windows + abaqus console output
def log(text):
	print >> sys.__stdout__, text
	print text

start_0 = time()
log('')
log('### --------------- SOFT ISLANDS COHESIVE --------------- ###')
log('')
log('Working directory: '+str(wdir))



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// VARIABLES \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

#						model name
model_name			=	'RVE_soft_islands'


# / GEOMETRY \ #

#						soft island large radius
#						(default: 10.0 um)
si_r1_list			=	[10.0, ]

#						soft island tip radius
#						(default: 1.0 um)
si_r2_list			=	[1.0, ]

#						soft island length, r1 center to r2 center
#						(default: 18 um)
si_l_list			=	[18.0, ]

#						crack length
#						(default: 29.0 um)
l_list				=	[29.0, ]

#						hard substrate thickness
#						(default: 0.95 um)
t_hard_list			=	[0.95, ]

#						soft substrate thickness
#						(default: 10.0 um)
t_soft_list			=	[10.0, ]

#						crack separation
#						(default: 10.0 um)
d_list				=	[10.0, ]

#						substrate base layer thickness
#						(default: 8.0 um)
t_base_list			=	[8.0, ]

#						cohesive layer thickness
#						(default: 0.001 um)
t_coh_list			=	[0.001, ]



# / PDMS \ #

#						soft substrate young's modulus
#						(default: 2 MPa)
E_soft_list			=	[2.0, ]

#						ratio between soft and hard young's modulus
#						(default: 100)
E_factor_list		=	[100.0, ]

#						soft substrate poisson's ratio
#						(default: 0.48)
nu_soft_list		=	[0.48, ]

#						hard substrate poisson's ratio
#						(default: 0.48)
nu_hard_list		=	[0.48, ]

#						soft cohesive material fracture energy
#						(default: 34.0 J/m^2)
ener_soft_list		=	[34.0, ]

#						hard cohesive material fracture energy
#						(default: 1.0 J/m^2)
ener_hard_list		=	[1.0, ]

#						soft cohesive material maxs damage nominal stress
#						(default: 1.0 MPa)
s_soft_list			=	[1.0, ]

#						cohesive element viscosity, convergence stability control
#						(default: 0.001)
visc_list			=	[0.001, ]



# / SILVER LAYER DELAMINATION \ #

#						additional metal layer and cohesive zone between metal and hardened polymer
metal_layer			=	True

#						metal thickness
#						(default: 0.1 um)
t_metal_list		=	[0.1, ] 


#						metal layer young's modulus 	
#						(default: 63000.0 MPa)
E_metal				=	63000.0								

#						metal layer poisson's ratio	
#						(default: 0.37)
nu_metal			=	0.37		

#						cohesive interface layer thickness
#						(default: 0.001 um)
t_inter				=	0.001

#                       interface fracture energy
#                       (default: 0.01 J/m^2)
ener_interface_list =   [0.01, ]

#                       interface maxs damage nominal stress
#                       (default: 2 MPa)
s_interface_list	=   [2.0, ]



# / ADDITIONAL 3D STRUCTURE \ #

#							additional secondary 3D islands
add_secondary_islands	=	True

# 							True: structure with cohesive zone; False: structure without cohesive zone
sec_islands_cohesive	=	True

#							only if sec_islands_coh = False. Choose from 'polymer_soft', 'polymer_hard', 'polymer_var', 'metal'.			
sec_islands_material	=	'polymer_hard'

#							x center position of additional structure, middle = 0.5
#							(default: 0.5)
add_x_fact_list 		= 	[0.5, ]

#							y center position of additional structure, middle = 0.5
#							(default: 0.5)
add_y_fact_list			=	[0.5, ]

#							additional structure radius
#							(default: 5 um)
add_r_list				= 	[5.0, ]

#							additional structure young's modulus ratio
#							(default: 100)
E_var_factor_list		=	[100.0, ]



# / MESH \ #
seed_size_list			=	[1.3, ]						# mesh size 							(default: 1.3 um)
elnum_ha				=	4							# hard layer thickness element amount	(default: 4)
elnum_so				=	4							# soft layer thickness element amount	(default: 4)
elnum_dz				=	2							# base layer thickness element amount	(default: 2)



# / BOUNDARY CONDITIONS \ #
E11						=	0.0							# x strain
E22						=	0.0							# y strain
E12						=	0.0							# xy strain



# / ABAQUS \ #
max_inc					=	0.01						# maximum increment size
cpus					=	16							# cpus used for calculations
gpus					=	0							# gpus used for calculations
mem						=	90							# memory usage
run						=	False						# run job immediately
newdir					=	True						# create new folder





# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// ITERATONS \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

# create new folder
if newdir == True:
	today = datetime.now()
	timestr = today.strftime('%Y-%m-%d_%H-%M-%S_')
	os.mkdir(wdir+timestr+model_name)

# parameters used
RVE_count = 0

param_lists = [
	si_r1_list,
	si_r2_list,
	si_l_list,
	l_list,
	t_hard_list,
	t_soft_list,
	d_list,
	E_soft_list,
	E_factor_list,
	nu_soft_list,
	nu_hard_list,
	ener_soft_list,
	ener_hard_list,
	s_soft_list,
	visc_list,
	t_base_list,
	t_coh_list,
	seed_size_list,
	ener_interface_list,
	s_interface_list,
	t_metal_list,
	add_r_list,
	add_x_fact_list,
	add_y_fact_list,
	E_var_factor_list,
	]

# generate list of all parameter combinations
parameters = list(itertools.product(*param_lists))

# clean and prepare parameters file
open(wdir+str(model_name)+'_parameters.txt', 'w').close()

# setup start.bat
open(wdir+str(model_name)+'_start.bat', 'w').close()

# print RVE count
log('RVE amount: '+str(len(parameters)))

# generate RVEs
for i in range(len(parameters)):
	si_r1			= parameters[i][0]
	si_r2			= parameters[i][1]
	si_l			= parameters[i][2]
	cr_l			= parameters[i][3]
	hard_t			= parameters[i][4]
	soft_t			= parameters[i][5]
	d_y				= parameters[i][6]
	soft_E			= parameters[i][7]
	E_factor		= parameters[i][8]
	soft_nu			= parameters[i][9]
	hard_nu			= parameters[i][10]
	soft_ener		= parameters[i][11]
	hard_ener		= parameters[i][12]
	soft_s			= parameters[i][13]
	visc			= parameters[i][14]
	d_z				= parameters[i][15]
	coh_t			= parameters[i][16]
	seed_size		= parameters[i][17]
	interface_ener	= parameters[i][18]
	interface_s		= parameters[i][19]
	silver_t		= parameters[i][20]
	cylndr_r		= parameters[i][21]
	cylndr_x_fact	= parameters[i][22]
	cylndr_y_fact	= parameters[i][23]
	E_var_factor	= parameters[i][24]
	# dependent parameters
	#cr_l			= 2*d_y+9.0 				# this can be used for proportional increase of crack length and crack separation
	hard_s			= 0.15*soft_E*E_factor
	soft_K			= soft_E/coh_t
	hard_K			= soft_E*E_factor/coh_t
	inter_K			= hard_K

	RVE_count = RVE_count+1

	log('')
	log('START: RVE '+str(RVE_count)+' / '+str(len(parameters)))
	start_rve = time()

	# write parameter values to txt file
	paramtxt = open(wdir+str(model_name)+'_parameters.txt', 'a')
	paramtxt.write(
		'RVE: '				+str(RVE_count)			+', '+
		'si_r1: '			+str(si_r1)				+', '+
		'si_r2: '			+str(si_r2)				+', '+
		'si_l: '			+str(si_l)				+', '+
		'cr_l: '			+str(cr_l)				+', '+
		'hard_t: '			+str(hard_t)			+', '+
		'soft_t: '			+str(soft_t)			+', '+
		'd_y: '				+str(d_y)				+', '+
		'd_z: '				+str(d_z)				+', '+
		'coh_t: '			+str(coh_t)				+', '+
		'soft_E: '			+str(soft_E)			+', '+
		'E_factor: '		+str(E_factor)			+', '+
		'soft_nu: '			+str(soft_nu)			+', '+
		'hard_nu: '			+str(hard_nu)			+', '+
		'soft_ener: '		+str(soft_ener)			+', '+
		'hard_ener: '		+str(hard_ener)			+', '+
		'soft_s: '			+str(soft_s)			+', '+
		'hard_s: '			+str(hard_s)			+', '+
		'soft_K: '			+str(soft_K)			+', '+
		'hard_K: '			+str(hard_K)			+', '+
		'inter_K: '			+str(inter_K)			+', '+
		'visc: '			+str(visc)				+', '+
		'metal_layer: '		+str(metal_layer)		+', '+
		'silver_t: '		+str(silver_t)			+', '+
		'E_metal: '			+str(E_metal)			+', '+
		'nu_metal: '		+str(nu_metal)			+', '+
		'seed_size: '		+str(seed_size)			+', '+
		't_inter: '			+str(t_inter)			+', '+
		'interface_ener: '	+str(interface_ener)	+', '+
		'interface_s: '		+str(interface_s) 		+', '+
		'add_secondary_islands: '	+str(add_secondary_islands)	+', '+
		'sec_islands_cohesive: '	+str(sec_islands_cohesive)	+', '+
		'sec_islands_material: '	+str(sec_islands_material)	+', '+
		'cylndr_r: '		+str(cylndr_r) 			+', '+
		'cylndr_x_fact: '	+str(cylndr_x_fact) 	+', '+
		'cylndr_y_fact: '	+str(cylndr_y_fact) 	+', '+
		'E_var_factor: '	+str(E_var_factor)		+', '+
		'\n')
	paramtxt.close()



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// MODEL \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	job_name = str(model_name)+'_'+str(RVE_count)
	open(wdir+str(job_name)+'.cae', 'w').close()

	try:
		# close all models
		mdb.close()

		# create model and delete original model
		RVEmodel = mdb.Model(name=model_name)
		del mdb.models['Model-1']
	except:
		log('+ + + MODEL ERROR + + +')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// MATERIAL \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	# dependent 
	hard_E = soft_E * E_factor
	hard_E_var = soft_E * E_var_factor

	try:
		# substrate materials
		RVEmodel.Material(name='polymer_soft')
		soft_C10 = soft_E/(4*(1+soft_nu))
		soft_D1 = (6*(1-2*soft_nu))/soft_E
		RVEmodel.materials['polymer_soft'].Hyperelastic(
			materialType=ISOTROPIC,
			testData=OFF,
			type=NEO_HOOKE,
			volumetricResponse=VOLUMETRIC_DATA,
			table=((soft_C10, soft_D1),
			))

		RVEmodel.Material(name='polymer_hard')
		hard_C10 = hard_E/(4*(1+hard_nu))
		hard_D1 = (6*(1-2*hard_nu))/hard_E
		RVEmodel.materials['polymer_hard'].Hyperelastic(
			materialType=ISOTROPIC,
			testData=OFF,
			type=NEO_HOOKE,
			volumetricResponse=VOLUMETRIC_DATA,
			table=((hard_C10, hard_D1),
			))

		RVEmodel.Material(name='polymer_var')
		var_C10 = hard_E_var/(4*(1+hard_nu))
		var_D1 = (6*(1-2*hard_nu))/hard_E_var
		RVEmodel.materials['polymer_var'].Hyperelastic(
			materialType=ISOTROPIC,
			testData=OFF,
			type=NEO_HOOKE,
			volumetricResponse=VOLUMETRIC_DATA,
			table=((var_C10, var_D1),
			))


		# metal
		RVEmodel.Material(name='metal')
		RVEmodel.materials['metal'].Elastic(table=((E_metal, nu_metal), ))


		# cohesive materials
		RVEmodel.Material(name='polymer_soft_coh')
		RVEmodel.materials['polymer_soft_coh'].MaxsDamageInitiation(table=((soft_s, soft_s, soft_s), ))
		RVEmodel.materials['polymer_soft_coh'].maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((soft_ener, ), ))
		RVEmodel.materials['polymer_soft_coh'].Elastic(type=TRACTION,
			table=((soft_K, soft_K, soft_K), ))

		RVEmodel.Material(name='polymer_hard_coh')
		RVEmodel.materials['polymer_hard_coh'].MaxsDamageInitiation(table=((hard_s, hard_s, hard_s), ))
		RVEmodel.materials['polymer_hard_coh'].maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((hard_ener, ), ))
		RVEmodel.materials['polymer_hard_coh'].Elastic(type=TRACTION,
			table=((hard_K, hard_K, hard_K), ))
		
		RVEmodel.Material(name='interface_coh')
		RVEmodel.materials['interface_coh'].MaxsDamageInitiation(table=((interface_s, interface_s, interface_s), ))
		RVEmodel.materials['interface_coh'].maxsDamageInitiation.DamageEvolution(type=ENERGY, table=((interface_ener, ), ))
		RVEmodel.materials['interface_coh'].Elastic(type=TRACTION,
			table=((inter_K, inter_K, inter_K), ))

	except:
		log('+ + + MATERIAL ERROR + + +')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// GEOMETRY \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	# / DRAW SOFT ISLAND \ #

	try:

		""" Use genDropletCurve to draw the curve of a soft island into a sketch.
			Arguments:
				sketch_name		:	name of sketch you want to draw in
				orientation		:	left hand curve or right hand curve looking from top view, use 'RIGHT' or 'LEFT'
				center_coord_x	:	x coordinate of big radius circle origin
				center_coord_y	:	y coordinate of big radius circle origin
				si_r1			:	(GLOBAL VAR) big circle radius, no need to change
				si_r2			:	(GLOBAL VAR) small circle radius, no need to change
				si_l			:	(GLOBAL VAR) soft island center length, no need to change
				closed			:	True to draw closed shape, False to draw only outer curve
		"""

		def genDropletCurve(sketch_name, orientation, center_coord_x, center_coord_y, si_r1, si_r2, si_l, closed):

			# given coords
			x1 = center_coord_x + si_l
			y1 = center_coord_y
			x2 = center_coord_x
			y2 = center_coord_y

			# calculate angles
			gamma	= -atan( ( y2 - y1 )/float( ( x2 - x1 ) ) )
			beta	= atan( ( si_r1 - si_r2 )/float( sqrt( ( x2 - x1)**2 + ( y2 - y1 )**2 ) ) )
			alpha	= gamma - beta

			# calculate coords
			x3 = x1 + si_r2*abs( cos( pi/2.0 - alpha ) )
			x4 = x2 + si_r1*abs( cos( pi/2.0 - alpha ) )

			if orientation == 'RIGHT':
				y3 = y1 + si_r2*abs( sin( pi/2.0 - alpha ) )
				y4 = y2 + si_r1*abs( sin( pi/2.0 - alpha ) )
				rot = CLOCKWISE
			if orientation == 'LEFT':
				y3 = y1 - si_r2*abs( sin( pi/2.0 - alpha ) )
				y4 = y2 - si_r1*abs( sin( pi/2.0 - alpha ) )
				rot = COUNTERCLOCKWISE

			# draw geometry
			RVEmodel.sketches[sketch_name].ArcByCenterEnds(
				center=(x2, y2),
				point1=(x2-si_r1, y2),
				point2=(x4, y4),
				direction=rot)

			RVEmodel.sketches[sketch_name].ArcByCenterEnds(
				center=(x1, y1),
				point1=(x3, y3),
				point2=(x1+si_r2, y1),
				direction=rot)

			RVEmodel.sketches[sketch_name].Line(
				point1=(x4, y4),
				point2=(x3, y3))

			if closed == True:
				RVEmodel.sketches[sketch_name].Line(
					point1=(x1+si_r2, y1),
					point2=(x2-si_r1, y2))

	except:
		log('+ + + GEOMETRY FUNCTION ERROR + + +')


	# / PARTS \ #

	# RVE dimensions
	x_l = cr_l + si_l + si_r1 + si_r2
	y_l = 2 * d_y + 4 * si_r1


	# parameters
	cylndr_x_center = cylndr_x_fact * cr_l * 0.5
	cylndr_y_center = cylndr_y_fact * (si_r1+0.5*d_y)

	try:

		# SKETCHES

		# layer sketch standard
		l_sketch_standard = RVEmodel.ConstrainedSketch(name='layer_sketch_standard', sheetSize=1.0)
		RVEmodel.sketches['layer_sketch_standard'].Line(
			point1=(0, 0),
			point2=(0.5*cr_l, 0)
			)
		genDropletCurve('layer_sketch_standard', 'RIGHT', 0.5*cr_l+si_r1, 0.0, si_r1, si_r2, si_l, False)
		RVEmodel.sketches['layer_sketch_standard'].Line(
			point1=(0.5*cr_l+si_r1+si_l+si_r2, 0),
			point2=(cr_l+si_r1+si_l+si_r2, 0)
			)
		RVEmodel.sketches['layer_sketch_standard'].Line(
			point1=(cr_l+si_r1+si_l+si_r2, 0),
			point2=(cr_l+si_r1+si_l+si_r2, 0.5*d_y+si_r1)
			)
		RVEmodel.sketches['layer_sketch_standard'].Line(
			point1=(cr_l+si_r1+si_l+si_r2, si_r1+0.5*d_y),
			point2=(0, si_r1+0.5*d_y)
			)
		RVEmodel.sketches['layer_sketch_standard'].Line(
			point1=(0, si_r1+0.5*d_y),
			point2=(0, 0)
			)

		# layer sketch with hole
		l_sketch_hole = RVEmodel.ConstrainedSketch(name='layer_sketch_hole', sheetSize=1.0)
		RVEmodel.sketches['layer_sketch_hole'].Line(
			point1=(0, 0),
			point2=(0.5*cr_l, 0)
			)
		genDropletCurve('layer_sketch_hole', 'RIGHT', 0.5*cr_l+si_r1, 0.0, si_r1, si_r2, si_l, False)
		RVEmodel.sketches['layer_sketch_hole'].Line(
			point1=(0.5*cr_l+si_r1+si_l+si_r2, 0),
			point2=(cr_l+si_r1+si_l+si_r2, 0)
			)
		RVEmodel.sketches['layer_sketch_hole'].Line(
			point1=(cr_l+si_r1+si_l+si_r2, 0),
			point2=(cr_l+si_r1+si_l+si_r2, 0.5*d_y+si_r1)
			)
		RVEmodel.sketches['layer_sketch_hole'].Line(
			point1=(cr_l+si_r1+si_l+si_r2, si_r1+0.5*d_y),
			point2=(0, si_r1+0.5*d_y)
			)
		RVEmodel.sketches['layer_sketch_hole'].Line(
			point1=(0, si_r1+0.5*d_y),
			point2=(0, 0)
			)
		RVEmodel.sketches['layer_sketch_hole'].CircleByCenterPerimeter(
			center=(cylndr_x_center, cylndr_y_center), 
			point1=(cylndr_x_center+cylndr_r, cylndr_y_center)
			)

		# soft island sketch
		si_sketch = RVEmodel.ConstrainedSketch(name='si_sketch', sheetSize=1.0)
		genDropletCurve('si_sketch', 'RIGHT', 0.0, 0.0, si_r1, si_r2, si_l, True)

		# crack sketch
		cr_si_sk = RVEmodel.ConstrainedSketch(name='cr_si_sketch', sheetSize=1.0)
		cr_si_sk.rectangle(point1=(0.0, 0.0), point2=(0.5*(si_r1+si_r2+si_l), coh_t))

		cr_layer_sk = RVEmodel.ConstrainedSketch(name='cr_layer_sketch', sheetSize=1.0)
		cr_layer_sk.rectangle(point1=(0.0, 0.0), point2=(0.5*cr_l, coh_t))

		# cylinder sketch
		cylndr_sk = RVEmodel.ConstrainedSketch(name='cylndr_sketch', sheetSize=1.0)
		cylndr_sk.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(cylndr_r, 0.0))


		# PARTS

		layer_name = ['base', 'soft', 'hard', 'interface', 'metal']

		if metal_layer == True:
			layer_t = [d_z, soft_t, hard_t, t_inter, silver_t]
		else:
			layer_t = [d_z, soft_t, hard_t]


		for i in range(len(layer_t)):

			if i in [0, 1]:
				RVEmodel.Part(name='layer_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['layer_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=l_sketch_standard, depth=layer_t[i])
	
				RVEmodel.Part(name='cr_layer_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['cr_layer_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=cr_layer_sk, depth=layer_t[i])
	
				RVEmodel.Part(name='si_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['si_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=si_sketch, depth=layer_t[i])
	
				RVEmodel.Part(name='cr_si_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['cr_si_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=cr_si_sk, depth=layer_t[i])

			if i in [2, 3, 4]:
				if add_secondary_islands == True:
					RVEmodel.Part(name='layer_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
					RVEmodel.parts['layer_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=l_sketch_hole, depth=layer_t[i])

					RVEmodel.Part(name='cylinder_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
					RVEmodel.parts['cylinder_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=cylndr_sk, depth=layer_t[i])
				else:
					RVEmodel.Part(name='layer_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
					RVEmodel.parts['layer_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=l_sketch_standard, depth=layer_t[i])

				RVEmodel.Part(name='cr_layer_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['cr_layer_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=cr_layer_sk, depth=layer_t[i])
	
				RVEmodel.Part(name='si_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['si_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=si_sketch, depth=layer_t[i])
	
				RVEmodel.Part(name='cr_si_part_'+str(layer_name[i]), dimensionality=THREE_D, type=DEFORMABLE_BODY)
				RVEmodel.parts['cr_si_part_'+str(layer_name[i])].BaseSolidExtrude(sketch=cr_si_sk, depth=layer_t[i])

 	except:
 		log('+ + + PARTS ERROR + + +')



# # ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# # /// SECTIONS \\\ #
# # ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	#try:

	# section assignment lists

	solid_soft_list = [
		'layer_part_base',
		'si_part_base',
		'cr_layer_part_base',
		'cr_si_part_base',
		'layer_part_soft',
		'si_part_soft',
		'si_part_hard',
		]
	
	solid_hard_list = [
		'layer_part_hard',
		]

	solid_var_list = [
		]
	
	coh_soft_list = [
		'cr_layer_part_soft',
		'cr_si_part_soft',
		'cr_si_part_hard',
		]
	
	coh_hard_list = [
		'cr_layer_part_hard',
		]
	
	solid_silver_list = [
		'layer_part_metal',
		]

	interface_list = [
		'layer_part_interface',
		]

	if ((add_secondary_islands == True) and (sec_islands_cohesive == False)):
		if metal_layer == False:
			if sec_islands_material == 'polymer_soft':
				solid_soft_list.append('cylinder_part_hard')

			if sec_islands_material == 'polymer_hard':
				solid_hard_list.append('cylinder_part_hard')

			if sec_islands_material == 'polymer_var':
				solid_var_list.append('cylinder_part_hard')

			if sec_islands_material == 'metal':
				solid_silver_list.append('cylinder_part_hard')

		if metal_layer == True:
			if sec_islands_material == 'polymer_soft':
				solid_soft_list.append('cylinder_part_hard')
				solid_soft_list.append('cylinder_part_interface')
				solid_soft_list.append('cylinder_part_metal')

			if sec_islands_material == 'polymer_hard':
				solid_hard_list.append('cylinder_part_hard')
				solid_hard_list.append('cylinder_part_interface')
				solid_hard_list.append('cylinder_part_metal')

			if sec_islands_material == 'polymer_var':
				solid_var_list.append('cylinder_part_hard')
				solid_var_list.append('cylinder_part_interface')
				solid_var_list.append('cylinder_part_metal')

			if sec_islands_material == 'metal':
				solid_silver_list.append('cylinder_part_hard')
				solid_silver_list.append('cylinder_part_interface')
				solid_silver_list.append('cylinder_part_metal')

	if ((add_secondary_islands == True) and (sec_islands_cohesive == True)):
		if metal_layer == False:
			solid_var_list.append('cylinder_part_hard')
		if metal_layer == True:
			interface_list.append('cylinder_part_interface')
			solid_silver_list.append('cylinder_part_metal')
			solid_var_list.append('cylinder_part_hard')


	# soft section
	RVEmodel.HomogeneousSolidSection(
		name='soft_section',
		material='polymer_soft',
		thickness=None)

	for part_name in solid_soft_list:
		part = RVEmodel.parts[part_name]
		part.Set(cells=part.cells, name=part_name+'_set')
		part.SectionAssignment(
			region=part.sets[part_name+'_set'],
			sectionName='soft_section',
			offset=0.0,
			offsetType=MIDDLE_SURFACE,
			offsetField='',
			thicknessAssignment=FROM_SECTION)


	# hard section
	RVEmodel.HomogeneousSolidSection(
		name='hard_section',
		material='polymer_hard',
		thickness=None)

	for part_name in solid_hard_list:
		part = RVEmodel.parts[part_name]
		part.Set(cells=part.cells, name=part_name+'_set')
		part.SectionAssignment(
			region=part.sets[part_name+'_set'],
			sectionName='hard_section',
			offset=0.0,
			offsetType=MIDDLE_SURFACE,
			offsetField='',
			thicknessAssignment=FROM_SECTION)


	# var section
	if add_secondary_islands == True:

		RVEmodel.HomogeneousSolidSection(
			name='var_section',
			material='polymer_var',
			thickness=None)
	
		for part_name in solid_var_list:
			part = RVEmodel.parts[part_name]
			part.Set(cells=part.cells, name=part_name+'_set')
			part.SectionAssignment(
				region=part.sets[part_name+'_set'],
				sectionName='var_section',
				offset=0.0,
				offsetType=MIDDLE_SURFACE,
				offsetField='',
				thicknessAssignment=FROM_SECTION)

	# soft cohesive section
	RVEmodel.CohesiveSection(
		name='soft_coh_section',
		material='polymer_soft_coh',
		response=TRACTION_SEPARATION,
		initialThicknessType=SPECIFY,
		outOfPlaneThickness=coh_t)

	for part_name in coh_soft_list:
		part = RVEmodel.parts[part_name]
		part.Set(cells=part.cells, name=part_name+'_set')
		part.SectionAssignment(
			region=part.sets[part_name+'_set'],
			sectionName='soft_coh_section',
			offset=0.0,
			offsetType=MIDDLE_SURFACE,
			offsetField='',
			thicknessAssignment=FROM_SECTION)
	
	
	# hard cohesive section
	RVEmodel.CohesiveSection(
		name='hard_coh_section',
		material='polymer_hard_coh',
		response=TRACTION_SEPARATION,
		initialThicknessType=SPECIFY,
		outOfPlaneThickness=coh_t)
	
	for part_name in coh_hard_list:
		part = RVEmodel.parts[part_name]
		part.Set(cells=part.cells, name=part_name+'_set')
		part.SectionAssignment(
			region=part.sets[part_name+'_set'],
			sectionName='hard_coh_section',
			offset=0.0,
			offsetType=MIDDLE_SURFACE,
			offsetField='',
			thicknessAssignment=FROM_SECTION)

	
# metal layer

	if metal_layer == True:

		# silver section
		RVEmodel.HomogeneousSolidSection(
			name='metal_section',
			material='metal',
			thickness=None)

		for part_name in solid_silver_list:
			part = RVEmodel.parts[part_name]
			part.Set(cells=part.cells, name=part_name+'_set')
			part.SectionAssignment(
				region=part.sets[part_name+'_set'],
				sectionName='metal_section',
				offset=0.0,
				offsetType=MIDDLE_SURFACE,
				offsetField='',
				thicknessAssignment=FROM_SECTION)

		# interface cohesive section
		RVEmodel.CohesiveSection(
			name='interface_section',
			material='interface_coh',
			response=TRACTION_SEPARATION,
			initialThicknessType=SPECIFY,
			outOfPlaneThickness=coh_t)

		for part_name in interface_list:
			part = RVEmodel.parts[part_name]
			part.Set(cells=part.cells, name=part_name+'_set')
			part.SectionAssignment(
				region=part.sets[part_name+'_set'],
				sectionName='interface_section',
				offset=0.0,
				offsetType=MIDDLE_SURFACE,
				offsetField='',
				thicknessAssignment=FROM_SECTION)


	#except:
	#	log('+ + + SECTIONS ERROR + + +')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// STEP \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	RVEmodel.StaticStep(
		name='Step_1',
		previous='Initial',
		matrixSolver=DIRECT,
		nlgeom=ON,
		maxNumInc=10000,
		initialInc=0.01,
		maxInc=max_inc,
		minInc=1E-009)

	RVEmodel.steps['Step_1'].control.setValues(
	allowPropagation=OFF,
	resetDefaultValues=OFF,
	timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 12.0, 20.0, 6.0, 3.0, 50.0))

	RVEmodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
		'S',		# All stress components.
		'LE',		# All logarithmic strain components.
		'U',		# All physical displacement components.
		'SDEG',		# Scalar stiffness degradation variable d.
		'DMICRT',	# All active components of the damage initiation criteria.
		'STATUS',	# Status of the element (0 or 1).
		'MAXSCRT',	# Maximum nominal stress damage initiation criterion.
		'EVOL',		# Element volume.
		'ENER',		# All energy magnitudes, explanations see below.
		))

	""" ### ENER ### OUTPUT VARIABLES
		# 'DMENER'		# Energy dissipated by damage, per unit volume.
		# 'CENER'		# Energy dissipated by creep, swelling, and viscoelasticity, per unit volume.
		# 'EENER'		# Electrostatic energy density.
		# 'JENER'		# Electrical energy dissipated as a result of the flow of current, per unit volume.
		# 'PENER'		# Energy dissipated by rate-independent and rate-dependent plasticity, per unit volume.
		# 'SENER'		# Elastic strain energy density (with respect to current volume).
		# 'VENER'		# Energy dissipated by viscous effects, per unit volume.
	"""

	mdb.Job(
		name=job_name,
		model=model_name,
		description='',
		type=ANALYSIS,
		atTime=None,
		waitMinutes=0,
		waitHours=0,
		queue=None,
		memory=mem,
		getMemoryFromAnalysis=True,
		explicitPrecision=SINGLE,
		nodalOutputPrecision=SINGLE,
		echoPrint=OFF,
		modelPrint=OFF,
		contactPrint=OFF,
		historyPrint=OFF,
		userSubroutine='',
		scratch='C:\\scratch',
		resultsFormat=ODB,
		multiprocessingMode=DEFAULT,
		)



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// ASSEMBLY & MESH \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	assembly = RVEmodel.rootAssembly

	# layers 
	if metal_layer == True:
		layer_name = ['base', 'soft', 'hard', 'interface', 'metal']
		layer_t = [d_z, soft_t, hard_t, t_inter, silver_t]
	else:
		layer_name = ['base', 'soft', 'hard']
		layer_t = [d_z, soft_t, hard_t]

	# assemble unit cell
	translate = 0.0
	for i in range(len(layer_t)):
		assembly.Instance(name='inst_layer_part_'+str(layer_name[i]), part=RVEmodel.parts['layer_part_'+str(layer_name[i])], dependent=ON)
		assembly.translate(instanceList=('inst_layer_part_'+str(layer_name[i]), ), vector=(0.0, 0.0, translate))
		if ((add_secondary_islands == True) and (i > 1)):
			assembly.Instance(name='inst_cylinder_part_'+str(layer_name[i]), part=RVEmodel.parts['cylinder_part_'+str(layer_name[i])], dependent=ON)
			assembly.translate(instanceList=('inst_cylinder_part_'+str(layer_name[i]), ), vector=(cylndr_x_center, cylndr_y_center, translate))
		else:
			pass
		if i <= 2:
			assembly.Instance(name='inst_si_part_'+str(layer_name[i]), part=RVEmodel.parts['si_part_'+str(layer_name[i])], dependent=ON)
			assembly.translate(instanceList=('inst_si_part_'+str(layer_name[i]), ), vector=(0.0, 0.0, translate))
			assembly.translate(instanceList=('inst_si_part_'+str(layer_name[i]), ), vector=(si_r1+0.5*cr_l, 0.0, 0.0))
		else:
			pass

		translate = translate + layer_t[i]

	# merge unit cell
	instances = []
	for key in assembly.instances.keys():
		instances.append(assembly.instances[key])
	inst_tup = tuple(instances)
	unit_cell_inst = assembly.InstanceFromBooleanMerge(name='unit_cell',
		instances=inst_tup,
		keepIntersections=True, originalInstances=DELETE, domain=BOTH)
	unit_cell_part = RVEmodel.parts['unit_cell']


	# assemble and cut unit cell
	cutter_sk = RVEmodel.ConstrainedSketch(name='cutter_sketch', sheetSize=1.0)
	cutter_sk.rectangle(point1=(0.0, 0.0), point2=(x_l, si_r1+d_y))
	RVEmodel.Part(name='cutter_part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	RVEmodel.parts['cutter_part'].BaseSolidExtrude(sketch=cutter_sk, depth=sum(layer_t))
	cutter_part = RVEmodel.parts['cutter_part']

	unit_cell_cutter = assembly.Instance(name='unit_cell-2', part=cutter_part, dependent=ON)
	
	assembly.translate(instanceList=('unit_cell-2', ), vector=(0.5*x_l, 0.0, 0.0))
	unit_cell_blunt = assembly.InstanceFromBooleanCut(
		cuttingInstances=(unit_cell_cutter, ),
		instanceToBeCut=unit_cell_inst, name='unit_cell_blunt', originalInstances=DELETE)

	unit_cell_cut = assembly.Instance(name='unit_cell-3', part=unit_cell_part, dependent=ON)
	
	unit_cell_cutter = assembly.Instance(name='unit_cell-4', part=cutter_part, dependent=ON)
	assembly.translate(instanceList=('unit_cell-4', ), vector=(-0.5*x_l, 0.0, 0.0))
	unit_cell_sharp = assembly.InstanceFromBooleanCut(
		cuttingInstances=(unit_cell_cutter, ),
		instanceToBeCut=unit_cell_cut, name='unit_cell_sharp', originalInstances=DELETE)


	# delete instances
	assembly.deleteFeatures(('unit_cell_blunt-1', 'unit_cell_sharp-1', ))


	# assemble cohesive cell layer
	translate = 0.0
	if metal_layer == True:
		layer_t = [d_z, soft_t, ]
	else:
		layer_t = [d_z, soft_t, hard_t]
	
	for i in range(len(layer_t)):
		assembly.Instance(name='inst_cr_layer_part_'+str(layer_name[i]), part=RVEmodel.parts['cr_layer_part_'+str(layer_name[i])], dependent=ON)
		assembly.translate(instanceList=('inst_cr_layer_part_'+str(layer_name[i]), ), vector=(0.0, 0.0, translate))
		translate = translate + layer_t[i]


	# merge coh cell layer
	instances = []
	for key in assembly.instances.keys():
		instances.append(assembly.instances[key])
	inst_tup = tuple(instances)
	coh_cell_layer_inst = assembly.InstanceFromBooleanMerge(name='coh_cell_layer',
		instances=inst_tup,
		keepIntersections=True, originalInstances=DELETE, domain=GEOMETRY)
	assembly.deleteFeatures(('coh_cell_layer-1', ))


	# assemble cohesive cell si
	translate = 0.0
	layer_t = [d_z, soft_t, hard_t]
	for i in range(len(layer_t)):
		assembly.Instance(name='inst_cr_si_part_'+str(layer_name[i]), part=RVEmodel.parts['cr_si_part_'+str(layer_name[i])], dependent=ON)
		assembly.translate(instanceList=('inst_cr_si_part_'+str(layer_name[i]), ), vector=(0.0, 0.0, translate))
		translate = translate + layer_t[i]

	layer_t = [d_z, soft_t, hard_t, silver_t, t_inter]

	# merge coh cell si
	instances = []
	for key in assembly.instances.keys():
		instances.append(assembly.instances[key])
	inst_tup = tuple(instances)
	coh_cell_si_inst = assembly.InstanceFromBooleanMerge(name='coh_cell_si',
		instances=inst_tup,
		keepIntersections=True, originalInstances=DELETE, domain=GEOMETRY)
	assembly.deleteFeatures(('coh_cell_si-1', ))


	# mesh parts
	unit_cell_blunt_part = RVEmodel.parts['unit_cell_blunt']
	unit_cell_sharp_part = RVEmodel.parts['unit_cell_sharp']
	coh_cell_layer_part = RVEmodel.parts['coh_cell_layer']
	coh_cell_si_part = RVEmodel.parts['coh_cell_si']
	start_mesh = time()


	# create instances
	assembly.Instance(name='inst_unit_cell_blunt', part=unit_cell_blunt_part, dependent=ON)
	assembly.Instance(name='inst_unit_cell_sharp', part=unit_cell_sharp_part, dependent=ON)
	assembly.Instance(name='inst_coh_cell_layer', part=coh_cell_layer_part, dependent=ON)
	assembly.Instance(name='inst_coh_cell_si', part=coh_cell_si_part, dependent=ON)


	# set to sweep mesh
	unit_cell_blunt_part.setMeshControls(regions=assembly.instances['inst_unit_cell_blunt'].cells, technique=SWEEP)
	unit_cell_sharp_part.setMeshControls(regions=assembly.instances['inst_unit_cell_sharp'].cells, technique=SWEEP)
	coh_cell_layer_part.setMeshControls(regions=assembly.instances['inst_coh_cell_layer'].cells, technique=SWEEP)
	coh_cell_si_part.setMeshControls(regions=assembly.instances['inst_coh_cell_si'].cells, technique=SWEEP)



	# seed parts
	for part in [unit_cell_blunt_part, unit_cell_sharp_part, coh_cell_layer_part, coh_cell_si_part]:
		# seed d_z
		dz_edges_seq = []
		for edge in part.edges:
			try:
				if round(edge.getSize(printResults = False), 3) == round(d_z, 3):
					vert = edge.getVertices()
					if part.vertices[vert[0]].pointOn[0][2] != part.vertices[vert[1]].pointOn[0][2]:
						dz_edges_seq.append(part.edges[edge.index])
			except:
				None
		part.seedEdgeByNumber(edges=dz_edges_seq, number=elnum_dz, constraint=FINER)

		# seed soft_t
		so_edges_seq = []
		for edge in part.edges:
			try:
				if round(edge.getSize(printResults = False), 3) == round(soft_t, 3):
					vert = edge.getVertices()
					if part.vertices[vert[0]].pointOn[0][2] != part.vertices[vert[1]].pointOn[0][2]:
						so_edges_seq.append(part.edges[edge.index])
			except:
				None
		for edge in so_edges_seq:
			verti = edge.getVertices()
			ver_z0 = part.vertices[verti[0]].pointOn[0][2]
			ver_z1 = part.vertices[verti[1]].pointOn[0][2]
			if ver_z0 == d_z+soft_t:
				part.seedEdgeByNumber(edges=so_edges_seq, number=elnum_so, constraint=FINER)
			elif ver_z0 == d_z:
				part.seedEdgeByNumber(edges=so_edges_seq, number=elnum_so, constraint=FINER)

		# seed hard_t
		ha_edges_seq = []
		for edge in part.edges:
			try:
				if round(edge.getSize(printResults = False), 3) == round(hard_t, 3):
					ha_edges_seq.append(part.edges[edge.index])
			except:
				None
		part.seedEdgeByNumber(edges=ha_edges_seq, number=elnum_ha, constraint=FINER)

		# seed coh_t
		coh_t_edges_seq = []
		for edge in part.edges:
			try:
				if round(edge.getSize(printResults = False), 3) in (round(coh_t, 3),):
					coh_t_edges_seq.append(part.edges[edge.index])
			except:
				None
		part.seedEdgeBySize(edges=coh_t_edges_seq, size=coh_t, constraint=FINER)

		# global seed and generate mesh
		part.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)
		try:
			part.generateMesh()
		except:
			log('+ + + MESH ERROR + + +')


	# remove from assembly to start fresh with orphan meshes
	assembly.deleteFeatures(('inst_unit_cell_blunt', 'inst_unit_cell_sharp', 'inst_coh_cell_layer', 'inst_coh_cell_si', ))


	# generate orphan mesh parts
	uc_blunt_1_part = unit_cell_blunt_part.PartFromMesh(name='unit_cell_blunt_mesh_1', copySets=True)
	uc_sharp_1_part = unit_cell_sharp_part.PartFromMesh(name='unit_cell_sharp_mesh_1', copySets=True)


	# mirror orphan mesh parts
	uc_blunt_2_part = RVEmodel.Part(name='unit_cell_blunt_mesh_2', objectToCopy=uc_blunt_1_part, compressFeatureList=ON, mirrorPlane=XZPLANE)
	uc_sharp_2_part = RVEmodel.Part(name='unit_cell_sharp_mesh_2', objectToCopy=uc_sharp_1_part, compressFeatureList=ON, mirrorPlane=XZPLANE)


	# assemble RVE
	assembly.Instance(name='unit_cell_blunt_1-1', part=uc_blunt_1_part, dependent=ON)
	assembly.Instance(name='unit_cell_sharp_1-2', part=uc_sharp_1_part, dependent=ON)
	assembly.Instance(name='unit_cell_blunt_1-3', part=uc_blunt_1_part, dependent=ON)
	assembly.Instance(name='unit_cell_sharp_1-4', part=uc_sharp_1_part, dependent=ON)
	assembly.Instance(name='unit_cell_blunt_2-5', part=uc_blunt_2_part, dependent=ON)
	assembly.Instance(name='unit_cell_sharp_2-6', part=uc_sharp_2_part, dependent=ON)
	assembly.Instance(name='unit_cell_blunt_2-7', part=uc_blunt_2_part, dependent=ON)
	assembly.Instance(name='unit_cell_sharp_2-8', part=uc_sharp_2_part, dependent=ON)

	assembly.translate(instanceList=('unit_cell_blunt_1-1', 'unit_cell_sharp_1-2',), vector=(0.0, coh_t, 0.0))
	assembly.translate(instanceList=('unit_cell_blunt_2-5', ), vector=(0.5*x_l, 0.5*y_l+coh_t, 0.0))
	assembly.translate(instanceList=('unit_cell_sharp_2-6', ), vector=(-0.5*x_l, 0.5*y_l+coh_t, 0.0))
	assembly.translate(instanceList=('unit_cell_blunt_1-3', ), vector=(0.5*x_l, 0.5*y_l+2*coh_t, 0.0))
	assembly.translate(instanceList=('unit_cell_sharp_1-4', ), vector=(-0.5*x_l, 0.5*y_l+2*coh_t, 0.0))

	assembly.Instance(name='coh_cell_layer_part-1', part=coh_cell_layer_part, dependent=ON)
	assembly.Instance(name='coh_cell_layer_part-2', part=coh_cell_layer_part, dependent=ON)
	assembly.Instance(name='coh_cell_layer_part-3', part=coh_cell_layer_part, dependent=ON)
	assembly.Instance(name='coh_cell_layer_part-4', part=coh_cell_layer_part, dependent=ON)
	assembly.Instance(name='coh_cell_si_part-1', part=coh_cell_si_part, dependent=ON)
	assembly.Instance(name='coh_cell_si_part-2', part=coh_cell_si_part, dependent=ON)
	assembly.Instance(name='coh_cell_si_part-3', part=coh_cell_si_part, dependent=ON)
	assembly.Instance(name='coh_cell_si_part-4', part=coh_cell_si_part, dependent=ON)

	assembly.translate(instanceList=('coh_cell_layer_part-2',), vector=(0.5*cr_l+si_r1+si_r2+si_l, 0.0, 0.0))
	assembly.translate(instanceList=('coh_cell_layer_part-3',), vector=(0.5*(si_r1+si_r2+si_l), 0.5*y_l+coh_t, 0.0))
	assembly.translate(instanceList=('coh_cell_layer_part-4',), vector=(0.5*(cr_l+si_r1+si_r2+si_l), 0.5*y_l+coh_t, 0.0))

	assembly.translate(instanceList=('coh_cell_si_part-1',), vector=(0.5*cr_l, 0.0, 0.0))
	assembly.translate(instanceList=('coh_cell_si_part-2',), vector=(0.5*cr_l+0.5*(si_r1+si_r2+si_l), 0.0, 0.0))
	assembly.translate(instanceList=('coh_cell_si_part-3',), vector=(0.0, 0.5*y_l+coh_t, 0.0))
	assembly.translate(instanceList=('coh_cell_si_part-4',), vector=(cr_l+0.5*(si_r1+si_r2+si_l), 0.5*y_l+coh_t, 0.0))


	# merge all instances
	instances = []
	for key in assembly.instances.keys():
		instances.append(assembly.instances[key])
	inst_tup = tuple(instances)
	rve_inst = assembly.InstanceFromBooleanMerge(
		name='RVE',
		instances=inst_tup,
		keepIntersections=True,
		mergeNodes=BOUNDARY_ONLY,
		nodeMergingTolerance=1e-05,
		originalInstances=DELETE,
		domain=MESH)
	rve_part = RVEmodel.parts['RVE']


	# assign full integration elements
	elem_full = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
	full_regions = []
	full_list = [
		'layer_part_hard',
		'si_part_hard',
		'layer_part_soft',
		'si_part_soft',
		]
	
	if add_secondary_islands == True:
		if metal_layer == True:
			full_list.append('layer_part_metal')
			if sec_islands_cohesive == True:
				full_list.append('cylinder_part_metal')
				full_list.append('cylinder_part_hard')
			if sec_islands_cohesive == False:
				full_list.append('cylinder_part_metal')
				full_list.append('cylinder_part_hard')
				full_list.append('cylinder_part_interface')
		if metal_layer == False:
			full_list.append('cylinder_part_hard')
	if add_secondary_islands == False:
		if metal_layer == True:
			full_list.append('layer_part_metal')


	for part in full_list:
		full_regions.append(rve_part.sets[part+'_set'].elements)
	full_tuple = tuple(full_regions)
	rve_part.setElementType(regions=full_tuple, elemTypes=(elem_full, ))


	# assign cohesive elements
	elem_coh = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD, viscosity=visc)
	coh_regions = []
	
	if metal_layer == True:
		for part in interface_list+coh_soft_list:
			coh_regions.append(rve_part.sets[part+'_set'].elements)
	else:
		for part in coh_soft_list+coh_hard_list:
			coh_regions.append(rve_part.sets[part+'_set'].elements)
	coh_tuple = tuple(coh_regions)
	rve_part.setElementType(regions=coh_tuple, elemTypes=(elem_coh, ))


	# assign stack direction to cohesive elements in crack
	bbox_e = rve_part.elements.getByBoundingBox(
		-0.4*seed_size,
		-0.4*seed_size-0.5*d_y-si_r1,
		-0.4*seed_size,
		1.4*seed_size,
		1.4*seed_size-0.5*d_y-si_r1,
		1.4*d_z*elnum_dz**(-1)
		)

	bbox_f = bbox_e[0].getElemFaces()

	for f in bbox_f:
		bbox_n = f.getNormal()
		if bbox_n == (0.0, -1.0, 0.0):
			bbox_i = f

	if metal_layer == True:
		for set in coh_soft_list:
			rve_part.orientElements(referenceRegion=bbox_i, pickedElements=rve_part.sets[set+'_set'].elements)
	else:
		for set in coh_soft_list+coh_hard_list:
			rve_part.orientElements(referenceRegion=bbox_i, pickedElements=rve_part.sets[set+'_set'].elements)
	

	log('Mesh time: '+str(round(time()-start_mesh, 1))+' seconds ('+str(len(rve_part.elements))+' elements).')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// PERIODIC BOUNDARY CONDITIONS \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	E33					=	0.0							# z strain
	E13					=	0.0							# xz strain
	E23					=	0.0							# yz strain

	try:
		# PBC NODE BASED

		# setup
		start_pbc = time()
		name_instance = 'RVE-1'
		name_step = 'Step_1'
		tol = 0.5*coh_t
		tol_digits = 7

		x_min = 0.0
		x_max = 0.0
		y_min = 0.0
		y_max = 0.0
		z_min = 0.0
		z_max = 0.0


		# get min/max coordinates of RVE
		for node in assembly.instances[name_instance].nodes:
			if (node.coordinates[0]<x_min):
				x_min = node.coordinates[0]
			if (node.coordinates[0]>x_max):
				x_max = node.coordinates[0]
			if (node.coordinates[1]<y_min):
				y_min = node.coordinates[1]
			if (node.coordinates[1]>y_max):
				y_max = node.coordinates[1]
			if (node.coordinates[2]<z_min):
				z_min = node.coordinates[2]
			if (node.coordinates[2]>z_max):
				z_max = node.coordinates[2]

		edge_length_x = x_max-x_min
		edge_length_y = y_max-y_min
		edge_length_z = z_max-z_min


		# node sequences
		x_neg_set_face_seq = []
		y_neg_set_face_seq = []
		z_neg_set_face_seq = []
		x_pos_set_face_seq = []
		y_pos_set_face_seq = []
		z_pos_set_face_seq = []

		set_edge_BF_seq = []
		set_edge_CG_seq = []
		set_edge_BC_seq = []
		set_edge_FG_seq = []
		set_edge_AE_seq = []
		set_edge_DH_seq = []
		set_edge_AD_seq = []
		set_edge_EH_seq = []
		set_edge_CD_seq = []
		set_edge_GH_seq = []
		set_edge_AB_seq = []
		set_edge_EF_seq = []

		set_vertex_A_seq = []
		set_vertex_B_seq = []
		set_vertex_C_seq = []
		set_vertex_D_seq = []
		set_vertex_E_seq = []
		set_vertex_F_seq = []
		set_vertex_G_seq = []
		set_vertex_H_seq = []


		# get node sequences
		for node in assembly.instances[name_instance].nodes:

			n0 = node.coordinates[0]
			n1 = node.coordinates[1]
			n2 = node.coordinates[2]

			n0r = round(n0, tol_digits)
			n1r = round(n1, tol_digits)
			n2r = round(n2, tol_digits)

			x_minr = round(x_min, tol_digits)
			y_minr = round(y_min, tol_digits)
			z_minr = round(z_min, tol_digits)
			x_maxr = round(x_max, tol_digits)
			y_maxr = round(y_max, tol_digits)
			z_maxr = round(z_max, tol_digits)

			# face x_neg seq
			if	(		n0r == x_minr
					and n1r != y_minr
					and n1r != y_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				x_neg_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# face y_neg seq
			if	(		n1r == y_minr
					and n0r != x_minr
					and n0r != x_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				y_neg_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# face z_neg seq
			if	(		n2r == z_minr
					and n0r != x_minr
					and n0r != x_maxr
					and n1r != y_minr
					and n1r != y_maxr
				):
				z_neg_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# face x_pos seq
			if	(		n0r == x_maxr
					and n1r != y_minr
					and n1r != y_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				x_pos_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# face y_pos seq
			if	(		n1r == y_maxr
					and n0r != x_minr
					and n0r != x_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				y_pos_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# face z_pos seq
			if	(		n2r == z_maxr
					and n0r != x_minr
					and n0r != x_maxr
					and n1r != y_minr
					and n1r != y_maxr
				):
				z_pos_set_face_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge AE seq
			if	(		n0r == x_minr
					and n1r == y_minr
					and n2r != z_minr
					and n2r != z_maxr
				):
				set_edge_AE_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge BF seq
			if	(		n0r == x_minr
					and n1r == y_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				set_edge_BF_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge CG seq
			if	(		n0r == x_maxr
					and n1r == y_maxr
					and n2r != z_minr
					and n2r != z_maxr
				):
				set_edge_CG_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge DH seq
			if	(		n0r == x_maxr
					and n1r == y_minr
					and n2r != z_minr
					and n2r != z_maxr
				):
				set_edge_DH_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge BC seq
			if	(		n1r == y_maxr
					and n2r == z_minr
					and n0r != x_minr
					and n0r != x_maxr
				):
				set_edge_BC_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge FG seq
			if	(		n1r == y_maxr
					and n2r == z_maxr
					and n0r != x_minr
					and n0r != x_maxr
				):
				set_edge_FG_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge EH seq
			if	(		n1r == y_minr
					and n2r == z_maxr
					and n0r != x_minr
					and n0r != x_maxr
				):
				set_edge_EH_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge AD seq
			if	(		n1r == y_minr
					and n2r == z_minr
					and n0r != x_minr
					and n0r != x_maxr
				):
				set_edge_AD_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge CD seq
			if	(		n0r == x_maxr
					and n2r == z_minr
					and n1r != y_minr
					and n1r != y_maxr
				):
				set_edge_CD_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge GH seq
			if	(		n0r == x_maxr
					and n2r == z_maxr
					and n1r != y_minr
					and n1r != y_maxr
				):
				set_edge_GH_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge AB seq
			if	(		n0r == x_minr
					and n2r == z_minr
					and n1r != y_minr
					and n1r != y_maxr
				):
				set_edge_AB_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# edge EF seq
			if	(		n0r == x_minr
					and n2r == z_maxr
					and n1r != y_minr
					and n1r != y_maxr
				):
				set_edge_EF_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex A seq
			if	(		n0r == x_minr
					and n1r == y_minr
					and n2r == z_minr
				):
				set_vertex_A_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex B seq
			if	(		n0r == x_minr
					and n1r == y_maxr
					and n2r == z_minr
				):
				set_vertex_B_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex C seq
			if	(		n0r == x_maxr
					and n1r == y_maxr
					and n2r == z_minr
				):
				set_vertex_C_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex D seq
			if	(		n0r == x_maxr
					and n1r == y_minr
					and n2r == z_minr
				):
				set_vertex_D_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex E seq
			if	(		n0r == x_minr
					and n1r == y_minr
					and n2r == z_maxr
				):
				set_vertex_E_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex F seq
			if	(		n0r == x_minr
					and n1r == y_maxr
					and n2r == z_maxr
				):
				set_vertex_F_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex G seq
			if	(		n0r == x_maxr
					and n1r == y_maxr
					and n2r == z_maxr
				):
				set_vertex_G_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])

			# vertex H seq
			if	(		n0r == x_maxr
					and n1r == y_minr
					and n2r == z_maxr
				):
				set_vertex_H_seq.append(assembly.instances[name_instance].nodes[node.label-1:node.label])


		# create face, edge, vertex sets
		x_neg_set_face_nodes = assembly.Set(nodes=x_neg_set_face_seq, name='NODES_FACE_X_NEG')
		y_neg_set_face_nodes = assembly.Set(nodes=y_neg_set_face_seq, name='NODES_FACE_Y_NEG')
		z_neg_set_face_nodes = assembly.Set(nodes=z_neg_set_face_seq, name='NODES_FACE_Z_NEG')
		x_pos_set_face_nodes = assembly.Set(nodes=x_pos_set_face_seq, name='NODES_FACE_X_POS')
		y_pos_set_face_nodes = assembly.Set(nodes=y_pos_set_face_seq, name='NODES_FACE_Y_POS')
		z_pos_set_face_nodes = assembly.Set(nodes=z_pos_set_face_seq, name='NODES_FACE_Z_POS')

		set_edge_BF_nodes = assembly.Set(nodes=set_edge_BF_seq, name='NODES_EDGES_BF')
		set_edge_CG_nodes = assembly.Set(nodes=set_edge_CG_seq, name='NODES_EDGES_CG')
		set_edge_BC_nodes = assembly.Set(nodes=set_edge_BC_seq, name='NODES_EDGES_BC')
		set_edge_FG_nodes = assembly.Set(nodes=set_edge_FG_seq, name='NODES_EDGES_FG')
		set_edge_AE_nodes = assembly.Set(nodes=set_edge_AE_seq, name='NODES_EDGES_AE')
		set_edge_DH_nodes = assembly.Set(nodes=set_edge_DH_seq, name='NODES_EDGES_DH')
		set_edge_AD_nodes = assembly.Set(nodes=set_edge_AD_seq, name='NODES_EDGES_AD')
		set_edge_EH_nodes = assembly.Set(nodes=set_edge_EH_seq, name='NODES_EDGES_EH')
		set_edge_CD_nodes = assembly.Set(nodes=set_edge_CD_seq, name='NODES_EDGES_CD')
		set_edge_GH_nodes = assembly.Set(nodes=set_edge_GH_seq, name='NODES_EDGES_GH')
		set_edge_AB_nodes = assembly.Set(nodes=set_edge_AB_seq, name='NODES_EDGES_AB')
		set_edge_EF_nodes = assembly.Set(nodes=set_edge_EF_seq, name='NODES_EDGES_EF')

		set_vertex_A_nodes = assembly.Set(nodes=set_vertex_A_seq, name='NODE_VERTEX_A')
		set_vertex_B_nodes = assembly.Set(nodes=set_vertex_B_seq, name='NODE_VERTEX_B')
		set_vertex_C_nodes = assembly.Set(nodes=set_vertex_C_seq, name='NODE_VERTEX_C')
		set_vertex_D_nodes = assembly.Set(nodes=set_vertex_D_seq, name='NODE_VERTEX_D')
		set_vertex_E_nodes = assembly.Set(nodes=set_vertex_E_seq, name='NODE_VERTEX_E')
		set_vertex_F_nodes = assembly.Set(nodes=set_vertex_F_seq, name='NODE_VERTEX_F')
		set_vertex_G_nodes = assembly.Set(nodes=set_vertex_G_seq, name='NODE_VERTEX_G')
		set_vertex_H_nodes = assembly.Set(nodes=set_vertex_H_seq, name='NODE_VERTEX_H')


		# reference points
		ref_point_normal = assembly.ReferencePoint(point=(0.0, 0.0, 0.0)).id
		ref_point_shear1 = assembly.ReferencePoint(point=(10.0, 0.0, 0.0)).id
		ref_point_shear2 = assembly.ReferencePoint(point=(20.0, 0.0, 0.0)).id
		ref_point_shear3 = assembly.ReferencePoint(point=(30.0, 0.0, 0.0)).id
		assembly.Set(referencePoints=(assembly.referencePoints[ref_point_normal],), name='Ref_Point_normal')
		assembly.Set(referencePoints=(assembly.referencePoints[ref_point_shear1],), name='Ref_Point_shear1')
		assembly.Set(referencePoints=(assembly.referencePoints[ref_point_shear2],), name='Ref_Point_shear2')
		assembly.Set(referencePoints=(assembly.referencePoints[ref_point_shear3],), name='Ref_Point_shear3')
		ref_point_normal = assembly.sets['Ref_Point_normal']
		ref_point_shear1 = assembly.sets['Ref_Point_shear1']
		ref_point_shear2 = assembly.sets['Ref_Point_shear2']
		ref_point_shear3 = assembly.sets['Ref_Point_shear3']


		# apply boundary conditions to reference points
		RVEmodel.DisplacementBC(name='BC_normal', createStepName='Step_1', region=ref_point_normal,
			u1=E11*edge_length_x, u2=E22*edge_length_y, u3=E33*edge_length_z, ur1=SET, ur2=SET, ur3=SET,
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		RVEmodel.DisplacementBC(name='BC_shear1', createStepName='Step_1', region=ref_point_shear1,
			u1=E12*edge_length_y, u2=E12*edge_length_x, u3=SET, ur1=SET, ur2=SET, ur3=SET,
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		RVEmodel.DisplacementBC(name='BC_shear2', createStepName='Step_1', region=ref_point_shear2,
			u1=E13*edge_length_z, u2=E13*edge_length_x, u3=SET, ur1=SET, ur2=SET, ur3=SET,
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		RVEmodel.DisplacementBC(name='BC_shear3', createStepName='Step_1', region=ref_point_shear3,
			u1=E23*edge_length_y, u2=E23*edge_length_z, u3=SET, ur1=SET, ur2=SET, ur3=SET,
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


		# PBC FACES X
		for i in range(len(x_pos_set_face_nodes.nodes)):
			x = x_pos_set_face_nodes.nodes[i].coordinates[0]
			y = x_pos_set_face_nodes.nodes[i].coordinates[1]
			z = x_pos_set_face_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(x_pos_set_face_nodes.nodes[i].label,)),),
				name='NODE_FACE_X_POS'+str(i))
			opposite_node=x_neg_set_face_nodes.nodes.getByBoundingBox(x_min-tol,y-tol,z-tol,x_min+tol,y+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_FACE_X_NEG'+str(i))
			RVEmodel.Equation(
				name='Equ_face_x'+str(i)+str(0),
				terms=((1.0,'NODE_FACE_X_POS'+str(i),1),(-1.0,'NODE_FACE_X_NEG'+str(i),1),(-1.0,'Ref_Point_normal',1)))
			RVEmodel.Equation(
				name='Equ_face_x'+str(i)+str(1),
				terms=((1.0,'NODE_FACE_X_POS'+str(i),2),(-1.0,'NODE_FACE_X_NEG'+str(i),2),(-1.0,'Ref_Point_shear1',2)))
			RVEmodel.Equation(
				name='Equ_face_x'+str(i)+str(2),
				terms=((1.0,'NODE_FACE_X_POS'+str(i),3),(-1.0,'NODE_FACE_X_NEG'+str(i),3),(-1.0,'Ref_Point_shear2',2)))

		# PBC FACES Y
		for i in range(len(y_pos_set_face_nodes.nodes)):
			x = y_pos_set_face_nodes.nodes[i].coordinates[0]
			y = y_pos_set_face_nodes.nodes[i].coordinates[1]
			z = y_pos_set_face_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(y_pos_set_face_nodes.nodes[i].label,)),),
				name='NODE_FACE_Y_POS'+str(i))
			opposite_node=y_neg_set_face_nodes.nodes.getByBoundingBox(x-tol,y_min-tol,z-tol,x+tol,y_min+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_FACE_Y_NEG'+str(i))
			RVEmodel.Equation(
				name='Equ_face_y'+str(i)+str(0),
				terms=((1.0,'NODE_FACE_Y_POS'+str(i),1),(-1.0,'NODE_FACE_Y_NEG'+str(i),1),(-1.0,'Ref_Point_shear1',1)))
			RVEmodel.Equation(
				name='Equ_face_y'+str(i)+str(1),
				terms=((1.0,'NODE_FACE_Y_POS'+str(i),2),(-1.0,'NODE_FACE_Y_NEG'+str(i),2),(-1.0,'Ref_Point_normal',2)))
			RVEmodel.Equation(
				name='Equ_face_y'+str(i)+str(2),
				terms=((1.0,'NODE_FACE_Y_POS'+str(i),3),(-1.0,'NODE_FACE_Y_NEG'+str(i),3),(-1.0,'Ref_Point_shear3',1)))

		# PBC EDGES DH CG
		for i in range(len(set_edge_DH_nodes.nodes)):
			x = set_edge_DH_nodes.nodes[i].coordinates[0]
			y = set_edge_DH_nodes.nodes[i].coordinates[1]
			z = set_edge_DH_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_DH_nodes.nodes[i].label,)),),
				name='NODE_EDGE_DH_1'+str(i))
			opposite_node=set_edge_CG_nodes.nodes.getByBoundingBox(x-tol,y_max-tol,z-tol,x+tol,y_max+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_CG'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_DH_CG'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_DH_1'+str(i),1),(-1.0,'NODE_EDGE_CG'+str(i),1),(1.0,'Ref_Point_shear1',1)))
			RVEmodel.Equation(
				name='Equ_edge_DH_CG'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_DH_1'+str(i),2),(-1.0,'NODE_EDGE_CG'+str(i),2),(1.0,'Ref_Point_normal',2)))
			RVEmodel.Equation(
				name='Equ_edge_DH_CG'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_DH_1'+str(i),3),(-1.0,'NODE_EDGE_CG'+str(i),3),(1.0,'Ref_Point_shear3',1)))

		# PBC EDGES BF AE
		for i in range(len(set_edge_BF_nodes.nodes)):
			x = set_edge_BF_nodes.nodes[i].coordinates[0]
			y = set_edge_BF_nodes.nodes[i].coordinates[1]
			z = set_edge_BF_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_BF_nodes.nodes[i].label,)),),
				name='NODE_EDGE_BF'+str(i))
			opposite_node=set_edge_AE_nodes.nodes.getByBoundingBox(x-tol,y_min-tol,z-tol,x+tol,y_min+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_AE_1'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_BF_AE'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_BF'+str(i),1),(-1.0,'NODE_EDGE_AE_1'+str(i),1),(-1.0,'Ref_Point_shear1',1)))
			RVEmodel.Equation(
				name='Equ_edge_BF_AE'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_BF'+str(i),2),(-1.0,'NODE_EDGE_AE_1'+str(i),2),(-1.0,'Ref_Point_normal',2)))
			RVEmodel.Equation(
				name='Equ_edge_BF_AE'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_BF'+str(i),3),(-1.0,'NODE_EDGE_AE_1'+str(i),3),(-1.0,'Ref_Point_shear3',1)))

		# PBC EDGES AE DH
		for i in range(len(set_edge_AE_nodes.nodes)):
			x = set_edge_AE_nodes.nodes[i].coordinates[0]
			y = set_edge_AE_nodes.nodes[i].coordinates[1]
			z = set_edge_AE_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_AE_nodes.nodes[i].label,)),),
				name='NODE_EDGE_AE_2'+str(i))
			opposite_node=set_edge_DH_nodes.nodes.getByBoundingBox(x_max-tol,y-tol,z-tol,x_max+tol,y+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_DH_2'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_AE_DH'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_AE_2'+str(i),1),(-1.0,'NODE_EDGE_DH_2'+str(i),1),(1.0,'Ref_Point_normal',1)))
			RVEmodel.Equation(
				name='Equ_edge_AE_DH'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_AE_2'+str(i),2),(-1.0,'NODE_EDGE_DH_2'+str(i),2),(1.0,'Ref_Point_shear1',2)))
			RVEmodel.Equation(
				name='Equ_edge_AE_DH'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_AE_2'+str(i),3),(-1.0,'NODE_EDGE_DH_2'+str(i),3),(1.0,'Ref_Point_shear2',2)))

		# PBC EDGES GH EF
		for i in range(len(set_edge_GH_nodes.nodes)):
			x = set_edge_GH_nodes.nodes[i].coordinates[0]
			y = set_edge_GH_nodes.nodes[i].coordinates[1]
			z = set_edge_GH_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_GH_nodes.nodes[i].label,)),),
				name='NODE_EDGE_GH'+str(i))
			opposite_node=set_edge_EF_nodes.nodes.getByBoundingBox(x_min-tol,y-tol,z-tol,x_min+tol,y+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_EF_1'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_GH_EF'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_GH'+str(i),1),(-1.0,'NODE_EDGE_EF_1'+str(i),1),(-1.0,'Ref_Point_normal',1)))
			RVEmodel.Equation(
				name='Equ_edge_GH_EF'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_GH'+str(i),2),(-1.0,'NODE_EDGE_EF_1'+str(i),2),(-1.0,'Ref_Point_shear1',2)))
			RVEmodel.Equation(
				name='Equ_edge_GH_EF'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_GH'+str(i),3),(-1.0,'NODE_EDGE_EF_1'+str(i),3),(-1.0,'Ref_Point_shear2',2)))

		# PBC EDGES AB CD
		for i in range(len(set_edge_AB_nodes.nodes)):
			x = set_edge_AB_nodes.nodes[i].coordinates[0]
			y = set_edge_AB_nodes.nodes[i].coordinates[1]
			z = set_edge_AB_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_AB_nodes.nodes[i].label,)),),
				name='NODE_EDGE_AB_2'+str(i))
			opposite_node=set_edge_CD_nodes.nodes.getByBoundingBox(x_max-tol,y-tol,z-tol,x_max+tol,y+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_CD'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_AB_CD'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_AB_2'+str(i),1),(-1.0,'NODE_EDGE_CD'+str(i),1),(1.0,'Ref_Point_normal',1)))
			RVEmodel.Equation(
				name='Equ_edge_AB_CD'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_AB_2'+str(i),2),(-1.0,'NODE_EDGE_CD'+str(i),2),(1.0,'Ref_Point_shear1',2)))
			RVEmodel.Equation(
				name='Equ_edge_AB_CD'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_AB_2'+str(i),3),(-1.0,'NODE_EDGE_CD'+str(i),3),(1.0,'Ref_Point_shear2',2)))

		# PBC EDGES AD BC
		for i in range(len(set_edge_AD_nodes.nodes)):
			x = set_edge_AD_nodes.nodes[i].coordinates[0]
			y = set_edge_AD_nodes.nodes[i].coordinates[1]
			z = set_edge_AD_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_AD_nodes.nodes[i].label,)),),
				name='NODE_EDGE_AD_2'+str(i))
			opposite_node=set_edge_BC_nodes.nodes.getByBoundingBox(x-tol,y_max-tol,z-tol,x+tol,y_max+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_BC_2'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_AD_BC'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_AD_2'+str(i),1),(-1.0,'NODE_EDGE_BC_2'+str(i),1),(1.0,'Ref_Point_shear1',1)))
			RVEmodel.Equation(
				name='Equ_edge_AD_BC'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_AD_2'+str(i),2),(-1.0,'NODE_EDGE_BC_2'+str(i),2),(1.0,'Ref_Point_normal',2)))
			RVEmodel.Equation(
				name='Equ_edge_AD_BC'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_AD_2'+str(i),3),(-1.0,'NODE_EDGE_BC_2'+str(i),3),(1.0,'Ref_Point_shear3',1)))

		# EH FG edges
		for i in range(len(set_edge_EH_nodes.nodes)):
			x = set_edge_EH_nodes.nodes[i].coordinates[0]
			y = set_edge_EH_nodes.nodes[i].coordinates[1]
			z = set_edge_EH_nodes.nodes[i].coordinates[2]
			assembly.SetFromNodeLabels(
				nodeLabels=((name_instance,(set_edge_EH_nodes.nodes[i].label,)),),
				name='NODE_EDGE_EH_2'+str(i))
			opposite_node = set_edge_FG_nodes.nodes.getByBoundingBox(x-tol,y_max-tol,z-tol,x+tol,y_max+tol,z+tol)
			j = opposite_node[0].label
			assembly.SetFromNodeLabels(nodeLabels=((name_instance,(j,)),),name='NODE_EDGE_FG_2'+str(i))
			RVEmodel.Equation(
				name='Equ_edge_EH_FG'+str(i)+str(0),
				terms=((1.0,'NODE_EDGE_EH_2'+str(i),1),(-1.0,'NODE_EDGE_FG_2'+str(i),1),(1.0,'Ref_Point_shear1',1)))
			RVEmodel.Equation(
				name='Equ_edge_EH_FG_'+str(i)+str(1),
				terms=((1.0,'NODE_EDGE_EH_2'+str(i),2),(-1.0,'NODE_EDGE_FG_2'+str(i),2),(1.0,'Ref_Point_normal',2)))
			RVEmodel.Equation(
				name='Equ_edge_EH_FG_'+str(i)+str(2),
				terms=((1.0,'NODE_EDGE_EH_2'+str(i),3),(-1.0,'NODE_EDGE_FG_2'+str(i),3),(1.0,'Ref_Point_shear3',1)))

		# PBC VERTICES F E
		RVEmodel.Equation(name='Equ_vertex_F_E'+str(0), terms=((
			1.0,'NODE_VERTEX_F',1),(-1.0,'NODE_VERTEX_E',1),(-1.0,'Ref_Point_shear1',1)))
		RVEmodel.Equation(name='Equ_vertex_F_E'+str(1), terms=((
			1.0,'NODE_VERTEX_F',2),(-1.0,'NODE_VERTEX_E',2),(-1.0,'Ref_Point_normal',2)))
		RVEmodel.Equation(name='Equ_vertex_F_E'+str(2), terms=((
			1.0,'NODE_VERTEX_F',3),(-1.0,'NODE_VERTEX_E',3),(-1.0,'Ref_Point_shear3',1)))

		# PBC VERTICES B A
		RVEmodel.Equation(name='Equ_vertex_B_A'+str(0), terms=((
			1.0,'NODE_VERTEX_B',1),(-1.0,'NODE_VERTEX_A',1),(-1.0,'Ref_Point_shear1',1)))
		RVEmodel.Equation(name='Equ_vertex_B_A'+str(1), terms=((
			1.0,'NODE_VERTEX_B',2),(-1.0,'NODE_VERTEX_A',2),(-1.0,'Ref_Point_normal',2)))
		RVEmodel.Equation(name='Equ_vertex_B_A'+str(2), terms=((
			1.0,'NODE_VERTEX_B',3),(-1.0,'NODE_VERTEX_A',3),(-1.0,'Ref_Point_shear3',1)))

		# PBC VERTICES C D
		RVEmodel.Equation(name='Equ_vertex_C_D'+str(0), terms=((
			1.0,'NODE_VERTEX_C',1),(-1.0,'NODE_VERTEX_D',1),(-1.0,'Ref_Point_shear1',1)))
		RVEmodel.Equation(name='Equ_vertex_C_D'+str(1), terms=((
			1.0,'NODE_VERTEX_C',2),(-1.0,'NODE_VERTEX_D',2),(-1.0,'Ref_Point_normal',2)))
		RVEmodel.Equation(name='Equ_vertex_C_D'+str(2), terms=((
			1.0,'NODE_VERTEX_C',3),(-1.0,'NODE_VERTEX_D',3),(-1.0,'Ref_Point_shear3',1)))

		# PBC VERTICES H G
		RVEmodel.Equation(name='Equ_vertex_H_G'+str(0), terms=((
			1.0,'NODE_VERTEX_H',1),(-1.0,'NODE_VERTEX_G',1),(1.0,'Ref_Point_shear1',1)))
		RVEmodel.Equation(name='Equ_vertex_H_G'+str(1), terms=((
			1.0,'NODE_VERTEX_H',2),(-1.0,'NODE_VERTEX_G',2),(1.0,'Ref_Point_normal',2)))
		RVEmodel.Equation(name='Equ_vertex_H_G'+str(2), terms=((
			1.0,'NODE_VERTEX_H',3),(-1.0,'NODE_VERTEX_G',3),(1.0,'Ref_Point_shear3',1)))

		# PBC VERTICES A D
		RVEmodel.Equation(name='Equ_vertex_A_D'+str(0), terms=((
			1.0,'NODE_VERTEX_A',1),(-1.0,'NODE_VERTEX_D',1),(1.0,'Ref_Point_normal',1)))
		RVEmodel.Equation(name='Equ_vertex_A_D'+str(1), terms=((
			1.0,'NODE_VERTEX_A',2),(-1.0,'NODE_VERTEX_D',2),(1.0,'Ref_Point_shear1',2)))
		RVEmodel.Equation(name='Equ_vertex_A_D'+str(2), terms=((
			1.0,'NODE_VERTEX_A',3),(-1.0,'NODE_VERTEX_D',3),(1.0,'Ref_Point_shear2',2)))

		# PBC VERTEX E H
		RVEmodel.Equation(name='Equ_vertex_E_H'+str(0),terms=((
			1.0,'NODE_VERTEX_E',1),(-1.0,'NODE_VERTEX_H',1),(1.0,'Ref_Point_normal',1)))
		RVEmodel.Equation(name='Equ_vertex_E_H'+str(1),terms=((
			1.0,'NODE_VERTEX_E',2),(-1.0,'NODE_VERTEX_H',2),(1.0,'Ref_Point_shear1',2)))
		RVEmodel.Equation(name='Equ_vertex_E_H'+str(2),terms=((
			1.0,'NODE_VERTEX_E',3),(-1.0,'NODE_VERTEX_H',3),(1.0,'Ref_Point_shear2',2)))

		# no displacement in z direction for z neg surface sets
		fix_z_sets = [
			'NODES_FACE_Z_NEG',
			'NODES_EDGES_AB',
			'NODES_EDGES_BC',
			'NODES_EDGES_CD',
			'NODES_EDGES_AD',
			'NODE_VERTEX_A',
			'NODE_VERTEX_B',
			'NODE_VERTEX_C',
			'NODE_VERTEX_D',
			]
		for fix_z_set in fix_z_sets:
			RVEmodel.DisplacementBC(name='BC_'+fix_z_set, createStepName='Initial',
				region=assembly.sets[fix_z_set], u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET,
				ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='',
				localCsys=None)

		# suppress rigid body modes
		rig_tol = 2*seed_size
		rig_nodes = assembly.instances[name_instance].nodes.getByBoundingBox(
			0.5*x_max - rig_tol,
			0.5*y_max - rig_tol,
			z_min - tol,
			0.5*x_max + rig_tol,
			0.5*y_max + rig_tol,
			z_min + tol
			)

		rig_node = rig_nodes[0].label
		assembly.SetFromNodeLabels(nodeLabels=((name_instance,(rig_node,)),),name='RIGID_NODE')
		RVEmodel.DisplacementBC(name='BC_rigid', createStepName='Initial',
			region=assembly.sets['RIGID_NODE'],
			u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
			amplitude=UNSET, distributionType=UNIFORM, fieldName='',
			localCsys=None)

		log('PBC time: '+str(round(time()-start_pbc, 1))+' seconds.')
	except:
		log('PBC ERROR')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// JOB \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	start_job = time()

	log('Writing .inp file...')
	mdb.jobs[job_name].writeInput(consistencyChecking=OFF)

	try:
		# run job
		if run == True:
			log('Run job '+str(job_name)+'...')
			mdb.jobs[job_name].submit(consistencyChecking=OFF)
			mdb.jobs[job_name].waitForCompletion()
			log('Job completed in '+str(round(time()-start_job, 1))+' seconds.')
		else:
			None
	except:
		log('+ + + RUN JOB ERROR + + +')



# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
# /// POST \\\ #
# ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #

	# write to start.bat
	startbat = open(wdir+str(model_name)+'_start.bat', 'a')
	startbat.write('call abaqus j='+str(job_name)+' interactive cpus='+str(cpus)+' ask_delete=OFF'+'\n')
	startbat.close()

	# save cae file, sleep prevents abaqus from crashing during save
	try:
		log('Saving .cae file...')
		mdb.saveAs(pathName=wdir+str(job_name))
		sleep(10)
	except:
		log('Could not save .cae file.')

	log('RVE '+str(RVE_count)+' duration: '+str(round((time()-start_rve)/60, 1))+' minutes.')
	log('END: RVE '+str(RVE_count)+' / '+str(len(parameters)))
	log('')


log('')
log('Total time: '+str(round((time()-start_0)/60, 1))+' minutes.')
log('')
log('### ----------------------------------------------------- ###')
log('')
log('')
log('')



# # ~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~ #
