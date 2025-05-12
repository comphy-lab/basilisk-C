CFLAGS +=

include $(BASILISK)/Makefile.defs

anisotropy.tst  					: CFLAGS += -fopenmp

atomisation_skeleton.tst  			: CFLAGS += -O2 -fopenmp

atomisation3D_skeleton.tst   		: CFLAGS += -O2

advection_LS.tst 					: CFLAGS += -catch

bump-vdw.tst : CLAGS += -fopenmp -O2

chen1D.tst 							: CFLAGS +=

cube.tst   							: CFLAGS += -fopenmp

crystal_growth3D.tst : CFLAGS += -fopenmp

crystal_growth_grid6.c 			: cube_template.c
					ln -sf cube_template.c crystal_growth_grid6.c

crystal_growth_grid6.tst 			: CFLAGS += -O2 -fopenmp -DGRIDLEVEL=6

crystal_growth_grid7.c 			: cube_template.c
					ln -sf cube_template.c crystal_growth_grid7.c

crystal_growth_grid7.tst 			: CFLAGS += -O2 -fopenmp -DGRIDLEVEL=7

crystal_growth_noGT.c 			: cube_template.c
					ln -sf cube_template.c crystal_growth_noGT.c

crystal_growth_noGT.tst 			: CFLAGS += -O2 -fopenmp -DGRIDLEVEL=7 -DEPSK=0. -DEPSV=0.

crystal_growth_grid8.c 			: cube_template.c
					ln -sf cube_template.c crystal_growth_grid8.c

crystal_growth_grid8.tst 			: CFLAGS += -O2 -DGRIDLEVEL=8

cube_aniso.tst						: CFLAGS += -fopenmp

cube_secondary.tst					: CFLAGS += -fopenmp -O2

cube_linear_solvability.tst			: CFLAGS += -fopenmp -O2

cube_sixfold.tst					: CFLAGS += -fopenmp -O2

cube_sixfold_10.tst					: CFLAGS += -fopenmp -O2

cube_template.tst 			: CFLAGS += -fopenmp -O2

cube3D.tst   						: CFLAGS += -fopenmp -Werror -O2

cube3D_dendrite.tst					: CFLAGS += -fopenmp -Werror -O2

crystal.tst  						: CFLAGS += -fopenmp -O2

crystal_smooth.tst  				: CFLAGS += -fopenmp

curvature_LS.tst                    : CFLAGS +=

dirichlet.tst 						: CFLAGS += -DDIRICHLET=1 -DMTRACE=3
	
dirichlet.c                         : neumann_double.c
					ln -s neumann_double.c dirichlet.c

distanceToEllipsoid.tst 			: CFLAGS += -O2

distance.tst 						: CFLAGS += -O2
distance_LS.tst 					: CFLAGS += -O2

double_diff.tst 					: CFLAGS += -fopenmp

dirichlet_expand_speed_calculated.tst : CFLAGS += -fopenmp

equilibrium_curvature.tst           : CFLAGS +=

extrapolation_ghigo.tst : CFLAGS +=

Favier_Ra-Be.tst 				: CFLAGS += -O2 -fopenmp

frank_spheres.tst                   : CFLAGS += -fopenmp

frank_3D.tst 		                : CFLAGS += -fopenmp

frank_nospeed.tst                   : CFLAGS += -fopenmp

gaussian_solidification.tst  		: CFLAGS += -fopenmp

grid_scalar.tst  : CFLAGS += -catch

granularCriticality.tst : CFLAGS += -O2

initial_speed.tst 					: CFLAGS += 

meanflow.tst 						: CFLAGS +=

meanflow3D.tst 						: CFLAGS += -fopenmp -DTRACE=2

mydirichlet_expand_1D.tst  			: CFLAGS += -fopenmp

mullins_amplitude.tst  				: CFLAGS += -fopenmp

mullins_sekerka.tst  				: CFLAGS += -fopenmp

mullins_timestep.tst  				: CFLAGS += -fopenmp

moving_embed.tst  					: CFLAGS += -fopenmp

Narrow_band.tst : CC = mpicc -D_MPI=2

phasesep.tst : CFLAGS += -O2

porous.tst        					: CFLAGS +=

Pinned_test.tst  					: CFLAGS += -fopenmp

recons_speed.tst					: CFLAGS += -fopenmp

reinit_LS.tst						: CFLAGS += -catch
	
reinit_circle.tst					: CFLAGS += -fopenmp

reinit_redistance.tst				: CFLAGS += -catch

rivulet.tst      					: CFLAGS += -fopenmp

rivulet2.tst      					: CFLAGS += -fopenmp

rivulet3D.tst      					: CFLAGS += -fopenmp

reversed.tst      					: CFLAGS += -fopenmp

reversed_skeleton.tst      			: CFLAGS += -fopenmp -O2

sessile-PF.tst : CFLAGS += -O2 -fopenmp

skeleton_snake.tst 					: CFLAGS += -O2

simili-RP-instab.tst			: CFLAGS += -O2

solidification_Gibbs_Thomson.tst 	: CFLAGS += -events

solidification_mwe.tst 				: CFLAGS += -fopenmp

solidification_sine.tst				: CFLAGS += -fopenmp

spinodal.tst 								: CFLAGS += -O2

substrate_dendrite.tst : CFLAGS += -O2 -fopenmp

test_advection.tst					: CFLAGS += -fopenmp

test_curve.tst						: CFLAGS += -fopenmp

test_recons.tst						: CFLAGS +=

test_recons3D.tst					: CFLAGS +=

test_recons_circle.tst				: CFLAGS +=

test_bicubic.tst					: CFLAGS +=

testRK3.tst 						: CFLAGS += 

update.tst     					   	: CFLAGS += -fopenmp

vof_rotation.tst 					: CFLAGS += -O2

zalesak.tst    					   	: CFLAGS += -fopenmp -O2

zalesak_centered.tst			   	: CFLAGS += -O2 -fopenmp

extrapolation_gaussian.tst : CFLAGS += -O2 -grid=multigrid

Favier_Ra-Be.tst: CFLAGS += -O2 -fopenmp