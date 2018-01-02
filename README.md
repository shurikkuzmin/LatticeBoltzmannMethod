# Lattice Boltzmann Method
Here some simple codes for the Lattice Boltzmann Method are presented. Along with the codes I as well present some simple Python scripts to analyze results produced by those codes.

1. Shan-Chen gas-liquid model with the droplet in the center of a domain.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/ShanChen/droplet.png "Droplet simulated by the Shan-Chen method")

2. Binary-liquid free energy model with the droplet on the surface.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/FreeEnergy/droplet_on_surface.png "Droplet simulated by the free energy method")

3. Shan-Chen gas-liquid model programmed on GPU with OpenCL technology.

4. Binary liquid model of the long gas bubbles in 3D microchannels (Bretherton flow). The implementation is parallelized in MPI. Two different implementations are covering full scale simulation and quarter channel simulations with symmetry conditions.

![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/Microchannel3D/benchmark.jpg "Long gas bubbles in the 3D microchannels")

5. The binary liquid droplets on the curved surface. The boundary condition for a phase field is established in a such way as to control the wetting angle.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/CurvedSolid/droplet_on_solid.png "Three droplets on the solid surface")

6. Deposition of the droplet to the curved solid.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/Deposition/deposition.gif "Deposition of the droplet to the surface")

7. Propagation of the Gaussian hill via advection-diffusion equation modeled by two-relaxation-times LB model.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/GaussianHill/gaussian_hill_theoretical.jpg "Theoretical profiles of Gaussian hills")

8. 2D multiple-relaxation-times codes with two different approaches: through Gram-Schmidt procedure and through Hermite polynomials approach.

9. The free-surface code that reads velocities and geometry of the bubble from files, and then imposes the velocity of the walls that bubble stays in the center (reference frame of the bubble). The bubble surface is the subject to the free surface boundary condition, that no normal velocity components do exist. That allows to cleanup the LB multiphase simulation to get rid of non-zero divergence at the bubble surface to be able to perform proper mass transfer calculations. Otherwise non-zero divergence at the surface of the bubble makes mass transfer simulations unstable leading to an infinite accumalation of mass.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/FreeSurface/free_surface.jpg "Velocity streamlines")

10. The diffusion from the boundary of 2D cylinder into the inside cylinder domain. Two boundary conditions of the constant concentration are presented: Inamuro and anti-bounce-back conditions.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/DiffusionCylinder/cylinder_profile.jpg "Concentration dependence on radius with different times")

11. The thermal layer development from the wall with the uniform velocity profile.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/DiffusionChannelUniform/uniform_profile_benchmark.jpg "Boundary conditions for the thermal layer development with the uniform velocity profile")

![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/DiffusionChannelUniform/uniform_profile.jpg "Concentration profiles for the thermal layer development with the uniform velocity profile")

12. The thermal layer development from the wall with the Poiseuille velocity profile.
![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/DiffusionChannelPoiseuille/poiseuille_profile_benchmark.jpg "Boundary conditions for the thermal layer development with the uniform velocity profile")

![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/DiffusionChannelPoiseuille/poiseuille_profile.jpg "Concentration profiles for the thermal layer development with the uniform velocity profile")

13. Drag for steady cylinder simulated by immersed boundary method.

![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/SteadyCylinderImmersedBoundary/steady_cylinder.png "Domain for the steady cylinder")

![](https://github.com/shurikkuzmin/LatticeBoltzmannMethod/blob/master/SteadyCylinderImmersedBoundary/steady_cylinder_drags.png "Drag coefficients from literature")

