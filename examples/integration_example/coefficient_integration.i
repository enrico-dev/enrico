#
# Testing a solution that is second order in space and first order in time
#

[Mesh]
  # This is a cylinder with r=0.5, z=(0,1)
  file = 3D_sideset.exo
  block_id = '1'
  block_name = 'interior'
  boundary_id = '100 200 300'
  boundary_name = 'top bottom wall'
[]

[Variables]
  [./temp]
  [../]
[]

[AuxVariables]
  [./heat_flux_scalar_f_0_l]
    family = SCALAR
    order = TENTH
  [../]
  [./heat_flux_scalar_f_1_l]
    family = SCALAR
    order = TENTH
  [../]
  [./heat_flux_scalar_f_2_l]
    family = SCALAR
    order = TENTH
  [../]
  [./heat_flux_scalar_f_3_l]
    family = SCALAR
    order = TENTH
  [../]
  [./heat_flux_scalar_f_4_l]
    family = SCALAR
    order = TENTH
  [../]
  [./temp_bc_scalar_f_0_l]
    family = SCALAR
    order = TENTH
  [../]
  [./temp_bc_scalar_f_1_l]
    family = SCALAR
    order = TENTH
  [../]
  [./temp_bc_scalar_f_2_l]
    family = SCALAR
    order = TENTH
  [../]
  [./temp_bc_scalar_f_3_l]
    family = SCALAR
    order = TENTH
  [../]
  [./temp_bc_scalar_f_4_l]
    family = SCALAR
    order = TENTH
  [../]
[]

[ICs]
  [./v_ic]
    type = ScalarComponentIC
    variable = 'temp_bc_scalar_f_0_l'
    values = '1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]
  [./temp_ic]
    type = FunctionIC
    variable = temp
    function = '0.0'
  [../]
[]

[Kernels]
  [./HeatSource]
    type = HeatSource
    function = '1.0'
    variable = temp
  [../]
  [./HeatDiff]
    type = HeatConduction
    variable = temp
  [../]
  [./HeatTdot]
    type = HeatConductionTimeDerivative
    variable = temp
  [../]
[]

[Functions]
  # BCFunction just returns 0.0 right now
  [./bc_func]
    type = ConstantFunction

  [../]
  [./legendre_function]
    type = LegendrePolynomial
    l_geom_norm = '0.0 1.0'
  [../]
  [./fourier_function]
    type = FourierPolynomial
  [../]
  [./fl_reconstruction]
    type = FourierLegendreReconstruction
    l_order = 10
    f_order = 5
    l_direction = 2
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    poly_scalars = 'temp_bc_scalar_f_0_l temp_bc_scalar_f_1_l temp_bc_scalar_f_2_l temp_bc_scalar_f_3_l temp_bc_scalar_f_4_l'
  [../]
[]

[BCs]
  [./wall]
    type = FunctionDirichletBC
    variable = temp
    boundary = 'wall'
    function = fl_reconstruction
  [../]
  #[./wall]
  #  type = FunctionDirichletBC
  #  variable = temp
  #  boundary = 'wall'
  #  function = bc_func
  #[../]
[]

[Materials]
  [./k]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '1.0'
    block = 'interior'
  [../]
  [./cp]
    type = GenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '1.0'
    block = 'interior'
  [../]
  [./rho]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.0'
    block = 'interior'
  [../]
[]

[UserObjects]
  # Legendre functions with Fourier order 0
  [./nek_f_0_l_0]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 0
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_1]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 1
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_2]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 2
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_3]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 3
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_4]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 4
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_5]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 5
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_6]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 6
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_7]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 7
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_8]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 8
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_0_l_9]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 9
    f_order = 0
    aux_scalar_name = heat_flux_scalar_f_0_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  # Legendre functions with Fourier order 1
  [./nek_f_1_l_0]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 0
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_1]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 1
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_2]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 2
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_3]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 3
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_4]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 4
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_5]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 5
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_6]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 6
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_7]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 7
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_8]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 8
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_1_l_9]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 9
    f_order = 1
    aux_scalar_name = heat_flux_scalar_f_1_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  # Legendre functions with Fourier order 2
  [./nek_f_2_l_0]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 0
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_1]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 1
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_2]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 2
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_3]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 3
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_4]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 4
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_5]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 5
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_6]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 6
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_7]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 7
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_8]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 8
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_2_l_9]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 9
    f_order = 2
    aux_scalar_name = heat_flux_scalar_f_2_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  # Legendre functions with Fourier order 3
  [./nek_f_3_l_0]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 0
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_1]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 1
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_2]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 2
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_3]
    type = NekSideIntegralVariableUserObject
    variable = temp     
    boundary = wall     
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 3
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_4]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 4
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_5]
    type = NekSideIntegralVariableUserObject
    variable = temp     
    boundary = wall     
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 5
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_6]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 6
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_7]
    type = NekSideIntegralVariableUserObject
    variable = temp     
    boundary = wall     
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 7
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_8]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 8
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_3_l_9]
    type = NekSideIntegralVariableUserObject
    variable = temp     
    boundary = wall     
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 9
    f_order = 3
    aux_scalar_name = heat_flux_scalar_f_3_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_0]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 0
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_1]
    type = NekSideIntegralVariableUserObject
    variable = temp     
    boundary = wall     
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 1
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_2]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 2
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_3]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 3
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_4]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 4
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_5]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 5
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_6]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 6
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_7]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 7
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_8]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 8
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
  [./nek_f_4_l_9]
    type = NekSideIntegralVariableUserObject
    variable = temp
    boundary = wall
    legendre_function_name = 'legendre_function'
    fourier_function_name = 'fourier_function'
    l_direction = 2
    l_order = 9
    f_order = 4
    aux_scalar_name = heat_flux_scalar_f_4_l
    diffusion_coefficient_name = 'thermal_conductivity'
    surface_area_pp = 'surf_area'
  [../]
[]

[Postprocessors]
  [./surf_area]
    type = AreaPostprocessor
    boundary = wall
    execute_on = timestep_begin
  [../]
[]

[Executioner]
  type = Transient
  scheme     = 'Explicit-Euler' # Others available: backward Euler, Crank-Nicholson, etc.
  dt         = 0.001      # Initial timestep size
  start_time = 0        # Starting time
  num_steps  = 50     # Number of Steps
  nl_rel_tol = 1e-6     # Nonlinear relative tolerance
  l_tol      = 1e-6     # Linear tolerance

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  csv = true
  print_perf_log = true
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = GiraffeApp
    sub_cycling = true
    positions = '0 0 0'
    input_files = picard_sub_subcycling.i
  [../]
[]

[Transfers]
  [./to_nek]
    type = MultiAppPolynomialToNek
    direction = to_multiapp
    multi_app = sub
    source_variable = 'heat_flux_scalar_f_0_l heat_flux_scalar_f_1_l heat_flux_scalar_f_2_l heat_flux_scalar_f_3_l heat_flux_scalar_f_4_l '
    to_aux_scalar = 'foo'
  [../]

  [./from_nek]
    type = MultiAppPolynomialToNek
    direction = from_multiapp
    multi_app = sub
    source_variable = 'foo'
    to_aux_scalar = 'temp_bc_scalar_f_0_l temp_bc_scalar_f_1_l temp_bc_scalar_f_2_l temp_bc_scalar_f_3_l temp_bc_scalar_f_4_l'
  [../]
[]
