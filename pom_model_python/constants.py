import numpy as np
from helper_constants import HOURS_PER_DAY, SECONDS_PER_HOUR
from helper_conversions import g_per_cm3_TO_kg_per_m3

"""
r=0.05;%cm 


%you can either set total cells per particle or cell density at formation
cell_per_particle=1.5E4; %initial cells particle-1 
cell_density_init=cell_per_particle/(4*pi*r*r*1E-4); %cells.m-2.particle-1

F0_max=5e10;%5e10; %cells/m3  concentration of free living bact community 

LMWOM = Low Molecular Weight Organic Matter

"""

"""
CONSTANTS
"""
# ymom -> 
AEROBIC_LMWOM_GROWTH_EFFICIENCY = 0.2 # mmolC_cells/mmolC_LMWOM - aerobic bacteria yield 
# cell_c_content ->  
CELL_CARBON = 20E-15/12*1E3 # mmolC.cell-1 %Lee and Fuhrman 1987, unit fg into mmol C FIXME write out conversion? commented as 'recruited cells ~ late comers' when referred to later in the code

# For OM density and sinking rate calculation
# D_Omand ->
FRACTAL_DIMENSION_AGGREGATE = 1.85  # r_i > r_s  
# D_s ->
FRACTAL_DIMENSTION_SMALL_PARTICLE = 3  # r_i < r_s 
# kmom ->
K_LMWOM = 7.2 # mmolC_MOM.m-3; Half saturation uptake rate for lmwom  [empirically determined]

mu = 1e-3  # kg/m/s
# convert_Alldredge_Omand ->
# non_organic_c_proportion ->
PROPORTION_NON_ORGANIC_C =0.0494 # hard code the B ~ regression between Omand C content and Alldredge C content to account for non-org C fraction in a particle
# Af ->
PARTICLE_GEOMETRY = 1.49 # unitless
REFERENCE_TEMP = 293.15 # K

# fluid_density_ref ->
SEAWATER_DENISTY = 1.028  # g/cm3, Omand 2020

# D ->
GLUCOSE_DIFFUSION_cm_per_second = 6.7E-6  # cm2.s-1 (Diffusion of monomer, D of glucose)
GLUCOSE_DIFFUSION = GLUCOSE_DIFFUSION_cm_per_second * SECONDS_PER_HOUR * HOURS_PER_DAY  # cm2.day-1 (Diffusion of monomer, D of glucose)



"""
DEFAULTS (TODO make configurable in  config)
"""
DEFAULT_BACTERIAL_CELL_RADIUS = 6.2e-5  # cm (for Vol_cell=1um3)
# nbft -> 
DEFAULT_BACTERIAL_GROUP_COUNT = 2 # number of bacterial groups so that initial colonizers and 'recruits' tracked seperately
# ndays ->
DEFAULT_DAYS_PER_TIMESTEP = 5 # day
# dt ->
DEFAULT_dt = 1/HOURS_PER_DAY/6 #  day/step ~= 10 min
# you can either set total cells per particle or cell density at formation
DEFUALT_INITIAL_CELLS_PER_PARTICLE=1.5E4; #initial cells per particle 
DEFAULT_BACTERIAL_MOTILE_DIFFUSIVITY = 5.5E-10 # m2s-1
DEFAULT_DEPTH_OF_FORMATION = 100

# F0_max -> 
DEFAULT_FREE_LIVING_BACTERIA_CONCENTRATION=5e10 #cells/m3 concentration of free living bact community 
# umax -> 
DEFAULT_MAX_GROWTH_RATE = 3.7 # day-1 max low molecular weight organic matter (lmwom) uptake rate [empirically determined]
# mlin -> 
DEFAULT_MORTALITY = 2.5 # day -1 bacterial mortality rate
# motility -> 
DEFAULT_MOTILITY = (1E9+1E10)/2*60 # m2/minute ** FROM CODE
# c_avg_pa_ref ->
DEFAULT_PARTICLE_DENSITY = 1.23  # g/cm3, McDonell 2010, Sicko-Goad 1980, Cram 2018, small particle denisty?

# mybeta ->
DEFAULT_PARTICLE_LABILITY = 150
# r -> 
DEFAULT_PARTICLE_RADIUS = 0.05 # cm  
DEFAULT_TEMP_COEFF = 0.8 # unitless ## Arrhenius equation
DEFAULT_TEMP_ARRHENIUS_MODIFICATION = -4E3 # unitless described as temperature modification regulator in notes ## Arrhenius equation
DEFAULT_VISCOSITY = 1E-2  # cm2.s-1 (Kioebe 2001) KINEMATIC?




"""
FIXME
"""
# L0-> 
LOW_BACKGROUND_DETACHMENT_RATE = 0.25 # day-1
DEFAULT_NUM_DAYS = 100

# vmom_max -> 
vmom_max=DEFAULT_MAX_GROWTH_RATE/AEROBIC_LMWOM_GROWTH_EFFICIENCY #day-1


# initba -> 
DEFAULT_INITIAL_BIOMASS = DEFUALT_INITIAL_CELLS_PER_PARTICLE*CELL_CARBON # mmolC.particle-1 commented as 'early colonizers' when referenced later in the code



initial_biomass = [DEFAULT_INITIAL_BIOMASS, CELL_CARBON] # FIXME -- move to model init # mmolC.particle-1 commented as 'early colonizers' when referenced later in the code
# bo_pas[0] = DEFAULT_INITIAL_BIOMASS
# bo_pas[1] = CELL_CARBON

# FIXME change to a bool
VARY_TEMP = True 


"""
Config
"""

def calculate_schmidt_number(viscosity, diffusion_coefficient):
    """
    diffusion coefficient: units of cm2.day-1
    Schmidt number, unitless, Bianche 2018 Nat Geo (Supp)
    kinematic viscosity of seawater divided by molecular diffusion coefficient for a chemical solute in water
    FIXME: label  
    """
    # converts from per_day diffusion_coefficient
    return viscosity / (diffusion_coefficient/SECONDS_PER_HOUR/HOURS_PER_DAY)  


# SCHMIDT_NUMBER = calculate_schmidt_number(DEFAULT_VISCOSITY, GLUCOSE_DIFFUSION)

# FIXME -- should be better named? 
timesteps = round(DEFAULT_DAYS_PER_TIMESTEP/DEFAULT_dt) # =720 if 5 days and default dt unchanged
PARTICLE_LABILITIES = np.tile(DEFAULT_PARTICLE_LABILITY, (1, DEFAULT_BACTERIAL_GROUP_COUNT)) # C_POC mmol C cell-1 day-1

# you can either set total cells per particle or cell density at formation
CONFIG = {
    "bacterial_group_count": DEFAULT_BACTERIAL_GROUP_COUNT,
    "dt": DEFAULT_dt,
    "initial_biomass": initial_biomass,
    "particle_density": DEFAULT_PARTICLE_DENSITY,
    "particle_radius": DEFAULT_PARTICLE_RADIUS,
    "timesteps": timesteps,
    "intial_cells_per_particle": DEFUALT_INITIAL_CELLS_PER_PARTICLE
}