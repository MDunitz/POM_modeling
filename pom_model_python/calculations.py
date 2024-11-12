
from scipy.integrate import solve_ivp

import numpy as np
from pom_model.constants import (
    CELL_CARBON, 
    DEFAULT_FREE_LIVING_BACTERIA_CONCENTRATION, 
    REFERENCE_TEMP, 
    DEFAULT_TEMP_ARRHENIUS_MODIFICATION, 
    DEFAULT_TEMP_COEFF, 
    PARTICLE_GEOMETRY, 
    GLUCOSE_DIFFUSION, 
    DEFAULT_VISCOSITY, 
    DEFAULT_DEPTH_OF_FORMATION, 
    DEFAULT_MOTILITY, 
    SEAWATER_DENISTY,  
    DEFAULT_PARTICLE_DENSITY, 
    mu, 
    FRACTAL_DIMENSION_AGGREGATE, 
    FRACTAL_DIMENSTION_SMALL_PARTICLE, 
    VARY_TEMP, 
    PROPORTION_NON_ORGANIC_C,
    K_LMWOM,
    CONFIG
)
from helper_constants import GRAVITY, HOURS_PER_DAY, SECONDS_PER_HOUR
from helper_conversions import SURFACE_AREA, C_TO_K, cm_per_s_TO_m_per_day

def calculate_schmidt_number(viscosity, diffusion_coefficient):
    """
    Schmidt number, unitless, Bianche 2018 Nat Geo (Supp)
    kinematic viscosity of seawater divided by molecular diffusion coefficient for a chemical solute in water
    FIXME: label  
    """
    # converts from per_day diffusion_coefficient
    return viscosity / (diffusion_coefficient / SECONDS_PER_HOUR / HOURS_PER_DAY)  

def calculate_reynolds_number(sinking_rate_FIXME, viscosity, particle_radius):
    """
    unitless 
    Reynold Number (Moradi 2018, describe how viscous the particle is compared to the environment, 
    depends on the diffusion and advection on the particle 
    Re=sinking*radius/water viscosity (high Re=High sinking)
    ratio of inertical force to viscous force
    """
    return sinking_rate_FIXME * particle_radius / viscosity

def calculate_sherwood_number(schmidt_number, reynolds_number):
    """
    Sherwood number, unitless, Bianche 2018 Nat Geo (Supp)
    relative contribution of diffusive and advective loss from the particle
    =1 => pure diffusive loss
    >> 1 => significant advective loss
    """
    return 1+0.619*(reynolds_number**0.412)*schmidt_number**(1/3)

def calculate_diffloss_FIXME(sherwood_number, particle_radius, diffusion):
    """
    m3/day/particle
    """
    return 4*np.pi*diffusion*sherwood_number*particle_radius*1e-6

# you can either set total cells per particle or cell density at formation
# cell_per_particle=1.5E4; %initial cells particle-1 
# cell_density_init=cell_per_particle/(4*pi*r*r*1E-4); %cells.m-2.particle-1
# RENAME calculate_cell_density -> 
def calculate_cell_density_on_particle_surface(particle_radius, cell_per_particle):
    # Todo update particle size to be in microns?
    surface_area = SURFACE_AREA(particle_radius)
    return cell_per_particle/(surface_area*1E-4)

# Depth and Particle Radius
# POM density and sinking rate calculation
# radius minus diffusion loss?
def calculate_rs_3_minus_d_FIXME(
        mu=mu, 
        Af=PARTICLE_GEOMETRY, 
        c_avg_pa_ref=DEFAULT_PARTICLE_DENSITY, 
        fluid_density=SEAWATER_DENISTY, 
        fractal_dimension=FRACTAL_DIMENSION_AGGREGATE):
    sinking_rate = 0.0012 # m/s 
    ref_particle_volume = 1e-3
    # FIXME pull out stokes law
    if fractal_dimension>1:
        hold = sinking_rate*18*mu
        density_diff = round(c_avg_pa_ref-fluid_density, 4)
        x_var = GRAVITY*Af*density_diff*1000
        fractal_volume = round(ref_particle_volume**(fractal_dimension-1), 8)
        something = hold/(x_var*fractal_volume) # rs is in meter; in which the reference particle is 1e-3 m with 0.0012 m/s sinking rate
    else:
        print(mu, Af, fluid_density, c_avg_pa_ref)

    # rs_3_minus_d=0.0012*18*mu/(g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(1e-3)^(D_Omand-1));%rs is in meter; in which the reference particle is 1e-3 m with 0.0012 m/s sinking rate

    return something

# FIXME rename
def calculate_f_mult_vol(new_c_avg_pa):
    # FIXME Separate function and conversion
    return new_c_avg_pa/PROPORTION_NON_ORGANIC_C/(DEFAULT_PARTICLE_DENSITY/12*1000);


# cm, aka radius of the sub-unit
def calculate_cutoff_FIXME(rs_3_minus_d, fractal_dimension):
    if fractal_dimension < 3:
        return 100*rs_3_minus_d**(1/(3-fractal_dimension))
    else:
        return 'bahhhh'


def calculate_particle_radius(f_mult_vol, cutoff, rs_3_minus_d):
    pseudo_r = (3/4*f_mult_vol/np.pi)**(1/3)# %r calculated if f=1 aka small particles
    if pseudo_r > cutoff:
        pseudo_r = (f_mult_vol/(4/3*np.pi*PARTICLE_GEOMETRY*rs_3_minus_d*100**(3-FRACTAL_DIMENSION_AGGREGATE)))**(1/FRACTAL_DIMENSION_AGGREGATE) # r calcualted if f<1  
    return pseudo_r

def update_particle_radius_values_myRK13(particle_radius, cut_off, rs_3_minus_d, schmidt_number, particle_density=None):
    if not particle_density:
        particle_density = CONFIG["particle_density"]
    # TODO break up into named functions (stokes law and unit change)
    w_aggregate = GRAVITY * PARTICLE_GEOMETRY * (particle_density - SEAWATER_DENISTY) * 1000 * (particle_radius / 100)**(FRACTAL_DIMENSION_AGGREGATE - 1) * rs_3_minus_d / (18 * mu) * 100
    w_smallpa = GRAVITY * PARTICLE_GEOMETRY * (particle_density - SEAWATER_DENISTY) * 1000 * (particle_radius / 100)**(FRACTAL_DIMENSTION_SMALL_PARTICLE - 1) / (18 * mu) * 100


    # can this be represented as a sigmoid function? would that change the performance?
    if particle_radius <= cut_off:
        Ws_cms = w_smallpa  # cm/s
    else:
        Ws_cms = w_aggregate  # cm/s

    Ws_rk = cm_per_s_TO_m_per_day(Ws_cms)
    rey = calculate_reynolds_number(Ws_cms, viscosity=DEFAULT_VISCOSITY, particle_radius=particle_radius) # unitless
    sherwood_nu = calculate_sherwood_number(schmidt_number, rey)
    diffloss = calculate_diffloss_FIXME(sherwood_nu, particle_radius=particle_radius, diffusion=GLUCOSE_DIFFUSION)  # m3/day/particle

    return Ws_rk, diffloss, sherwood_nu



# FIXME renamee?
def calculate_attachment(depth):
    return np.exp(-3e-3 * depth) / np.exp(-3e-3 * DEFAULT_DEPTH_OF_FORMATION)


def stokes_law():
    pass


# FIXME WHATS HAPPENING HERE? 
# needs to be able to process an array? a matrix?
def calculate_temp_based_on_depth(depth_of_pom_formation=DEFAULT_DEPTH_OF_FORMATION):
    return 12*np.exp(-depth_of_pom_formation/150)+12*np.exp(-depth_of_pom_formation/500)+2 


INITIAL_Temp=calculate_temp_based_on_depth() # is this supposed to be an array? if so update TempFun


# FIXME
def TempFun_FIXME(varTemp=VARY_TEMP, TempCoeffArr=DEFAULT_TEMP_COEFF, TempAeArr=DEFAULT_TEMP_ARRHENIUS_MODIFICATION, ARRHENIUS_REF_TEMP=REFERENCE_TEMP, temp=INITIAL_Temp):
    if not varTemp==0:
        return 1
    temp_k = C_TO_K(temp)
    return TempCoeffArr * np.exp(TempAeArr*(1/temp_k-1/ARRHENIUS_REF_TEMP))


def calculate_unattached_cell_carbon_per_cubic_meter(cells_per_ml=DEFAULT_FREE_LIVING_BACTERIA_CONCENTRATION, cell_carbon=CELL_CARBON):
    return cells_per_ml * cell_carbon

# FIXME why is this new attached cells?
# newAttachedCells = calculate_carbon_per_cubic_meter()


def calculate_depth_based_temp_FIXME(depth, varTemp=VARY_TEMP):
    """
    Update values that are depth-dependent for the Carbon dynamics Model per PARTICLE base
    """
    # Calculate TempFunction (temperature punishment)
    Temp = calculate_temp_based_on_depth(depth)
    Temp_100 = calculate_temp_based_on_depth(100)  # Temp at 100 m (original depth "depth")
    TempFun_100 = TempFun_FIXME(varTemp=varTemp, temp=Temp_100)
    # fixme this is different than how it is defined in model init/above? --rename?
    TempFunction = TempFun_FIXME(varTemp=varTemp, temp=Temp) / TempFun_100  # to make sure tempfun is 1 at 100m
    return TempFunction

def calculate_depth_based_attatchment(depth, particle_radius, sherwood_num, new_attached_cels):
    # why is the attachment rate depth dependent??
    # Calculate attachment rate
    free_cell_function = calculate_attachment(depth)
    iAttach = np.zeros(2)
    # TODO break up this function
    iAttach[1] = 4 * np.pi * DEFAULT_MOTILITY * particle_radius * sherwood_num / 100 * free_cell_function * new_attached_cels * 60 * 24  # cells/mins*24*60=cells/day, added to bopa2 pool
    iAttach[0] = 0  # no new cells for original ba ?? 

    return iAttach



def calculate_pseudo_steady_state_mom(diffloss, consumption, production, kmom=K_LMWOM):
    if diffloss==0:
        x=kmom*production/(consumption-production)
    else:
        x=np.sqrt((kmom*production+diffloss*((diffloss*kmom+consumption-production)/(2*diffloss))**2)/diffloss)-(diffloss*kmom+consumption-production)/(2*diffloss)# Pseudo-Steady state: mom(i+1)
    return x


def calculate_density_based_sinking_rate(particle_radius, rs_3_minus_d):
    fixme = particle_radius/100
    density_diff = DEFAULT_PARTICLE_DENSITY-SEAWATER_DENISTY
    ## FIXME BREAK UP
    w_aggregate=GRAVITY*PARTICLE_GEOMETRY*density_diff*1000*fixme**(FRACTAL_DIMENSION_AGGREGATE-1)*rs_3_minus_d/(18*mu)*100 # cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
    return cm_per_s_TO_m_per_day(w_aggregate)
                                                      

def calculate_sinking_rate_for_solid_particle(particle_radius, particle_density=None, fluid_density=SEAWATER_DENISTY):
    if not particle_density:
        particle_density = CONFIG["particle_density"]
    fixme = particle_radius/100
    density_diff = particle_density-fluid_density
    ## FIXME BREAK UP
    w_smallpa=GRAVITY*PARTICLE_GEOMETRY*density_diff*1000*fixme**(FRACTAL_DIMENSTION_SMALL_PARTICLE-1)/(18*mu)*100 #cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
    return cm_per_s_TO_m_per_day(w_smallpa)

def calculate_particle_sinking_rate(cutoff, particle_radius, rs_3_minus_d): 
    """
    particle_radius: cm
    Output as m/day
    small particles are calculated as solid/with default particle density, larger aggregates are calculated based on rs_3_minus_d ??
    """
    if particle_radius <= cutoff:
        return calculate_sinking_rate_for_solid_particle(particle_radius)
    else:
        return calculate_density_based_sinking_rate(particle_radius, rs_3_minus_d)
    

# # f -> calculate_fraction_POM_in_aggregate
# def calculate_fraction_POM_in_aggregate(rs_3_minus_d, cut_off, Af=PARTICLE_GEOMETRY, fractal_dimension=FRACTAL_DIMENSION_AGGREGATE, particle_radius=DEFAULT_PARTICLE_RADIUS):
#     # previously had an element wise operator indicator prior raising to d_omand
#     val = Af*rs_3_minus_d*(particle_radius/100)**(fractal_dimension-3) # unitless, Omand 2020, Fraction of POM in an aggregate
#     if val > 1 or val<=cut_off:
#         val = 1
#     return val

# # NO IDEA WHATS HAPPENING HERE -- Is it actually used?
# def calculate_c_avg_pa_omand_FIXME(f, particle_density, fluid_density=SEAWATER_DENISTY):
#     return f*(particle_density-fluid_density)+fluid_density #  g/cm3, Omand 2020, POM Note4 Trang
