import numpy as np
import pandas as pd
from functools import partial
from scipy.integrate import solve_ivp


from pom_model.calculations import (
    calculate_particle_sinking_rate,
    update_particle_radius_values_myRK13,
    calculate_depth_based_temp_FIXME,
    calculate_depth_based_attatchment,
    calculate_f_mult_vol,
    calculate_particle_radius, 
    calculate_cutoff_FIXME, 
    calculate_rs_3_minus_d_FIXME,
    calculate_schmidt_number,
    calculate_unattached_cell_carbon_per_cubic_meter,
    calculate_pseudo_steady_state_mom,
    calculate_cell_density_on_particle_surface

)
from pom_model.constants import (
    DEFAULT_BACTERIAL_GROUP_COUNT,
    DEFAULT_DEPTH_OF_FORMATION,
    DEFAULT_PARTICLE_LABILITY,
    AEROBIC_LMWOM_GROWTH_EFFICIENCY,
    DEFAULT_VISCOSITY,
    DEFAULT_MORTALITY,
    K_LMWOM,
    LOW_BACKGROUND_DETACHMENT_RATE,
    vmom_max,
    PARTICLE_GEOMETRY, 
    FRACTAL_DIMENSION_AGGREGATE, 
    DEFAULT_PARTICLE_DENSITY, 
    SEAWATER_DENISTY, 
    mu,
    GLUCOSE_DIFFUSION,
    CONFIG
)
# REname? update_structure?
def create_structure(NPATCH, dCdt, timestep, bacterial_group_count=DEFAULT_BACTERIAL_GROUP_COUNT):
    """
    
    """
    NPATCH["pom"][timestep] = dCdt[0] # Total pool
    NPATCH["mom"][timestep] = dCdt[1] # Total pool
    NPATCH["depth"][timestep] = dCdt[2] # depth over time, inside myRK

    for i in range(bacterial_group_count):
        NPATCH["bacterial_biomass"][i, timestep] = dCdt[3+i]

    return NPATCH

def create_vector(NPATCH, timestep, bacterial_group_count=DEFAULT_BACTERIAL_GROUP_COUNT):
    dCdt = np.zeros(3 + bacterial_group_count)  # Pre-allocate the array
    dCdt[0] = NPATCH["pom"][timestep]  # total pom pool should be zero when created/set here?
    dCdt[1] = NPATCH["mom"][timestep]  # total mom pool
    dCdt[2] = NPATCH['depth'][timestep]    # depth over time inside myRK

    for i in range(0, bacterial_group_count):
        dCdt[3 + i] = NPATCH['bacterial_biomass'][i, timestep]

    return dCdt

def init_model(varTemp, ymom, mortality, kmom, 
               timesteps=None):
    

    particle_density = CONFIG["particle_density"]
    particle_radius=CONFIG["particle_radius"]
    cells_per_particle=CONFIG["intial_cells_per_particle"]
    initial_biomass=CONFIG["initial_biomass"]
    bacterial_group_count=CONFIG["bacterial_group_count"]
    if not timesteps:
        timesteps=CONFIG["timesteps"]

    new_attached_cells  = calculate_unattached_cell_carbon_per_cubic_meter() # mmol C of cell/m3
    cell_density_initial = calculate_cell_density_on_particle_surface(particle_radius, cells_per_particle)
    schmidt_number = calculate_schmidt_number(DEFAULT_VISCOSITY, GLUCOSE_DIFFUSION)

    PATCH = {
        'pom': np.zeros(timesteps + 1),  # Creates array of zeros
        "mom": np.zeros(timesteps + 1),  #  %mmol.m-3.particle-1,
        "depth": np.zeros(timesteps + 1), 
        "time":np.zeros(timesteps + 1),
        "particle_radius": np.zeros(timesteps + 1),  # FIXME FROM JULIA MODEL
        "bacterial_biomass": np.zeros((len(initial_biomass), timesteps + 1)), # Preallocate a matrix for with number of init_biomass arrays as rows and `timestep + 1 columns
        "cell_density": np.full(bacterial_group_count, cell_density_initial / bacterial_group_count),
        "pool": bacterial_group_count,
        "particle_lability": np.full((1, bacterial_group_count), DEFAULT_PARTICLE_LABILITY), # create 1 by bacterial_group_count matrix,
        "new_cells": new_attached_cells,
        "varTemperature": varTemp, # called a vector in the julia notes but is actually a single value?
        "vmax":vmom_max*ymom-mortality, # FIXME vmom_max = max_growth / ymom ... why multiply by ymom?
        "km": kmom, #half sat rate for lmwom
        "mortality": mortality,
        "schmidt_number": schmidt_number
    }
    PATCH["pom"][0]=particle_density #mmolC.particle-1
    PATCH["mom"][0]=0  #mmolC.m-3
    PATCH["depth"][0] = DEFAULT_DEPTH_OF_FORMATION
    PATCH["particle_radius"][0] = particle_radius
    PATCH["bacterial_biomass"][:,0] = initial_biomass ## FIXME THIS SEEMS OFF
    return PATCH

def particle_based_carbon_dynamics(
        timespan, 
        carbon_conc_array, 
        new_attached_cells,
        schmidt_number,
        kmom=K_LMWOM, 
        tot_bact_off=LOW_BACKGROUND_DETACHMENT_RATE, #Bacterial rate of leaving particle
        mortality=DEFAULT_MORTALITY
    ):
    """
    Carbon dynamics Model per PARTICLE base
    This is run for a particle, for a specific timespan
    ## TODO parallelize for loops if handling larger number of bacterial groups
    """
    # SETUP KEY VARS 

    bacterial_group_count=CONFIG["bacterial_group_count"]
    particulate_organic_matter = carbon_conc_array[0] # mmolC.particle-1
    microbial_organic_matter = carbon_conc_array[1] # mmolC.m-3
    depth = carbon_conc_array[2] # m
    biomass_by_group = np.zeros(bacterial_group_count)
    change_in_biomass_per_group = np.zeros(bacterial_group_count)
    change_in_pom_per_group = np.zeros(bacterial_group_count)
    change_in_c_per_group_FIXME = np.zeros(3+bacterial_group_count)
    vmom_per_group = np.zeros(bacterial_group_count)

    # POM dynamics
    # particle_labilities = np.tile(DEFAULT_PARTICLE_LABILITY, (1, bacterial_group_count)) # FIXME should this be recreated or passed forward?
    particle_labilities = np.full(bacterial_group_count, DEFAULT_PARTICLE_LABILITY)
    rs_3_minus_d = calculate_rs_3_minus_d_FIXME(mu=mu, Af=PARTICLE_GEOMETRY, c_avg_pa_ref=DEFAULT_PARTICLE_DENSITY, fluid_density=SEAWATER_DENISTY, fractal_dimension=FRACTAL_DIMENSION_AGGREGATE)
    f_mult_vol = calculate_f_mult_vol(particulate_organic_matter) # dropped rs_3_minus_d??
    cutoff = calculate_cutoff_FIXME(rs_3_minus_d, FRACTAL_DIMENSION_AGGREGATE)
    particle_radius = calculate_particle_radius(f_mult_vol, cutoff, rs_3_minus_d)
    # use default particle density?
    # is this the only place the sherwood number is used?
    Ws_rk, diffloss, sherwood_nu = update_particle_radius_values_myRK13(particle_radius, cutoff, rs_3_minus_d, schmidt_number)
    TempFunction = calculate_depth_based_temp_FIXME(depth)
    iAttach =calculate_depth_based_attatchment(depth, particle_radius, sherwood_nu, new_attached_cells)# iAttach is cell/day  

    for i in range(0, bacterial_group_count):
        # FIXME potential off by one error -- can 0-3 be stored as metadata?
        biomass_by_group[i] = carbon_conc_array[3+i] # mmolC.particle-1
        change_in_pom_per_group[i] = - particle_labilities[i] * biomass_by_group[i]
    
    total_biomass=np.sum(biomass_by_group)     
    
    # FIXME are these supposed to be arrays??
    tot_vmom=vmom_max*microbial_organic_matter/(microbial_organic_matter+kmom)*TempFunction # Rate of monomer uptake
    tot_u_bo_pa=tot_vmom*AEROBIC_LMWOM_GROWTH_EFFICIENCY; # Bacteria growth rates

    # Calculate pseudo steady state of MOM
    consumption=total_biomass*vmom_max # Consumption
    # where is this coming from? changeed from 1 to i
    production=total_biomass*particle_labilities[i]# beta(1) and beta(2) are the same
    # Solve for MOM when 0=Production-Consumption-Diffusion 
    # a=diffloss # m3/day/particle %Diffusion
    x = calculate_pseudo_steady_state_mom(diffloss, consumption, production)

    # Derivative equations ??
    for i in range(0, bacterial_group_count):
        vmom_per_group[i]=tot_vmom
        change_in_biomass_per_group[i]= biomass_by_group[i]*(tot_u_bo_pa - mortality*TempFunction - tot_bact_off)
        change_in_c_per_group_FIXME[3+i]= change_in_biomass_per_group[i]+iAttach[i]


    dmom_dt=(x-microbial_organic_matter)/CONFIG["dt"] # %Place holder, not used - just back calculated from Steady state approximation
    

    dpom_dt=sum(change_in_pom_per_group[i] for i in range(0, bacterial_group_count))

    carbon_conc_array[0]= dmom_dt # total change to pool of POM
    carbon_conc_array[1]= dpom_dt # total change to pool of MOM
    carbon_conc_array[2]= Ws_rk # depth
    carbon_conc_array = np.array(carbon_conc_array.T) # Why transpose here?
    return carbon_conc_array

def run_model(PATCH, timesteps=None):
    dt=CONFIG["dt"]
    bacterial_group_count=CONFIG["bacterial_group_count"]
    if not timesteps:
        timesteps=CONFIG["timesteps"]
    i = 0
    while i < timesteps: 
        if i % 10 == 0:
            print(f"Step {i} starting")
        
        carbon_conc_array = create_vector(NPATCH=PATCH, bacterial_group_count=bacterial_group_count, timestep=i)
        timespan = (PATCH["time"][i], PATCH["time"][i] + dt) 
        
        # my_func = partial(particle_based_carbon_dynamics, timespan, carbon_conc_array, new_attached_cells=PATCH["new_cells"][i], bacterial_group_count=bacterial_group_count, kmom=K_LMWOM)
        my_func = partial(particle_based_carbon_dynamics, new_attached_cells=PATCH["new_cells"], kmom=K_LMWOM, schmidt_number=PATCH["schmidt_number"])


        solution = solve_ivp(
            my_func, 
            timespan, 
            carbon_conc_array, 
            method='RK45', 
            rtol=1e-6, 
            atol=1e-6
        )

        # Most similar to Tsit5/ode45:
        # method='RK45'
        
        # for perf improvements
        # 1. Radau (stiff systems):
        # method='Radau'

        # # 2. RK23 (if lower accuracy is acceptable):
        # method='RK23'

        # Access solution:
        solution_timepoints = solution.t  # time points not used in julia code?
        y = solution.y  # solution array

        ## ??
        PATCH["time"][i+1]=PATCH["time"][i]+ dt
        
        # Get last row of y (equivalent to y[end,:] in Julia)
        ## FIXME double check the transpose wasnt fucked
        # Store results on PATCH
        results_to_store = y[:,-1] 
        PATCH = create_structure(PATCH, results_to_store, i+1, bacterial_group_count)
        
        # Update sinking rate, compute particle radius for next time point (TODO -- does this cause an issue with last round?)
        curr_particle_radius = PATCH["particle_radius"][i]
        rs_3_minus_d = calculate_rs_3_minus_d_FIXME(mu=mu, Af=PARTICLE_GEOMETRY, c_avg_pa_ref=DEFAULT_PARTICLE_DENSITY, fluid_density=SEAWATER_DENISTY, fractal_dimension=FRACTAL_DIMENSION_AGGREGATE)
        cutoff = calculate_cutoff_FIXME(rs_3_minus_d, FRACTAL_DIMENSION_AGGREGATE)
        sinking_rate = calculate_particle_sinking_rate(curr_particle_radius, cutoff, rs_3_minus_d)
        curr_pom = PATCH["pom"][i] # ; %uM.particle-1 
        f_mult_vol = calculate_f_mult_vol(curr_pom) # f*vol

        PATCH["particle_radius"][i+1] =calculate_particle_radius(f_mult_vol, cutoff, rs_3_minus_d)
        PATCH["depth"][i+1] = PATCH["depth"][i] + sinking_rate * dt # meters sunk in the timestep + curr depth
        i+=1

    return PATCH
