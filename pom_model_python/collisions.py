
import math
def calculate_collision_probability(microbe_density, particle_concentration, 
                                    microbe_radius, particle_radius, 
                                    temperature=293, viscosity=0.001, time=1):
    """
    Calculate the probability of microbes colliding with particles in a fluid.
    
    :param microbe_density: Number of microbes per ml
    :param particle_concentration: Concentration of particles in mg/m^3
    :param microbe_radius: Radius of microbes in meters
    :param particle_radius: Radius of particles in meters
    :param temperature: Temperature in Kelvin (default 293 K, which is ~20°C)
    :param viscosity: Viscosity of the fluid in Pa·s (default 0.001 for water)
    :param time: Time interval in seconds (default 1s)
    :return: Probability of collision for a single microbe in the given time interval
    """
    # Boltzmann constant
    k_B = 1.380649e-23  # J/K
    
    # Convert particle concentration to number density
    # Assuming spherical particles with density close to water (1g/cm^3)
    particle_volume = (4/3) * math.pi * particle_radius**3
    particle_mass = particle_volume * 1000  # kg
    particles_per_m3 = (particle_concentration / 1000) / particle_mass
    particles_per_ml = particles_per_m3 / 1e6
    
    # Calculate diffusion coefficients
    D_m = k_B * temperature / (6 * math.pi * viscosity * microbe_radius)
    D_p = k_B * temperature / (6 * math.pi * viscosity * particle_radius)
    
    # Calculate collision radius
    collision_radius = microbe_radius + particle_radius
    
    # Calculate collision rate constant (Smoluchowski equation)
    collision_rate = 4 * math.pi * (D_m + D_p) * collision_radius
    
    # Calculate probability of collision
    probability = 1 - math.exp(-collision_rate * particles_per_ml * time)
    
    return probability

# Example usage
microbe_density = 1e6  # per ml
particle_concentration = 54  # mg/m^3
microbe_radius = 1e-6  # 1 µm
particle_radius = 1e-7  # 100 nm

probability = calculate_collision_probability(microbe_density, particle_concentration, 
                                              microbe_radius, particle_radius)
print(f"Probability of collision in 1 second: {probability:.2e}")