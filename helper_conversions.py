import numpy as np
from helper_constants import METERS_PER_KM, SECONDS_PER_HOUR, HOURS_PER_DAY
km_per_hour_to_m_per_s = lambda x: round(x * METERS_PER_KM / SECONDS_PER_HOUR)
def C_TO_K(T):
    # Convert Celsius to Kelvin
    return T + 273.15 

def DAY_TO_SECONDS(day_count):
    return day_count * HOURS_PER_DAY * SECONDS_PER_HOUR


def OD_TO_CELLS_PER_mL(OD, species_conversion_factor):
    return OD * species_conversion_factor
    
def VOLUME(radius):
    return 4/3*np.pi*radius**3

def SURFACE_AREA(radius):
    return  4 * np.pi * radius**2

def g_per_cm3_TO_kg_per_m3(density):
    # g/cm3 * 1kg/1000 g * (100 cm / 1 m)^3 => 100^3/1000 = 1000
    return density * 1000

def cm_per_s_TO_m_per_day(speed):
    # cm/s * 1m/100 cm * 3600s/1 hr * 24 hr/1 day => * 864
    return speed * 864