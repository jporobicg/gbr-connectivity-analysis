## File: connectivity.py
## date : 2022-11-10
## creator: Javier Porobic
## email: javier.porobicgarate@csiro.au
## Description : This code calculates the connectivity matrices by spawning season
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~ ##
## ~     Libraries  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##
import geopandas as gpd
import numpy as np
import math
import matplotlib.pyplot as plt


## ~~~~~~~~~~~~~~~~~~ ##
## ~     Functions  ~ ##
## ~~~~~~~~~~~~~~~~~~ ##
def bathtub_curve(lmbda, v, sigma, t) :
    u = (lmbda * v * pow((lmbda * t), v - 1)) / (1 - sigma * pow((lmbda * t), v))
    return(u)

def piecewise_decay(age, Tcp, lmbda1, lmbda2, v1, v2):
    """
    Piecewise Weibull–Weibull survival model
    Corals (Moneghetti, et al. 2019)
    """
    decay = []
    ## Step 1
    hour = 1/24
    cot = 1
    for age_i in range(0, age):
        age_days = age_i * hour
        if(age_days <= Tcp):
            bathtub_c = bathtub_curve(lmbda1, v1, 0, age_days)
            decay.append(bathtub_c)
            decay_step1 = decay[-1]
            last_Day = age_days
        else :
            #bathtub_c = decay_step1 + bathtub_curve(lmbda2, v2, 0, age_days)
            bathtub_c = bathtub_curve(lmbda2, v2, 0, age_days)
            if(cot == 1) :
                print('lmbda2 ' + str(lmbda2))
                print('v2 ' + str(v2))
                print('age_days ' + str(age_days))
                print('step1 ' + str(math.exp(-decay_step1)))
                print('step2 ' + str(math.exp(-bathtub_c)))
                cot = cot + 2
            decay.append(bathtub_c)
    #print(decay)
    return(math.exp(-decay[-1]))


def piecewise_competency(age, tc, alpha, lmbda, gamma, sigma):
    """
    Piecewise Weibull–exponential competence model

    """
    ## Step 1
    hour = 1/24
    competence = []
    for i_age in range(0, age):
        age_days = i_age * hour

    if age_days < tc:
         competence.append(0)
    else:
        bathtub_c = bathtub_curve(lmbda, gamma, sigma, age_days)
        competence.append(alpha * math.exp(-alpha * (age_days - tc)) * math.exp(- bathtub_c))
    return(competence[-1])
## ~~~~~~~~~~ ##
## ~   Main ~ ##
## ~~~~~~~~~~ ##

parent_dir = "/datasets/work/oa-coconet/work/OceanParcels_outputs"
# Path
path = os.path.join(parent_dir, release_start_day)

## Read files
folderShape="Shape_files/"
originalfile = "gbr1_coral_1m_merged_buffer0p001.shp"


shapefile = folderShape+originalfile
data_shape = gpd.read_file(shapefile_name)

for source_i in range(0, 3805):
    nfile = path + "/GBR1_H2p0_Coral_Release_" + release_start_day + "_Polygon_" +  str(file_id) + '_Wind_3_percent_displacement_field.nc'
    ## lat and lon for all the particles
    oparcels = xr.open_dataset(OP_output)
    lat_part = oparcels['lat']
    lon_part = oparcels['lon']
    particles = gpd.GeoSeries([Point(x, y) for x, y in zip(lon_part, lat_part)])
    for sink_i in range(0, 3805):
        value_index = list(data_shape.loc[data_shape['FID'] == sink_i].index)
        value_index = int("".join(map(str, value_index)))
        polygon = data_shape['geometry'][value_index]
        base_polygon = gpd.GeoSeries([polygon])
        # Are these points inside the polygon?
        p_inside = base_polygon.apply(lambda x: particles.within(x))
        # m = p.to_numpy().reshape(num_particles, 1)
        # valid_indices = np.argwhere(m == True)
        # valid_indices = valid_indices[:, 0]

    piecewise_decay(274, 2.583, 0.4, 0.019, 2.892, 1.716)
    piecewise_competency(274, 3.333, 1.295, 0.002, 0.392, 0.364)



# ## read netcdf file
competency_vector = []
decay_vector = []
settlement_vector = []
days = []


for age_hours in range(1, round(80/ (1/24))):
    days.append(age_hours * (1 / 24))
    decay_vector.append(piecewise_decay(age_hours, 2.583, 0.4, 0.019, 2.892, 1.716))
    competency_vector.append(piecewise_competency(age_hours, 3.333, 1.295, 0.002, 0.392, 0.364))
    settlement_vector.append(decay_vector[-1] * competency_vector[-1])


ax1 = plt.subplot(311)
plt.plot(days, decay_vector)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylabel('rate')
# share x only
ax2 = plt.subplot(312, sharex=ax1)
plt.plot(days, competency_vector)
plt.ylabel('rate')
# make these tick labels invisible
plt.setp(ax2.get_xticklabels(), visible=False)

# share x and y
ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
plt.plot(days, settlement_vector)
plt.xlim(0, 80)
ax1.set_title('Decay rate')
ax2.set_title('Competency rate')
ax3.set_title('Settlement rate [Decay x Competency]')
plt.xlabel('Days')
plt.ylabel('rate')
plt.show()















# -loop for release area
# -- calculate decay
# -- loop trgouth each polygon
# --- if particle > minimum age (CoTS=14 2.5 days for corals) then:
# --- in_polygon :
# ---- umber of particles in the poligon multiplied by their decay rate
