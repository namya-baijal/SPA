import numpy as np
from .find_tracer import *


def get_launch_vector(data, velocity_data, tracer_idx, time_idx):
    longitude = data['phi'][tracer_idx, time_idx]
    theta = data['theta'][tracer_idx, time_idx]
    if theta <= 90:
        latitude = 90 - theta
        direction = 'N'
    else:
        latitude = -1 * (theta - 90)
        direction = 'S'
    v_r = velocity_data['v_result'][tracer_idx, time_idx]
    return v_r, (longitude, 'E'), (latitude, direction)


def get_launch_vector_vec(data, velocity_data, tracer_idxs, time_idxs):
    longitude = data['phi'][tracer_idxs, time_idxs]
    theta = data['theta'][tracer_idxs, time_idxs]
    N = np.where(theta <= 90)[0]
    S = np.where(theta > 90)[0]

    latitude = np.full(theta.shape, 90)
    latitude = latitude - theta

    latitude_direction = np.chararray(theta.shape)
    longitude_direction = np.chararray(theta.shape)

    latitude_direction[N] = 'N'
    latitude_direction[S] = 'S'

    longitude_direction[:] = 'E'
    v_r = velocity_data['v_result'][tracer_idxs, time_idxs]
    return v_r, (longitude, longitude_direction), (latitude, latitude_direction)


def get_v_tilda(velocity_data, tracer_idx, R, g, time_idx):

    """
    Function to calculate the launch velocity relative to the "first cosmic velocity" or the resultant velocity

    """
    v_tilda = velocity_data['v_result'][tracer_idx,time_idx ] / (np.sqrt(R * g))  # in km/s

    return v_tilda


def get_v_tilda_vec(velocity_data, tracer_idxs, time_idxs, R, g):

    """
    Function to calculate the launch velocity relative to the "first cosmic velocity" or the resultant velocity

    """
    v_tilda = velocity_data['v_result'][tracer_idxs, time_idxs] / (np.sqrt(R * g))  # in km/s

    return v_tilda


def get_ballistic_vectors(data, velocity_data, tracer_idx, time_idx, h_launch):
    R0 = 1750  # radius of the moon in km

    g0 = 0.00162  # gravity of the moon in km/s^2
    g = g0* R0**2 / (R0 + h_launch) **2
    #return launch_angle, height and distance of a given tracer
    v_tilda = get_v_tilda(velocity_data, tracer_idx, R0+h_launch, g, time_idx)
    beta = np.arccos(data['zmark'] / data['r'])[:, :-1]  # angle at the centre of the moon
    beta_deg = np.rad2deg(beta)
    # print(beta_deg.shape)
    # print('beta =',beta_deg)
    alpha = np.arcsin(velocity_data['v_z'] / (velocity_data['v_result'] + 1e-12))
    alpha_deg = np.rad2deg(alpha)


    theta_launch = alpha_deg[tracer_idx, time_idx] + beta_deg[tracer_idx, time_idx]
    theta_launch_rad = np.deg2rad(theta_launch)

    #calc bearing
    if velocity_data['v_x'][tracer_idx, time_idx] > 0:
        azimuth = np.arcsin(velocity_data['v_y']/ (velocity_data['v_cylindrical'] + 1e-12))

    else:
        azimuth = np.pi - np.arcsin(velocity_data['v_y'] / (velocity_data['v_cylindrical'] + 1e-12))

    azimuth_deg = np.rad2deg(azimuth[tracer_idx, time_idx])

    # print("theta launch", theta_launch_calc)
    # range between the launch and impact
    d = abs((((velocity_data['v_result'][tracer_idx, time_idx] ** 2 * np.sin(2 * theta_launch_rad)) / g) / (
        np.sqrt(1 - (2 - v_tilda ** 2) * v_tilda ** 2 * (np.cos(theta_launch_rad)) ** 2))) + (h_launch / np.tan(theta_launch_rad)))

    # print('d=', d)

    # maximum height of the projectile above the planetary surface
    height = (velocity_data['v_result'][tracer_idx, time_idx] ** 2 * (np.sin(theta_launch_rad)) ** 2 / g) / (
                1 - v_tilda ** 2 + np.sqrt(1 - (2 - v_tilda ** 2) * v_tilda ** 2 * (np.cos(theta_launch_rad)) ** 2))

    return theta_launch, height, d, azimuth_deg, alpha, beta


def get_ballistic_vectors_vec(data, velocity_data, tracer_idxs, time_idxs, h_launch, R0, g0):
    g = g0* R0**2 / (R0 + h_launch) **2
    #return launch_angle, height and distance of a given tracer
    v_tilda = get_v_tilda_vec(velocity_data, tracer_idxs, time_idxs, R0+h_launch, g) # select tracers and select times
    beta = np.arccos(data['zmark'] / (data['r']+ 1e-12))[:, :-1]  # angle at the centre of the moon
    beta_deg = np.rad2deg(beta) # beta degrees for all tracers at all times
    # print(beta_deg.shape)
    # print('beta =',beta_deg)
    alpha = np.arcsin(velocity_data['v_z'] / (velocity_data['v_result'] + 1e-12))
    alpha_deg = np.rad2deg(alpha) # alpha for all tracers at all times

    theta_launch = alpha_deg[tracer_idxs, time_idxs] + beta_deg[tracer_idxs, time_idxs]
    theta_launch_rad = np.deg2rad(theta_launch) # for select tracers at select timees

    azimuth = np.arcsin(velocity_data['v_y'][tracer_idxs, time_idxs] / (velocity_data['v_cylindrical'][tracer_idxs, time_idxs] + 1e-12)) # all tracers all time_steps
    condition = np.where(velocity_data['v_x'][tracer_idxs, time_idxs] <= 0)[0]
    azimuth[condition] = -1 * (azimuth[condition] - np.pi)
    azimuth_deg = np.rad2deg(azimuth)

    # range between the launch and impact
    # namya made a change
    d = abs((((np.square(velocity_data['v_result'][tracer_idxs, time_idxs]) * np.sin(2 * theta_launch_rad)) / g) / (
        np.sqrt(1 - (2 - np.square(v_tilda)) * np.square(v_tilda) * (np.square(np.cos(theta_launch_rad))))))
        + (np.full(theta_launch_rad.shape, h_launch) / np.tan(theta_launch_rad)))

    # maximum height of the projectile above the planetary surface
    height = (np.square(velocity_data['v_result'][tracer_idxs, time_idxs]) * (np.square(np.sin(theta_launch_rad))) / g) / (
                1 - np.square(v_tilda) + np.sqrt(1 - (2 - np.square(v_tilda)) * np.square(v_tilda) * (np.square(np.cos(theta_launch_rad)))))

    return theta_launch, height, d, azimuth_deg

def track(lat, lon, bearing, distance):
    """Update longitude and latitude data given an initial
    bearing and distance travelled.

    Parameters
    ----------
    lat: arraylike
        Initial latitude (as length n 1d arraylike).
    lon: arraylike
        Initial longitude (as length n 1d arraylike).
    bearing: arraylike
        Initial bearing (forward azimuth) in degrees (as length n 1d arraylike).
    distance: arraylike
        Distance travelled  (as length n 1d arraylike or length n x s 2d arraylike)

    Returns
    -------
    ndarray
        Final latitudes (as length n 1d arraylike or length n x s arraylike).
    ndarray
        Final longitude (as length n 1d arraylike or length n x s arraylike).


    Examples
    --------


         print(track([54.0], [0.0], [45.0], [1.0e5.0]))
    (array(5.509e+01), array(3.981e+00))
    """

    d = np.asarray(distance)/1750
    b = np.pi*np.asarray(bearing).reshape((-1, 1))/180.0

    rlat = np.pi*np.asarray(lat).reshape((-1, 1))/180.0
    rlon = np.pi*np.asarray(lon).reshape((-1, 1))/180.0

    final_lat = np.arcsin(np.sin(rlat)*np.cos(d)
                          + (np.cos(b)*np.cos(rlat))*np.sin(d))
    final_lon = (rlon + np.arctan2(np.sin(b)*np.cos(rlat)*np.sin(d),
                                   np.cos(d)-np.sin(rlat)*np.sin(final_lat)))

    return (np.squeeze(final_lat*180/np.pi),
            np.squeeze(final_lon*180/np.pi))



def track_vec(lat, lon, bearing, distance):
    """Update longitude and latitude data given an initial
    bearing and distance travelled.

    Parameters
    ----------
    lat: arraylike
        Initial latitude (as length n 1d arraylike).
    lon: arraylike
        Initial longitude (as length n 1d arraylike).
    bearing: arraylike
        Initial bearing (forward azimuth) in degrees (as length n 1d arraylike).
    distance: arraylike
        Distance travelled  (as length n 1d arraylike or length n x s 2d arraylike)

    Returns
    -------
    ndarray
        Final latitude of landing position (as length n 1d arraylike or length n x s arraylike).
    ndarray
        Final longitude of landing position (as length n 1d arraylike or length n x s arraylike).


    Examples
    --------

         print(track([54.0], [0.0], [45.0], [1.0e5.0]))
    (array(5.509e+01), array(3.981e+00))
    """

    d = np.asarray(distance).reshape(-1, 1) / 1750 #change this to the radius of the moon
    b = np.pi * np.asarray(bearing).reshape((-1, 1)) / 180.0

    rlat = np.pi * np.asarray(lat).reshape((-1, 1)) / 180.0
    rlon = np.pi * np.asarray(lon).reshape((-1, 1)) / 180.0

    final_lat = np.arcsin(np.sin(rlat) * np.cos(d)
                          + (np.cos(b) * np.cos(rlat)) * np.sin(d))
    final_lon = (rlon + np.arctan2(np.sin(b) * np.cos(rlat) * np.sin(d),
                                   np.cos(d) - np.sin(rlat) * np.sin(final_lat)))

    return (np.squeeze(final_lat * 180 / np.pi),
            np.squeeze(final_lon * 180 / np.pi))

def parabolic_path(tracer_idx, time_idx, theta_launch, beta, data, velocity_data, g):
    t = np.linspace(0, 1000, 501)
     # convert from time idx to timestep
    V_r = velocity_data['v_result'][tracer_idx, time_idx]
    sx_i = - 1 * V_r * np.cos(np.deg2rad(theta_launch)) * t
    sz_i = V_r * np.sin(np.deg2rad(theta_launch)) * t - 0.5 * g * np.square(t)


    b = beta[tracer_idx, time_idx]

    sx = np.cos(b) * sx_i - np.sin(b) * sz_i
    sz = np.sin(b) * sx_i + np.cos(b) * sz_i
    return sx, sz


def haversine_moon(lat1, long1, lat2, long2):

    """
    Calculates the great circle distance between two points on the moon. Also calculates the bearing from point 1 to point 2.
    Takes the input starting and ending lat, lons in degrees.

    """
    lat1 = np.deg2rad(lat1)
    long1 = np.deg2rad(long1)
    lat2 = np.deg2rad(lat2)
    long2 = np.deg2rad(long2)
    radius = 1750  #radius of the moon
    distance = 2 * radius * np.arcsin(np.sqrt(np.square(np.sin(lat2-lat1))+ ((1-np.square(np.sin(lat2-lat1))) - np.square(np.sin(lat2+lat1)))* np.square((long2-long1)/2)))
    num = (np.cos(lat2)* np.sin(long2 -long1))
    denom = (np.cos(lat1)* np.sin(lat2)) - (np.sin(lat1)*np.cos(lat2)* np.cos(long2-long1))
    bearing = np.arctan(num/denom)

    return distance, bearing




