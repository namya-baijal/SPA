import numpy as np
from scipy.signal import medfilt

def load_ballistic_data(file, end = -1):
    """
    Loads the numpy arrays as dictionaries. Returns time, cartesian, and spherical coordinates for all the tracers at every time step.

    :param file: data file path
       """
    data = np.load(file)
    # store the data as dictionaries
    xmark = data['xmark'][:, :end]  # x coordinate
    ymark = data['ymark'][:, :end]  # y coordinate
    zmark = data['zmark'][:, :end] + 1775  # z coordinate
    trm = np.round(data['TrM'], 0)
    trp = data['TrP']
    time = data['time'][:end]
    r = np.sqrt(xmark ** 2 + ymark ** 2 + zmark ** 2)
    phi = np.arctan2(-xmark, np.maximum(ymark, 1e-12)) + np.pi / 2
    phi_deg = np.rad2deg(phi)
    theta = np.arccos(zmark / r)
    theta_deg = np.rad2deg(theta)

    return {'xmark': xmark, 'ymark': ymark, 'zmark': zmark, 'trm': trm, 'time': time,
            'r': r, 'phi': phi_deg, 'theta': theta_deg , 'trp': trp}


def get_initial_depth(data, tracer_idxs):

    r0 = data['r'][tracer_idxs, 0]
    depth = 1750 - r0

    return {'depth': depth}

def get_average_time(data):
    """
        Finds the average time between any two time steps in the data file

        :param data: Takes the input of the data file

        :return: average_time
        """
    time_avg = (data['time'][1:] + data['time'][:-1]) / 2

    return {'average_time': time_avg}


def get_velocity(data):
    """
    Finds the velocity in 3D cartesian coordinates for the all tracers at every time step.

    :param data: Takes the input of the data file

    :return: v_x ,
             v_y ,
             v_z ,
             v_result
    """
    x_dist = (np.diff(data['xmark'], axis=1))
    y_dist = (np.diff(data['ymark'], axis=1))
    z_dist = (np.diff(data['zmark'], axis=1))
    # calculating each time step
    time_diff = (data['time'][1:] - data['time'][:-1])
    v_x = x_dist / (time_diff + 1e-12)
    v_y = y_dist / (time_diff + 1e-12)
    v_z = z_dist / (time_diff + 1e-12)
    v_resultant = np.sqrt(v_x ** 2 + v_y ** 2 + v_z ** 2)
    v_cylindrical = np.sqrt(v_x**2 +v_y**2)
    return {'v_x': v_x, 'v_y': v_y, 'v_z': v_z, 'v_result': v_resultant, 'v_cylindrical': v_cylindrical}


def get_radial_distance(data):
    """
        Finds the radial distance cumulatively for all the time steps

        :param data: Takes the input of the data file

        :return: radial_distance
        """
    r_cumulative = data['r'] - data['r'][:, 0][:, None]

    return {'r_cumul': r_cumulative}


def get_radial_velocity(data):
    """
        Finds the velocity in 3D speherical coordinates for the all tracers at every time step.

        :param data: Takes the input of the data file

        :return: v_r ,
                 v_phi ,
                 v_theta

        """
    time_diff = (data['time'][1:] - data['time'][:-1])
    r_dist = np.diff(data['r'], axis=1)
    phi_dist = np.diff(data['phi'], axis=1)
    theta_dist = np.diff(data['theta'], axis=1)

    # calculating the spherical coordinate velocities
    v_radial = r_dist / (time_diff+ 1e-12)  # in m/s
    v_phi = phi_dist / (time_diff+ 1e-12)  # in deg/s
    v_theta = theta_dist / (time_diff+ + 1e-12) # in deg/s
    radial_velocity_data = {'v_r': v_radial, 'v_phi': v_phi, 'v_theta': v_theta}
    return radial_velocity_data


def get_radial_acceleration(data, radial_velocity_data):
    """
        Finds the acceleration in 3D spherical coordinates for the all tracers at every time step.

        :param data: Takes the input of the data file

        :return: acc_r: radial acceleration ,
                 acc_phi: phi acceleration ,
                 acc_theta: theta acceleration

        """
    time_diff = (data['time'][1:] - data['time'][:-1])
    acceleration_radial = (radial_velocity_data['v_r'][1:]- radial_velocity_data['v_r'][:-1])/ (time_diff+ 1e-12)
    acceleration_phi = (radial_velocity_data['v_phi'][1:] - radial_velocity_data['v_phi'][:-1]) / (time_diff+ 1e-12)
    acceleration_theta = (radial_velocity_data['v_theta'][1:] - radial_velocity_data['v_theta'][:-1]) / (time_diff+ 1e-12)
    return {'acc_r': acceleration_radial, 'acc_phi': acceleration_phi, 'acc_theta':acceleration_theta}


def get_smooth_acceleration(radial_acceleration_data, tracer_idx, first_step, threshold):

    acc_data_tracer = radial_acceleration_data['acc_r'][tracer_idx].reshape(-1, 1)
    prediction_med_tracer = medfilt(acc_data_tracer.reshape(-1), 51)
    abs_prediction_med_tracer = np.where(np.abs(prediction_med_tracer)[first_step:] < threshold)[0][0]

    return {'abs_pred_med': abs_prediction_med_tracer, 'pred_med': prediction_med_tracer}




