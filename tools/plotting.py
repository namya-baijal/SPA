import numpy as np
import matplotlib.pyplot as plt
from .find_tracer import *


def TrajectoryPlotCartesian(data, tracer_idx):
    trajectory_x = data['xmark'][tracer_idx, :]
    trajectory_z = data['zmark'][tracer_idx, :]
    still_in_mesh = trajectory_x < 2000
    fig, ax = plt.subplots()
    ax.scatter(trajectory_x[still_in_mesh], trajectory_z[still_in_mesh], marker='x', label=tracer_idx, c='m', zorder=2)
    Moon = plt.Circle((0, 0), 1750, color='gray', zorder=1, alpha=0.6, label='Lunar surface')
    ax.add_patch(Moon)
    ax.set_aspect('equal')
    plt.legend()
    plt.xlabel('x [km]')
    plt.xlim(np.min(trajectory_x) - 10, np.max(trajectory_x) + 10)
    plt.ylim(np.min(trajectory_z) - 60, np.max(trajectory_z) + 60)
    plt.ylabel('z [km]')
    plt.title('Ballistic Trajectory')
    plt.grid(True)


def ParabolicPlot(sx, sz, tracer_idx):
    trajectory_x = sx
    trajectory_z = sz
    still_in_mesh = trajectory_x < 2000
    fig, ax = plt.subplots(figsize = (10, 10))
    ax.scatter(trajectory_x[still_in_mesh], trajectory_z[still_in_mesh], marker='x', label=tracer_idx, c='m', zorder=2)
    Moon = plt.Circle((0, 0), 1750, color='gray', zorder=1, alpha=0.6, label='Lunar surface')
    ax.add_patch(Moon)
    ax.set_aspect('equal')
    plt.legend()
    plt.xlabel('x [km]')
    plt.xlim(np.min(trajectory_x) - 10, np.max(trajectory_x) + 10)
    plt.ylim(np.min(trajectory_z) - 60, np.max(trajectory_z) + 60)
    plt.ylabel('z [km]')
    plt.title('Ballistic Trajectory')
    plt.grid(True)


def MotionGraphs(average_time, radial_distance_data, radial_velocity_data, radial_acceleration_data, smooth_acceleration_data ,tracer_idx, cross_times):
    fig, ax = plt.subplots(3, sharex=True)
    ax[0].plot(average_time['average_time'], radial_distance_data['r_cumul'][tracer_idx, :-1], label=tracer_idx, c='m')
    for cross_time in cross_times:
        ax[0].scatter(average_time['average_time'][cross_time], radial_distance_data['r_cumul'][tracer_idx, cross_time], marker="x")
    ax[0].set_ylabel('d (km)')
    ax[0].grid(b=True)
    ax[0].legend()
    ax[1].plot(average_time['average_time'], radial_velocity_data['v_r'][tracer_idx], c='m')
    ax[1].set_ylabel('v (km/s)')
    ax[1].grid(b=True)
    ax[2].plot(average_time['average_time'], radial_acceleration_data['acc_r'][tracer_idx], c='pink')
    ax[2].plot(average_time['average_time'], smooth_acceleration_data, c='m')
    # ax[2].scatter([bmp.time_data[bmp.abs_prediction_med_tracer + bmp.first_step]], [0], zorder=3, marker='x', c='k')
    ax[2].grid(b=True)
    ax[2].set_ylabel('a ($km/s^{2}$)')
    ax[2].set_xlabel('t (s)')

    plt.suptitle('Graphs: Radial Displacement, Velocity, Acceleration')
    plt.show()


