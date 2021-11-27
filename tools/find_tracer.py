import numpy as np
import argparse
from .data_utils import load_ballistic_data


def findTracer3D(data, time_step, xtarg, ytarg, ztarg):
    """
    Find the index of the closest tracer(s) to the given x, y and z
    coordinates. If more than one tracer matches that condition, then
    return list of all matching tracers.

    :param data: Path to the data file
    :param time_step: Time step at which tracer is closest to the input x,y,z points
    :param xtarg: X coordinate
    :param ytarg: Y coordinate
    :param ztarg: Z coorindate

    :return: tracer_idx: Index of the tracer(s) closest to the given coordinates
    """
    # distances between all tracers and the target coordinate at a given time
    rmxt = np.sqrt((data['xmark'][:, time_step] - xtarg) ** 2. + (data['ymark'][:, time_step] - ytarg) ** 2. + (
            data['zmark'][:, time_step] - ztarg) ** 2.)

    # # # # list all tracers that are closest to the point
    tracer_idx = list(np.where(rmxt == rmxt.min())[0])

    # # return the tracer index. if more than one, return a list of indices.
    return tracer_idx[0] if len(tracer_idx) == 1 else tracer_idx


def BallisticCheck(data, tracer_idx, h_launch):
    """
    Function to determine if the ejecta follows a ballistic trajectory. If the particle gets to a radial distance  > launch height above the lunar surface,
    it is classified as ballistic. Returns True/False and the first time step index at which the condition is satisfied.

    :param data:
    :param tracer_idx:
    :param h_launch:

    """

    r_val = data['r'][tracer_idx, :]
    # return True/False time_index
    if not np.any(r_val > 1750 + h_launch):
        return False, None
    else:
        time_idx = np.where(r_val > 1750 + h_launch)
        time_idx = time_idx[0][0]
        return True, time_idx

def BallisticCheckVec(data, h_launch):
    r_vals = data['r'][:, :]

    r_vals[:, :8] = 0

    tracer_idxs, _ = np.where((r_vals > 1750 + h_launch)* (data['xmark'] <2000))
    tracer_idxs = list(set(tracer_idxs))
    time_steps = np.argmax(r_vals[tracer_idxs, :] > 1750 + h_launch, axis=1)
    return np.array(tracer_idxs), np.array(time_steps)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Finding the Tracer Index',
                                     description='Find the index of the closest tracer(s) for a given x,y,z input')
    parser.add_argument('--data', type=str, help='path to data file')
    parser.add_argument('--xtarg', type=float, help='x coordinate of choice')
    parser.add_argument('--ytarg', type=float,
                        help='y coordinate of choice. This input should be 0 for tracers on the symmetry plane')
    parser.add_argument('--ztarg', type=float, help='z coordinate of choice')
    parser.add_argument('--time', type=int, help='timestep')

    args = parser.parse_args()
    xtarg = args.xtarg
    ytarg = args.ytarg
    ztarg = args.ztarg
    time = args.time
    datafile = args.data

    data = load_ballistic_data(datafile)
    data['zmark'] += 1775

    tracer_idx = findTracer3D(data, time, xtarg, ytarg, ztarg)
    # print(ballisticCheck(tracer_idx, 100))
    print('tracer index=', tracer_idx)

