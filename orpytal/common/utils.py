import numpy as np
import scipy as sp
import math
import os
import pickle


def newton_raphson(fx, dfx, ic, rtol, max_iterations, extra_params=(), attenuation=1.0, atol=1e-6):
    val_list = [ic]

    val = val_list[0]

    for n in range(max_iterations):
        # Find f(x) value
        fx_val = fx(val, *extra_params)

        # Find f'(x) value
        dfdx_val = dfx(val, *extra_params)

        # Do not continue if f'(x)=0
        if not dfdx_val:
            raise ValueError('Cannot continue for f\'(x)=0')

        # Perform Newton-Raphson step
        new_val = val - attenuation * (fx_val / dfdx_val)

        # Add to iterative list (in case we want to see entire list at some point)
        val_list.append(new_val)

        # Stop if we are within the tolerance, else continue iterations
        if np.abs(new_val - val) < rtol and np.abs(fx_val) < atol:
            # print 'Solution found in %i iterations' % len(val_list)
            return val_list[-1]
        else:
            val = new_val

    # Return None if nothing is found
    raise ValueError('No solution found in %i iterations' % max_iterations)


def pickle_dict(filename, pdict):
    file_object = open(filename, 'wb')
    pickle.dump(pdict, file_object)
    file_object.close()


def load_pickle(filename):
    if os.path.getsize(filename) == 0:
        print('{} is an empty file!'.format(filename))
    filename_bytes = filename.encode()
    try:
        file_object = open(filename_bytes, 'rb')
        return pickle.load(file_object)
    except UnicodeDecodeError as ude:
        file_object = open(filename_bytes, 'rb')
        return pickle.load(file_object, encoding='latin1')
