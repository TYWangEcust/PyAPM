import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from math import factorial

def model_fuct_lj126(x, a, b):
	return 4*a * ((b/x)**12-(b/x)**6)

def model_fuct_harmonic(x, a, b):
	return a*(x-b)**2

def gaussian(r, a, w, r0):
    constant = np.sqrt(np.pi / 2)
    return a / (w * constant) * np.exp(-2 * ((r - r0) / w) ** 2)

def multi_gaussian_model(r, *params):
    return sum(gaussian(r, a, w, r0) for a, w, r0 in zip(params[::3], params[1::3], params[2::3]))

def fitting_lj(x, y):
	is_inf = np.isinf(y)
	y = y[~is_inf]
	x = x[~is_inf]
	is_nan = np.isnan(y)
	y = y[~is_nan]
	x = x[~is_nan]
	initial_guess = [1, np.min(x[np.where(y < 0)])]
	params, params_covariance = curve_fit(model_fuct_lj126, x, y, p0=initial_guess, maxfev = 100000)
	return params

def fitting_harmonic_bond(x, y):
	is_inf = np.isinf(y)
	y = y[~is_inf]
	x = x[~is_inf]
	is_nan = np.isnan(y)
	y = y[~is_nan]
	x = x[~is_nan]
	b0 = np.mean(x[np.where(y==0)])
	if np.isnan(b0) == True:
		b0 = np.mean(x[np.isnan(np.gradient(y))])
	initial_guess = [1, b0]
	params, params_covariance = curve_fit(model_fuct_harmonic, x, y, p0=initial_guess, bounds = ([0,b0-0.1],[1000,b0+0.1]), maxfev = 100000)
	return params

def fitting_harmonic_angle(x, y):
    is_inf = np.isinf(y)
    y = y[~is_inf]
    x = x[~is_inf]
    is_nan = np.isnan(y)
    y = y[~is_nan]
    x = x[~is_nan]
    b0 = np.mean(x[np.where(y==0)])
    initial_guess = [1, b0]
    params, params_covariance = curve_fit(model_fuct_harmonic, x, y, p0=initial_guess, bounds = ([0,b0-30],[5,b0+30]), maxfev = 100000)
    params[0] = params[0]*(180/np.pi)**2
    return params

def fitting_gaussian_bond(x, y):
    is_inf = np.isinf(y)
    y = y[~is_inf]
    x = x[~is_inf]
    is_nan = np.isnan(y)
    y = y[~is_nan]
    x = x[~is_nan]
    peaks, _ = find_peaks(y, height=0.25, distance=5)
    n = len(peaks)
    initial_guess = []
    bound_l = []
    bound_u = [] 
    for i, peak_index in enumerate(peaks):
        peak_x = x[peak_index]
        initial_guess.extend([0.1, 0.1, peak_x])
        bound_l.extend([0.0001, 0.0001, peak_x - 0.5])
        bound_u.extend([1, 1, peak_x + 0.5])
    params, params_covariance = curve_fit(multi_gaussian_model, x, y, p0=initial_guess, bounds=(bound_l, bound_u), maxfev=100000)
    return n, params 

def fitting_gaussian_angle(x, y):
    is_inf = np.isinf(y)
    y = y[~is_inf]
    x = x[~is_inf]
    is_nan = np.isnan(y)
    y = y[~is_nan]
    x = x[~is_nan]
    peaks, _ = find_peaks(y, height=0.25, distance=10)
    n = len(peaks)
    initial_guess = []
    bound_l = []
    bound_u = []
    x = np.deg2rad(x)
    for i, peak_index in enumerate(peaks):
        peak_x = x[peak_index]
        initial_guess.extend([0.1, 0.1, peak_x])
        bound_l.extend([0.0001, 0.0001, peak_x - 0.1])
        bound_u.extend([1, 1, peak_x + 0.1])
    params, params_covariance = curve_fit(multi_gaussian_model, x, y, p0=initial_guess, bounds=(bound_l, bound_u), maxfev=100000)
    for i in range(n):
        params[i*3+2] = np.rad2deg(params[i*3+2])
    return n, params

def Boltzmann_inversion(x , temp):
	y =  -8.31/4186*temp*np.log(x)
	return y

def iteration_Boltzmann(pot_func, dist_func, target_func, temp, alpha = 0.05):
	predict_pot = pot_func + alpha * 8.31/4186*temp * np.log(dist_func/target_func)
	return predict_pot

def match(rdf, target_rdf):
	fit_tunc = 1- np.sum(np.abs(rdf-target_rdf))/np.sum(np.abs(rdf)+np.abs(target_rdf))
	return fit_tunc

def diff(y, x):#Central difference method
	is_inf = np.isinf(y) #Remove infinities
	y[is_inf] = y[np.argmin(x[~is_inf])]
	is_nan = np.isnan(y) #Remove null items
	y[is_nan] = 0
	n = len(x)
	y[-1] = 0
	dy_dx = np.zeros(n)
	for i in range(1, n-1):
		dx = x[i+1] - x[i-1]
		dy = y[i+1] - y[i-1]
		dy_dx[i] = dy / dx
	dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])#Handle boundary points
	dy_dx[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2])
	return x, y, -dy_dx

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
    Savitzky-Golay filtering is a signal smoothing method that removes noise or reduces noise interference in a signal.
    It uses a set of sliding windows and estimates smoothed values by performing polynomial fitting on the data within each window.
    The smoothing effect of window_length on the curve: the smaller the window_length value, the closer the curve is to the true curve; the larger the window_length value, 
    the stronger the smoothing effect (note: this value must be a positive odd integer).
    The smoothing effect of the k value on the curve: the larger the k value, the closer the curve is to the true curve; the smaller the k value, the stronger the smoothing.
    In addition, when the k value is large, limited by the window length, the fitting may become problematic, and high-frequency curves may become straight lines.
    """
    if not (isinstance(window_size, int) and isinstance(order, int)):#判断是否为int类型
        raise ValueError("window_size and order must be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError('window_size must be a positive odd number')
    if window_size < order + 2:
        raise TypeError('window_size is too small for the polynomials order')

    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window +1)])#np.mat将一个数组转换为矩阵。
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

def find_nearest(array, target):
    """Find array component whose numeric value is closest to 'target'. """
    idx = np.abs(array - target).argmin()
    return idx, array[idx]

def tail_correction(r, V, r_switch):
    """Apply a tail correction to a potential making it go to zero smoothly.

    Parameters
    ----------
    r : np.ndarray, shape=(n_points,), dtype=float
        The radius values at which the potential is given.
    V : np.ndarray, shape=r.shape, dtype=float
        The potential values at each radius value.
    r_switch : float, optional, default=pot_r[-1] - 5 * dr
        The radius after which a tail correction is applied.

    References
    ----------
    .. [1] https://codeblue.umich.edu/hoomd-blue/doc/classhoomd__script_1_1pair_1_1pair.html

    """
    r_cut = r[-1]
    idx_r_switch, r_switch = find_nearest(r, r_switch)

    S_r = np.ones_like(r)
    r = r[idx_r_switch:]
    S_r[idx_r_switch:] = ((r_cut ** 2 - r ** 2) ** 2 *
                          (r_cut ** 2 + 2 * r ** 2 - 3 * r_switch ** 2) /
                          (r_cut ** 2 - r_switch ** 2) ** 3)
    is_nan = np.isnan(V)
    V[is_nan] = 0.0
    return V * S_r


def head_correction(r, V, previous_V, form='linear'):
    """Apply head correction to V making it go to a finite value at V(0).

    Parameters
    ----------
    r : np.ndarray, shape=(n_points,), dtype=float
        The radius values at which the potential is given.
    V : np.ndarray, shape=r.shape, dtype=float
        The potential values at each radius value.
    previous_V : np.ndarray, shape=r.shape, dtype=float
        The potential from the previous iteration.
    form : str, optional, default='linear'
        The form of the smoothing function used.

    """
    if form == 'linear':
        correction_function = linear_head_correction
    elif form == 'exponential':
        correction_function = exponential_head_correction
    else:
        raise ValueError('Unsupported head correction form: "{0}"'.format(form))

    for i, pot_value in enumerate(V[::-1]):
        # Apply correction function because either of the following is true:
        #   * both current and target RDFs are 0 --> nan values in potential.
        #   * current rdf > 0, target rdf = 0 --> +inf values in potential.
        if np.isnan(pot_value) or np.isposinf(pot_value):
            last_real = V.shape[0] - i - 1
            if last_real > len(V) - 2:
                raise RuntimeError('Undefined values in tail of potential.'
                                   'This probably means you need better '
                                   'sampling at this state point.')
            return correction_function(r, V, last_real)
        # Retain old potential at small r because:
        #   * current rdf = 0, target rdf > 0 --> -inf values in potential.
        elif np.isneginf(pot_value):
            last_neginf = V.shape[0] - i - 1
            for i, pot_value in enumerate(V[:last_neginf+1]):
                V[i] = previous_V[i]
            return V
    else:
        # TODO: Raise error?
        #       This means that all potential values are well behaved.
        pass


def linear_head_correction(r, V, cutoff):
    """Use a linear function to smoothly force V to a finite value at V(0). """
    slope = ((V[cutoff+1] - V[cutoff+2]) / (r[cutoff+1] - r[cutoff+2]))
    V[:cutoff + 1] = slope * (r[:cutoff + 1] - r[cutoff + 1]) + V[cutoff + 1]
    return V


def exponential_head_correction(r, V, cutoff):
    """Use an exponential function to smoothly force V to a finite value at V(0)

    Parameters
    ----------
    r : np.ndarray
        Separation values
    V : np.ndarray
        Potential at each of the separation values
    cutoff : int
        The last real value of V when iterating backwards

    This function fits the small part of the potential to the form:
    V(r) = A*exp(-Br)
    """
    dr = r[cutoff+2] - r[cutoff+1]
    B = np.log(V[cutoff+1] / V[cutoff+2]) / dr
    A = V[cutoff+1] * np.exp(B * r[cutoff+1])
    V[:cutoff+1] = A * np.exp(-B * r[:cutoff+1])
    return V

def mass_center(m, x, size):
    if len(m) != len(x):  
        print("The length of A and B is not equal")
        return None
    l = len(m)
    d_x = np.abs(np.subtract.outer(x,x).flatten())
    if np.any(d_x > 0.5*size):
        x[x>0.5*size] -= size
    mc = np.sum(m*x)/np.sum(m)
    if mc > size:
        mc -= size
    elif mc < 0:
        mc += size
    return mc 