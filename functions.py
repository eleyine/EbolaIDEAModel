import numpy as np
from scipy import integrate
import plots
import json
import data as Data
from scipy.optimize import curve_fit, leastsq
import collections

CSV = 'data/country_timeseries.csv'
JSON = 'data/raw_data.json'

def cumI(T, R0, d):
    def f(t):
        return I(t, R0, d)

    if isinstance(T, collections.Iterable):
        y = [integrate.quad(f, 0, t)[0] for t in T]
        return np.array(y, dtype=float)
    else:
        return integrate.quad(f, 0, T)[0] # remember to +0.5

def I(t, R0, d):
    '''
    Returns incidence case counts for a given R0 and discount factor d at time t.
    '''
    # print t
    # print R0
    # print d
    # print 'I(%i) ~ (%f, %f)' % (t, R0, d)
    # t = t/7.0
    return (R0 / ((1.0 + d)**t))**t

def RMSD_fit(x, y):
    '''
    Fit R0 and d params of function f to minimise the root-mean-squared differences between
    generation-specific case counts.
    Input
        x: time
        y: cumulative incidence
    Output
        R0
        d
    '''
    def residual(p, x, y):
        return y - cumI(x,*p)

    guess = (1.8, 0.05)
    popt, pcov = leastsq(residual, guess, args=(x, y))
    return popt

def RMSD(y0, y):
    '''
    Return root mean square deviation
    '''
    # print 'y0', y0
    # print 'y', y
    # print 'y0 - y', y0 - y
    # print 'y0 - y**2', ((y0 - y) ** 2)
    # print 'y0 - y**2 mean', ((y0 - y) ** 2).mean()
    return np.sqrt(((y0 - y) ** 2).mean())

def leastsq_fit(x, y):
    '''
    Fit R0 and d params of function f to minimise the root-mean-squared differences between
    generation-specific case counts.
    Input
        x: time
        y: cumulative incidence
    Output
        R0
        d
    '''
    guess = (2.0, 0.1)
    popt, pcov = curve_fit(cumI, x, y, guess)
    return popt
    # pass

def get_t_max(R0, d):
    return np.log(R0) / np.log(1+d)

def _get_I_total(R0, d, t):
    A = np.ln(R0)
    B = np.ln(1+ d) 
    mu = 0.5 * A / B
    I = 0.5 * np.exp( 0.25 * A**2 / B * np.sqrt( np.pi / B))
    I = I * np.sqrt(B) * (np.erf(t - mu)- np.erf(-mu))
    return I


def get_I_total(R0, d):
    def f(t):
        return I(t, R0, d)

    return integrate.quad(f, 0, np.inf)[0]



def get_SIR(t_max, R0, RR, n, N):
    '''
    Returns cumulative incidence count series.
    '''

    x = range(0, t_max)
    y = [None] * len(x)

    # initialise arrays
    Re = [None] * len(x)
    S = [None] * len(x)
    R = [None] * len(x)
    I = [None] * len(x)

    for t in x:
        if t == 0:
            R[t] = 0
            S[t] = N
            I[t] = 1
            Re[t] = R0 / float(N)
            continue
        S[t] = S[t-1] - Re[t-1] * I[t-1]
        Re[t] = R0 * RR**(t**n) * S[t] / float(N)
        I[t] = Re[t-1] * I[t-1]
        R[t] = R[t-1] + I[t]
        print '(t=%i) S=%f, I=%f, Re=%f' % (
            t, S[t], I[t], Re[t])
        y[t] = I[t]
    return (np.array(x), np.array(y))

def test():
    import matplotlib.pyplot as plt

    all_series = Data.get_xy_series(incidence_type='cases') #, fill_time=True)
    x, y = all_series['SierraLeone']
    x, y = Data.to_SI(x, y, 15)
    x = np.array(x)
    y = np.array(y)

    plt.scatter(x, y, label="Real", color="blue")
    R0 = 1.7
    d = 0.38
    y_fit = cumI(x, R0, d)
    plt.scatter(x, y_fit, label="Initial guess", color="red")

    R0, d = RMSD_fit(x, y)
    y_fit = cumI(x, R0, d)
    plt.plot(x, y_fit, label="RMSD Fit (R0=%f, d=%f)" % (R0, d), color="green")

    R0, d = leastsq_fit(x, y)
    y_fit = cumI(x, R0, d)
    plt.plot(x, y_fit, label="Least Squares Fit", color="orange")


    plt.legend(loc=2)
    # plt.show()

    R0, d = (1.7, 0.38)
    print 'RMSD for (%f, %f)' % (R0, d), RMSD(y, cumI(x, R0, d))
    R0, d = (2.221053, 0.022105)
    print 'RMSD for (%f, %f)' % (R0, d), RMSD(y, cumI(x, R0, d))

if __name__ == '__main__':
    # get_SIR(t_max, R0, RR, n, N0):
    import matplotlib.pyplot as plt
    RR_range = np.linspace(1.1, 1.5, 10)

    for RR in RR_range:
        x, y = get_SIR(5, 3.0, 3, 1, 100000)
        plt.title("RR=%f" % (RR))
        plt.plot(x, y)
        plt.show()
        break





