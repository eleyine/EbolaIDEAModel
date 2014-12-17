import matplotlib.pyplot as plt
import numpy as np


OUT = 'fig/'

def real_series(x, y, country, title=None):
    '''
    x: date 
    y: cumulative infections 
    trace for various countries
    '''
    plt.figure()
    plt.scatter(x, y)
    plt.plot(x, y)
    if title is None:
        title = "%s Ebola Cases vs Time" % (country)
    plt.title(title)
    plt.xlabel("Days (starting from 3/22/2014)")
    plt.ylabel("Cumulative Incidence")
    fn = OUT + 'real/' + title.lower().replace(' ', '_') + '.png'
    plt.savefig(fn)

def aggregate_real_series(all_series, title=None):
    '''
    x: date 
    y: cumulative infections 
    trace for various countries
    '''
    plt.figure()
    colors = get_color_map(len(all_series.keys()))
    for country, series in all_series.iteritems():
        x, y = series
        color = next(colors)
        plt.scatter(x, y, color=color, alpha=0.3)
        plt.plot(x, y, color=color, label=country)
    if title is None:
        title = "Ebola Cases vs Time"
    plt.title(title)
    plt.legend(loc=2)
    plt.xlabel("Days (starting from 3/22/2014)")
    plt.ylabel("Cumulative Incidence")
    fn = OUT + 'real/' + 'aggregate.png'
    plt.savefig(fn)


def parameter_heatmap(x, y, R0, d, model, error, title=None, fn=None):
    '''
    Constructs heatmap of RMSD obtained via error(y, y_fit) for 
    different values of R0 and d
    '''
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import LogFormatter


    # get RMSD
    Z = np.zeros((len(R0), len(d)))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            # print 'R0', R0[i]
            # print 'd', d[j]
            y_fit = model(x, R0[i], d[j])
            Z[j][i] = error(y, y_fit) # np.log10(error(y, y_fit))
    # print Z

    fig=plt.figure()
    ax1 = fig.add_subplot(111)
    # print 'Z max', Z.max()
    # print int(np.log10(Z.max()))
    max_logZ = np.ceil(np.log10(Z.max()))
    lvls = np.logspace(0,max_logZ,30)
    CF = ax1.contourf(R0, d, Z,
         norm = LogNorm(),
         levels = lvls
        )
    l_f = LogFormatter(10, labelOnlyBase=False)
    cbar = plt.colorbar(CF, ticks=lvls, format=l_f)
    # lvls = np.linspace(0, np.ceil(Z.max()), 20)
    # plt.contourf(R0, d, Z)
    # plt.colorbar()
    # lvls = np.linspace(0, np.ceil(Z.max()), 20)
    # plt.contourf(R0, d, Z,levels=lvls)
    # plt.colorbar(ticks=lvls)



    plt.ylabel('Control Parameter (d)')
    plt.xlabel('Basic Reproductive Number (R0)')


    # plot the fitted R0 and d values 
    from functions import RMSD_fit, leastsq_fit
    R0, d = RMSD_fit(x, y)
    plt.scatter(R0, d, c='yellow', marker="*", s=80)
    if title:
        # title += 'R0=%f, d=%f' % (R0, d)
        plt.title(title)
    if fn:
        plt.savefig(fn)
    else:
        plt.show()


def order_control_vs_RMSD():
    '''
    x: order of control [0, 3]
    y: Normalized RMSD
    trace curves with varying R0
    '''
    pass

def idea_vs_SIR():
    '''
    Trace cumulative infections SIR vs IDEA for given 
    R0, d and order of control
    x: number of generations
    y: cumulative infections
    trace curves for SIR and IDEA
    '''
    pass

def case_IDEA():
    '''
    x: fraction of cases identified
    y: Normalized RMSD
    trace curves for varying R0
    '''
    pass

def projection_evaluation():
    '''
    x: actual case 
    y: projected case
    trace y = x
    at every t, get projected case based on [0, t-1] 
    and compare with actual case
    '''
    pass

def deltaD():
    '''
    x: t
    y: d, deltaD
    '''
    pass

def interval_effect():
    '''
    x: t
    y: cumulative incidence IDEA fit
    Trace curves with different SI
    '''
    pass

# helpers

def get_color_map(n):
    import matplotlib.cm as cm
    colors = ['red', 'blue', 'green', 'orange', 'gray', 'purple', 'black']
    if n > len(colors):
        colors = iter(cm.rainbow(np.linspace(0, 1, n)))
    else: 
        colors = iter(colors)
    return colors

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2,(y1-y2)/2,v2)
    adjust_yaxis(ax1,(y2-y1)/2,v1)

def adjust_yaxis(ax,ydif,v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny>maxy or (-miny==maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)
    
if __name__ == '__main__':
    import functions as f
    import numpy as np
    import data as Data

    R0_range = np.linspace(1.4, 2.6, 10)
    d_range = np.linspace(0, 0.07, 10)
    all_series = Data.get_xy_series(incidence_type='cases') #, fill_time=True)
    x, y = all_series['SierraLeone']
    x, y = Data.to_SI(x, y, 15)

    parameter_heatmap(x, y, R0_range, d_range, f.cumI, f.RMSD)