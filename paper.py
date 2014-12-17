import plots as p
import functions as f
import data as Data
import numpy as np
import matplotlib.pyplot as plt
import os

class Analysis:
    def __init__(self,
        incidence_type='cases', 
        countries=['Guinea', 'Liberia', 'SierraLeone'],
        serial_interval=15,
        output_folder='fig/',
        fisman_data=False):

        self.incidence_type = incidence_type
        self.countries = ['Guinea', 'Liberia', 'SierraLeone']
        self.fisman_data=fisman_data
        self.fisman_limit=190
        self.serial_interval=serial_interval
        self.root = output_folder

    # ==========================================================================
    # DATA FUNCTIONS
    # ==========================================================================

    def get_SI_series(self, *args, **kwargs):
        x, y = self.get_series(*args, **kwargs)
        if 'scale' in kwargs.keys():
            x, y = Data.to_SI(x, y, self.serial_interval, scale=kwargs['scale'])
        else:
            x, y = Data.to_SI(x, y, self.serial_interval)
        return (x, y)

    def get_series(self, country, fill_y=False, fill_x=False, 
            skip_none=False, fisman_data=None):

        # special keyword
        if country == 'All Countries':
            return self.get_all_series(fill_y=fill_y, fisman_data=fisman_data)

        if fisman_data is None:
            fisman_data = self.fisman_data

        all_series = Data.get_xy_series(incidence_type=self.incidence_type, 
            countries=[country], 
            fill=fill_y, fill_time=fill_x)
        x, y = all_series[country]

        if fisman_data:
            # find i where x_i > fisman_limit
            i = 0
            while i < len(x) and x[i] < self.fisman_limit:
                i += 1
            x = x[:i]
            y = y[:i]

        if skip_none:
            tuples = [(x_i, y_i) for x_i, y_i in zip(x, y) if y_i is not None]
            x, y = zip(*tuples)
        return (x, y)

    def get_all_series(self, fill_y=False, fisman_data=None):

        if fisman_data is None:
            fisman_data=self.fisman_data

        series = {}
        for i, c in enumerate(self.countries):
            X, Y = self.get_series(c, fill_y=fill_y, fisman_data=fisman_data)
            for x, y in zip(X, Y):
                if x not in series.keys():
                    series[x] = [None] * len(self.countries)
                series[x][i] = y

        to_delete = []
        for k, v in series.iteritems():
            if None in v:
                to_delete.append(k)
            else:
                series[k] = sum(v)
        
        # exclude data points that have at least one missing country 
        for k in to_delete:           
            del series[k]

        x, y = zip(*sorted(series.items()))
        return x, y

    # ==========================================================================
    # MISC HELPERS
    # ==========================================================================

    def get_color_iter(self, n=None):
        if not n:
            n = len(self.countries)
        return Analysis.get_color_map(n)

    def get_colors(self, n=None):
        i = self.get_color_iter(n=n)
        return [c for c in i]

    @staticmethod
    def check_dir(d):
        if not os.path.exists(d):
            os.makedirs(d)

    @staticmethod
    def get_color_map(n):
        import matplotlib.cm as cm
        colors = ['red', 'blue', 'green', 'orange', 'gray', 'purple', 'black']
        if n > len(colors):
            colors = iter(cm.rainbow(np.linspace(0, 1, n)))
        else: 
            colors = iter(colors)
        return colors

    @staticmethod
    def to_latex_table(d, keys, countries=None, caption='Caption'):
        content = ''
        if countries is None:
            countries = sorted(d.keys())
        content += '\\begin{table}[hbt]\n'
        content += '\caption{%s}\n' % (caption)
        content += '\centering\n'
        content += '\\begin{tabular}{l ' + ' '.join(['c' for i in countries]) + '}\n'
        content += '\\toprule \n'
        content += '     \ & ' + ' & '.join([c for c in countries]) +  '\\\\ \n'
        content += '\\midrule \n'
        for k in keys:
            content += '    \\textbf{' + k + '} & ' + ' & '.join(['%.3f' % d[c][k] for c in countries]) + ' \\\\ \n'
        content += '\\bottomrule \n'
        content += '\\end{tabular}\n'
        content += '\\end{table}\n'
        return content


    # ==========================================================================
    # MODEL FITTING
    # ==========================================================================

    def contour_plot(self, out='contour/', save=False):
        countries = self.countries + ['All Countries']
        intervals = [12, 15, 18]
        R0_range = np.linspace(1.4, 3.6, 100)
        d_range = np.linspace(0, 0.07, 100)
        data = {}
        for i, c in enumerate(countries):
            print 'Contour Plot for %s' % (c)
            data[c] = {}
            x, y = self.get_SI_series(c, fill_y=True)
            x = x + 5
            R0, d = f.RMSD_fit(x, y)
            data[c]['R0'] = R0
            data[c]['d'] = d
            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + "param_contour_%s.png" % (c.lower().replace(' ', '_'))
                p.parameter_heatmap(x, y, R0_range, d_range, f.cumI, f.RMSD,
                fn=fn)
            else:
                p.parameter_heatmap(x, y, R0_range, d_range, f.cumI, f.RMSD)
        if save:
            d = self.root + out
            Analysis.check_dir(d)
            keys=['R0', 'd']
            with open(d + 'contour.tex', 'w') as fw:
                fw.write(Analysis.to_latex_table(data, keys, caption='Best-fit R0 and d by Country'))

    # Observed vs Fit Cumulative Incidence, Overall + Country specific
    def observed_vs_model(self, out='simple_fit/', save=False):

        countries = self.countries + ['All Countries']
        data = {}
        for i, c in enumerate(countries):
            data[c] = {}

            # first reported cases assumed to have been reported in generation 5
            start_generation = 5

            fig = plt.figure()
            ax = fig.add_subplot(111)
            width = 0.35

            # up-to-date data
            x, y = self.get_SI_series(c, fill_y=True)
            x = x + 5
            R0, d = f.RMSD_fit(x, y)
            y_rmsd = f.cumI(x, R0, d)
            ax.bar(x, y, width, color='b', alpha=0.3)
            ax.plot(x+width, y_rmsd, color='b', label='Latest Data Fit')

            data[c]['R0'] = R0
            data[c]['d'] = d
            data[c]['RMSD'] = f.RMSD(y, y_rmsd)

            # fisman data
            x, y = self.get_SI_series(c, fill_y=True, fisman_data=True)
            x = x + 5
            R0, d = f.RMSD_fit(x, y)
            y_rmsd = f.cumI(x, R0, d)
            ax.bar(x+width, y, width, color='r', alpha=0.3)
            ax.plot(x+width, y_rmsd, color='r', label='Fisman Data Fit')

            data[c]['Fisman R0'] = R0
            data[c]['Fisman d'] = d
            data[c]['Fisman RMSD'] = f.RMSD(y, y_rmsd)


            plt.legend(loc=2)
            # plt.title(c)
            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'simple_fit_%s.png' % (c.lower().replace(' ', '_'))
                plt.savefig(fn)

                # save latex table
                keys=['R0', 'Fisman R0', 'd', 'Fisman d', 'RMSD', 'Fisman RMSD']
                with open(d + 'simple_fit.tex', 'w') as fw:
                    fw.write(Analysis.to_latex_table(data, keys))
            else:
                plt.show()

    # ==========================================================================
    # PROJECTION
    # ==========================================================================

    def duration_and_size_projection(self, out='projected/', save=False):
        countries = self.countries + ['All Countries']
        data = {}
        for i, c in enumerate(countries):
            data[c] = {}

            # first reported cases assumed to have been reported in generation 5
            start_generation = 5

            # up-to-date data
            x, y = self.get_SI_series(c, fill_y=True)
            x = x + start_generation

            start_generation_index = 2
            generations = range(start_generation_index, len(x))
            I_total = [None] * len(x)
            t_max= [None] * len(x)
            for i in generations:
                _x = x[:i]
                _y = y[:i]
                R0, d = f.RMSD_fit(_x, _y)
                try:
                    I_total[i] = f.get_I_total(R0, d)
                except:
                    I_total[i] = None
                if I_total[i] > 1e6:
                    I_total[i] = None
                t_max[i] = f.get_t_max(R0, d)

            colors = ['red', 'blue', 'green']
            series = [t_max, I_total, y]
            labels = ['Projected Outbreak Duration', 'Projected Outbreak Size', 'Cumulative Incidence']

            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Generations Available for Fitting')
            ax2 = ax1.twinx()
            # ax2.set_ylim(0, 200000)

            ax3 = ax1.twinx()
            ax3.spines['right'].set_position(('outward', 55))
            plt.subplots_adjust(right=0.8)
            axes = [ax1, ax2, ax3]

            width=1
            for i, s in enumerate(series):
                ax = axes[i]
                if i == 2:
                    ax.bar(x-width*0.5, s, width, color=colors[i], alpha=0.3)
                else:
                    ax.plot(x, s, color=colors[i])

                    scatter_tuples = [ (n, m) for n, m in zip(x, s) if m is not None]
                    xS, sS = zip(*scatter_tuples)
                    print xS, sS
                    ax.scatter(xS, sS, color=colors[i], alpha=0.5)
                ax.set_ylabel(labels[i], color=colors[i])
                for tl in ax.get_yticklabels():
                    tl.set_color(colors[i])

            # p.align_yaxis(ax1, 0, ax2, 0)
            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'projected_size_and_duration_%s.png' % (c.lower().replace(' ', '_'))
                plt.savefig(fn)
            else:
                plt.show()

            print c
            print 'I total:', I_total[-1]
            print 't max:', t_max[-1]

    # ==========================================================================
    # CONTROL ASSESSMENT
    # ==========================================================================

    def control_assessment(self, out='control/', save=False):
        countries = self.countries + ['All Countries']
        data = {}
        for i, c in enumerate(countries):
            data[c] = {}

            # first reported cases assumed to have been reported in generation 5
            start_generation = 5

            # up-to-date data
            x0, y0 = self.get_SI_series(c, fill_y=True)
            x0 = x0 + start_generation

            fig, (ax1, ax2) = plt.subplots(1, 2)#, sharey=True)
            fig.set_size_inches(15,5)
            #.plot(x, y) # projection
            # axarr[0].set_title('Sharing X axis')
            # axarr[1].scatter(x, y) # error

            for i in range(2, len(x0)):
                x = x0[0:i]
                y = y0[0:i]

                actual = y0[i]

                R0, d = f.RMSD_fit(x, y)
                SI = x[-1]
                projected = f.cumI(SI+1, R0, d)

                # projection vs actual
                if projected < 2 * actual:
                    ax1.scatter(actual, projected, alpha=0.3)
                    if actual > 500:
                        ax1.annotate('%i' % (SI), (actual-100, projected+150))


                # error
                error = np.abs(projected-actual) / float(actual)
                if error > 1.5:
                    print 'Error for SI=%i is %f !' % (x[-1], error)
                else:
                    plt.scatter(x0[i], error)

            ax1.plot(y0, y0, 'b-')
            ax1.set_ylabel('Projected Cases')
            ax1.set_xlabel('Actual Cases')
            ax2.set_ylim([0, 1])
            ax2.set_ylabel('Percent Error')
            ax2.set_xlabel('Serial Inteval')

            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'control_%s.png' % (c.lower().replace(' ', '_'))
                plt.savefig(fn)
            else:
                plt.show()

    # ==========================================================================
    # MULTI-WAVE EPIDEMICS
    # ==========================================================================
  
    def deltaD(self, out='deltaD/', save=False):
        countries = self.countries + ['All Countries']
        data = {}
        for i, c in enumerate(countries):
            data[c] = {}

            # first reported cases assumed to have been reported in generation 5
            start_generation = 5

            # up-to-date data
            x0, y0 = self.get_SI_series(c, fill_y=True)
            x0 = x0 + start_generation


            fig, ax1 = plt.subplots()

            intervals = range(5, len(x0))
            R0_list = [None] * len(x0)
            d_list = [None] * len(x0)
            deltaD_list = [None] * len(x0)

            colors = self.get_color_map(n=len(intervals))  
            prev_d = None 
            for i in intervals:     
                x = x0[:i]
                y = y0[:i] 
                R0, d = f.RMSD_fit(x, y)
                R0_list[i] = R0
                d_list[i] = d

                if prev_d is None:
                    prev_d = d
                else:
                    deltaD = d - prev_d
                    deltaD_list[i] = deltaD
                    prev_d = d

            ax1.set_xlabel("Generations Available For fitting")

            # plot d
            ax1.plot(x0, d_list, 'b--', label='d')
            d_tuples = [i for i in zip(x0, d_list) if i[1] is not None]
            ax1.scatter(*zip(*d_tuples), alpha=0.5)
            ax1.plot(x0, deltaD_list, 'b-', label='delta d')
            ax1.set_ylabel('', color='b')
            for tl in ax1.get_yticklabels():
                tl.set_color('b')

            # plot deltaD
            deltaD_tuples = [i for i in zip(x0, deltaD_list) if i[1] is not None]
            _x, _y = zip(*deltaD_tuples)
            ax1.fill_between(_x, 1e-6, _y, facecolor='blue', alpha=0.5)
            plt.legend(loc=2)

            # plot R0
            ax2 = ax1.twinx()
            # ax2.set_yscale('log')
            ax2.plot(x0, R0_list, 'r-', label='R0')
            R0_tuples = [i for i in zip(x0, R0_list) if i[1] is not None]
            ax2.scatter(*zip(*R0_tuples), color='red', alpha=0.5)
            ax2.set_ylabel('', color='r')
            for t2 in ax2.get_yticklabels():
                t2.set_color('r')
            plt.legend()

            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'deltaD_%s.png' % (c.lower().replace(' ', '_'))
                plt.savefig(fn)
            else:
                plt.show()          
    
    # ==========================================================================
    # SENSITIVITY ANALYSIS
    # ==========================================================================
    def progressive_projections(self, out='sensitivity/', save=False):
        countries = self.countries + ['All Countries']
        data = {}
        for i, c in enumerate(countries):
            data[c] = {}

            # first reported cases assumed to have been reported in generation 5
            start_generation = 5

            # up-to-date data
            x0, y0 = self.get_SI_series(c, fill_y=True)
            x0 = x0 + start_generation

            intervals = range(2, len(x0))

            plt.figure()
            fig = plt.figure()
            ax = plt.subplot(111)

            colors = self.get_colors(n=len(intervals))
            for ind, i in enumerate(intervals):
                x = x0[:i+1]
                y = y0[:i+1]
                R0, d = f.RMSD_fit(x, y)

                y_fit = f.cumI(x0, R0, d)
                ax.plot(x, y_fit[:i+1], color=colors[ind], label='SI=%i' % (i))
                ax.plot(x0, y_fit, '--', color=colors[ind]) #, ls='--')
                ax.scatter(x[-1], y_fit[i], color=colors[ind])

            width = 0.8
            ax.bar(x0-0.5*width, y0, width, color='blue', alpha=0.1)
            ax.set_ylim(0, y0[-1]*1.5)
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
            ncol=5, fancybox=True, shadow=True)
            ax.set_ylabel('Cumulative Incidence Count')
            ax.set_xlabel('Generations Available for Fitting')
            
            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'progressive_projections_%s.png' % (c.lower().replace(' ', '_'))
                plt.savefig(fn)
            else:
                plt.show()

    def sensitivity_analysis(self, out='sensitivity/', save=False):
        countries = self.countries + ['All Countries']
        keys = [
            'Base Case', 
            '12 day generation time', 
            '18 day generation time',
            'Outbreak recognized generation 3',
            'Outbreak recognized generation 7',
            'Outbreak 50% underreported',
            'Outbreak 99% underreported',
            'Deaths only'
            ]

        for c in countries:
            series = {'up-to-date': {}, 'Fisman': {}}

            for k, xy in series.iteritems():

                if k == 'Fisman':
                    self.fisman_data = True
                else:
                    self.fisman_data = False

                x0, y0 = self.get_SI_series(c, fill_y=True)
                xy['Base Case'] = (x0+5, y0)

                self.serial_interval = 12
                x0, y0 = self.get_SI_series(c, fill_y=True)
                xy['12 day generation time'] = (x0+5, y0)

                self.serial_interval = 18
                x0, y0 = self.get_SI_series(c, fill_y=True)
                xy['18 day generation time'] = (x0+5, y0)

                #reset
                self.serial_interval = 15

                x0, y0 = self.get_SI_series(c, fill_y=True)
                xy['Outbreak recognized generation 3'] = (x0+3, y0)
                xy['Outbreak recognized generation 7'] = (x0+7, y0)

                xy['Outbreak 50% underreported'] = ((x0+5)*0.5, y0)
                xy['Outbreak 99% underreported'] = ((x0+5)*0.01, y0)

                self.incidence_type = 'deaths'
                x0, y0 = self.get_SI_series(c, fill_y=True)
                xy['Deaths only'] = (x0+5, y0)

                # reset
                self.incidence_type = 'cases'

            data = {}
            categories = ['R0', 'R0 (Fisman)', 'd', 'd (Fisman)']
            for category in categories:
                if 'Fisman' in category:
                    sub_series = series['Fisman']
                else:
                    sub_series = series['up-to-date']

                column = {}
                for k in keys:
                    x, y = sub_series[k]
                    R0, d = f.RMSD_fit(x, y)

                    if 'R0' in category:
                        column[k] = R0
                    else:
                        column[k] = d
                data[category] = column
            table = Analysis.to_latex_table(data, keys, countries=categories, 
                caption=c)
            if save:
                d = self.root + out
                Analysis.check_dir(d)
                fn = d + 'sensitivity_%s.tex' % (c.lower().replace(' ', '_'))
                with open(fn, 'w') as fw:
                    fw.write(table)
            else:
                print table

    # ==========================================================================
    # APPENDIX
    # ==========================================================================

    def demonstrate_extrapolation(self, out='extrapolation/', save=False):
        '''
        Observed vs Extrapolated Cumulative Incidence, Overall + Country specific.
        Local linear extrapolation points are circles, real data points are disks.
        '''
        plt.figure()
        ALPHA = 0.3
        colors = self.get_colors(n=len(self.countries)+1)
        for i, c in enumerate(self.countries):
            # plot observed data
            oX, oY = self.get_series(c, fill_y=False, skip_none=True)
            plt.scatter(oX, oY, color=colors[i], alpha=ALPHA)
            # plot interpolated data
            iX, iY = self.get_series(c, fill_y=True)
            for x, y in zip(iX, iY):
                if x not in oX:
                    plt.scatter(x, y, color=colors[i], facecolors='none')
            plt.plot(iX, iY, label=c, color=colors[i], alpha=ALPHA)

        # add all countries data
        oX, oY = self.get_all_series(fill_y=False)
        plt.scatter(oX, oY, color=colors[-1], alpha=ALPHA)
        iX, iY = self.get_all_series(fill_y=True)
        plt.plot(iX, iY, label='All Countries', color=colors[-1], alpha=ALPHA)
        for x, y in zip(iX, iY):
            if x not in oX:
                plt.scatter(x, y, color=colors[-1], facecolors='none')

        plt.legend(loc=2)
        plt.xlabel("Days (starting from 03/22/2014)")
        plt.ylabel("Cumulative Incidence")
        if save:
            d = self.root + out
            Analysis.check_dir(d)
            plt.savefig(d + 'demo.png')
        else:
            plt.show()

    def demonstrate_integration(self, R0=3.91, d=0.12, out='integration/', save=False):
        '''
        Compares simulated data for a given R0 and d with integration.
        '''
        import functions as f

        x = np.arange(0, 30)
        y = np.array([f.I(i, R0, d) for i in x], dtype=float)
        y_int = [f.cumI(i, R0, d) for i in x]

        import matplotlib.pyplot as plt
        plt.plot(x, y, label="Non-Cumulative", color='green')
        plt.plot(x, np.cumsum(y), label="Cumulative", color='red')
        plt.plot(x, y_int, label="Integral", color='blue', lw=10, alpha=0.3)

        plt.legend(loc=4)
        if save:
            d = self.root + out
            Analysis.check_dir(d)
            plt.savefig(d + 'demo.png')
        else:
            plt.show()


if __name__ == '__main__':
    a = Analysis(incidence_type='cases')
    # a.contour_plot(save=True)
    # a.observed_vs_model(save=True)
    a.duration_and_size_projection()
    # a.control_assessment(save=True)
    # a.deltaD(save=True)
    # a.progressive_projections(save=True)
    # a.sensitivity_analysis(save=True)
    # a.demonstrate_extrapolation()
    # a.demonstrate_integration()