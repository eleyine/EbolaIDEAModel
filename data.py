import numpy as np
import json

CSV = 'data/country_timeseries.csv'
JSON = 'data/raw_data.json'

def write_raw_json(output=JSON, is_reversed=True):
    '''
    Convert CSV to json.
    '''
    def int_or_None(i):
        if i == '':
            return None
        elif i is None:
            return None
        else:
            return int(i)

    with open(CSV, 'r') as f:
        lines = f.readlines()
    header = [i.strip() for i in lines[0].split(',')]
    content = lines[1:]

    csv_d = {}
    for row in content:
        row = row.strip()
        for i, t in enumerate(row.split(',')):
            h = header[i]  
            if h not in csv_d.keys():
                csv_d[h] = []
            csv_d[h].append(t)

    time = {'date': [], 'day': []}
    country_data = {}
    for k, v in csv_d.iteritems():
        if k in ['Date', 'Day']:
            if k == 'Day':
                v = [int_or_None(i) for i in v]
            if is_reversed:
                v.reverse()
            time[k.lower()] = v
        else:
            incidence, country = k.split('_')
            if country not in country_data.keys():
                country_data[country] = {
                    'cases': [],
                    'deaths': []
                }
            l = [int_or_None(i) for i in v]
            if is_reversed:
                l.reverse()
            country_data[country][incidence.lower()] = l
    data = {'country_data': country_data, 'time': time}

    with open(output, 'w') as f:
        json.dump(data, f, indent=2)

def _get_data(cumulative=True, use_date=False, fill=False, fill_time=False, 
    countries=None):

    with open(JSON, 'r') as f:
        data = json.load(f)

    # filter countrydata used
    if countries:
        for c in data['country_data'].keys():
            if c not in countries:
                del data['country_data'][c]

    if use_date:
        # convert to date objects
        import datetime
        for i, d in enumerate(data['time']['date']):
            data['time']['date'][i] = datetime.datetime.strptime(
                d, "%m/%d/%Y").date()
    else:
        del data['time']['date']

    if fill_time:
        new_data = {
            'time': {'day': []},
            'country_data': {}
        }
        # make new data dictionary where time series has no gaps
        time = data['time']['day']
        new_time = range(0, data['time']['day'][-1]+1) 
        # TODO: do the same for dates
        country_data = {}
        for c, cdata in data['country_data'].iteritems():
            country_data[c] = {
                'cases': [None] * len(new_time), 
                'deaths': [None] * len(new_time)
            }
            for inc, series in cdata.iteritems():
                for i, x in enumerate(series):
                    t = time[i]
                    country_data[c][inc][t] = x
        data = {'time': {'day': new_time}, 'country_data': country_data}

    if fill:
        for country, country_data in data['country_data'].iteritems():
            for inc, series in country_data.iteritems():
                for i, x in enumerate(series):
                    if x is None:
                        if i == 0:
                            series[i] = 0
                        else:
                            # get index of previous non-null element
                            prev_i = i
                            while series[prev_i] is None and prev_i > 0:
                                prev_i -= 1
                                # print '(prev=(%i, %f)' % (
                                #     prev_i, series[prev_i] )

                            # get index of next non-null element
                            next_i = i
                            while next_i < len(series) and series[next_i] is None:
                                next_i += 1

                            if next_i >= len(series):
                                series[i] = series[i-1]
                                continue
                                # print '(next=(%i, %f)' % (
                                #     next_i, series[next_i] if not None else 0 )
                            p = series[prev_i]; n = series[next_i]
                            if p and n:
                                series[i] = p + (n - p)/ float(next_i - prev_i) * (i - prev_i)
                                # print 'x=%i, y=%f (prev=(%i, %f), next=(%i, %f)' % (
                                #     i, series[i], prev_i, p, next_i, n)
                            else:
                                series[i] = series[i-1]
    return data

# def get_xy_tuples(cumulative=True, use_date=False, fill=True, fill_time=False, incidence_type='cases'):
#     series = get_xy_series(
#         cumulative=cumulative, 
#         use_date=use_date, 
#         fill=fill, 
#         fill_time=fill_time, 
#         incidence_type=incidence_type)

#     tuples = {}
#     for c, s in series.iteritems():
#         tuples[c] = zip(*s)
#     return tuples

def get_xy_series(cumulative=True, 
    use_date=False, 
    fill=False, 
    fill_time=False, 
    incidence_type='cases',
    countries=None):
    '''
    Read csv and return dictionary of the form:
    'country_data':
        <country>:
            'cases': list of int
            'deaths': list of int
    'time':
        'date': list of Dates (if use_date is True)
        'day': list of int
    '''
    data = _get_data(cumulative=cumulative, use_date=use_date, fill=fill, 
        fill_time=fill_time, countries=countries)

    if use_date:
        time = data['time']['date']
    else:
        time = data['time']['day']

    series = {}
    all_y = [0] * len(time)
    for c, cdata in data['country_data'].iteritems():
        series[c] = [time, cdata[incidence_type]]
        for i, y in enumerate(series[c][1]):
            if y:
                all_y[i] += y

    # now add extra category 'all' (only works if fill was selected)
    if fill:
        series['All Countries'] = [time, all_y]
    return series

def to_SI(x, y, si, is_cumulative=True, scale=False):
    # first make sure x has no gaps
    new_x = range(0, x[-1]+1)
    new_y = [None] * len(new_x)
    for x_i, y_i in zip(x, y):
        new_y[x_i] = y_i
    x = new_x; y = new_y

    # parition arrays according to si
    sub_x = np.array(x[::si])
    sub_y = np.array(y[::si])
    # print 'Y', sub_y
    for sub_i, y_i in enumerate(sub_y):
        # try to substitute by any value that precedes the ith
        # interval, but is within the interval
        if y_i is None:
            # print '    %ith position is None' % (sub_i)
            # substitute by earliest non-null value
            i = sub_i * si
            for i in reversed(range(i - si, i+1)):
                if i <= 0:
                    y_i = 0
                    break
                if y[i] is not None:
                    # print '      > sub by y[%i] = %i' % (i, y[i])
                    y_i = y[i]
                    break

        # if this fails, try to substitute by any preceding value
        if y_i is None:
            i = sub_i * si
            while i > 0 and y[i] is None:
                i -= 1
            y_i = y[i]
        # if y_i still is None, I don't know what you did wrong
        if y_i is None:
            raise Exception('wat')
        sub_y[sub_i] = y_i

    if not scale:
        sub_x = np.array(sub_x / si, dtype=float)
    else:
        sub_x = np.array(sub_x, dtype=float)
    sub_y = np.array(sub_y, dtype=float)
    return (sub_x, sub_y)

def to_tuples(x, y):
    return zip(*s)

if __name__ == '__main__':
    write_raw_json()
    # print _get_data(fill=True, fill_time=True)
    x, y = get_xy_series()['Liberia']
    # print 'x', x
    print 'y', y
    x, y = to_SI_interval(x, y, 7)
    # print 'x', x
    print 'y', y
    # print get_xy_tuples()['All Countries']

