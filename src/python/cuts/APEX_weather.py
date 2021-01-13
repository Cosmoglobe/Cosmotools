#! /usr/bin/env python

import re, time, os, sys
import urllib
URL = 'http://archive.eso.org/wdb/wdb/eso/meteo_apex/query'
SKIP_KEYWORDS = ['<', '>', 'alt', 'nbsp']


def getMonth(key):
    month_map = {'jan':  1, 'feb':  2, 'mar':  3, 'apr':  4, 'may':  5, 'jun':  6, 'jul':  7, 'aug':  8, 'sep':  9, 'oct': 10, 'nov': 11, 'dec': 12}
    key = key.lower()
    if month_map.has_key(key): return month_map[key]
    return 0

APEX_weather_keyList = ['PWV', 'dew_point', 'humidity', 'pressure', 'temperature', 'wind_direction', 'wind_speed']


APEX_weather_unit = {
    'PWV': 'mm',
    'dew_point': 'C',
    'humidity': '%', 
    'pressure': 'hPa', 
    'temperature': 'C', 
    'wind_direction': 'deg', 
    'wind_speed': 'm/s'
    }

def APEX_weather(begin_date, end_date):
    params = urllib.urlencode({
            'start_date': '%s..%s' %(begin_date, end_date),
            'tab_start_date': True,
            'tab_pwv': True, 
            #'pwv': '0.0..100.0'
            'tab_dewpoint': True,
            #'dewpoint': '-273.0..100',
            'tab_humidity': True,
            'tab_pressure': True,
            'tab_temperature': True,
            'tab_winddirection': True,
            'tab_windspeed': True,
            'tab_export': True,
            'force_tabular_mode': True
            })
    
    fr = urllib.urlopen(URL, params)

    APEX_weatherList = []
    for line in fr:
        is_it_skip = False
        for skip_keyword in SKIP_KEYWORDS: 
            if not re.search(skip_keyword, line): continue
            is_it_skip = True
            pass
        if is_it_skip: continue
        fieldList = line.strip().split()
        if not len(fieldList) >0: continue
        if re.search('N/A', line): 
            #print line.strip()
            continue

        # JZ updated to test for malformed query
        if not fieldList[2].isdigit():
            print 'The variable contains no digit'
            print 'f2',fieldList[2]
            print fieldList
            return {}

        time_str = '%4d.%02d.%02d %s %s' %(int(fieldList[2]), getMonth(fieldList[0]), int(fieldList[1]), fieldList[3][0:len(fieldList[3])-6], fieldList[3][-2:len(fieldList[3])]) 
        tm_tuple = time.strptime(time_str, '%Y.%m.%d %I:%M:%S %p')
        #date_str = time.strftime('%Y-%m-%d %H:%M:%S', tm_tuple)
        pwv = float(fieldList[4]) # mm
        dew_point = float(fieldList[5]) # C
        humidity = float(fieldList[6]) # %
        pressure = float(fieldList[7]) # hPa
        temperature = float(fieldList[8]) # C
        wind_direction = float(fieldList[9]) # deg
        wind_speed = float(fieldList[10]) # m/s
        APEX_weatherList.append( 
            {'tm_tuple': tm_tuple,
             'PWV': pwv,
             'dew_point': dew_point,
             'humidity': humidity,
             'pressure': pressure,
             'temperature': temperature,
             'wind_direction': wind_direction,
             'wind_speed': wind_speed
             }
            )
        pass

    return APEX_weatherList

if __name__ == '__main__':
    dataList = APEX_weather(sys.argv[1], sys.argv[2])

    line = '# Time (UTC)'
    for key in APEX_weather_keyList:    
        line += '  %s (%s)' %(key, APEX_weather_unit[key])
        pass
    print line

    for data in dataList:
        line = time.strftime('%Y-%m-%d %H:%M:%S', data['tm_tuple'])
        for key in APEX_weather_keyList:
            #line += '\t%6.2f [%s]' %(data[key], APEX_weather_unit[key])
            line += '\t%6.2f' %(data[key])
            pass
        print line
        pass
    pass
