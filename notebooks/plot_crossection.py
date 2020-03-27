import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import pyart
from pyart.io.sband_archive import read_sband_archive
from pyart.graph.gridmapdisplay import GridMapDisplay


def get_radar_info(filename, radar_info):
    """根据雷达名称找雷达的经纬度信息
    :param filename(str): 雷达数据文件名
    :param radar_info(json): 存储雷达站点信息数据
    :return lat(deg), lon(deg), elev(m)
    """
    import os

    name = os.path.basename(filename)
    station_id = [int(name[idx:idx+4]) for idx in range(len(name)-4) if name[idx:idx+4].isdigit()][0]

    if station_id not in radar_info.index:
        raise ValueError(f'没有检索到相关雷达站点信息! {station_id}')

    return radar_info.loc[station_id, 'Latitude'], radar_info.loc[station_id, 'Longitude'], radar_info.loc[station_id, 'Elevation']

if __name__ == '__main__':

    # 雷达文件和底图
    filename = 'Z_RADR_I_Z9515_20160623063100_O_DOR_SA_CAP.bin'
    #shpfn = 'data/shps/cnhimap.shp'

    # 读取雷达反射率
    sweep = 0    # 体扫仰角序号
    field = 'reflectivity'
    radar_info = pd.read_json('radar_info.json')

    radar = read_sband_archive(filename, station=get_radar_info(filename, radar_info))

    x = (-400000.0, 400000.0)
    y = (-400000.0, 400000.0)
    z = (10, 20000)
    gshape = (50, 801, 801)
    # exclude masked gates from the gridding
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_masked(field)

    grid = pyart.map.grid_from_radars((radar,),
                                    gatefilters=(gatefilter,),
                                    grid_shape=gshape,
                                    grid_limits=(z, y, x),
                                    fields=[field]
                                    )

    display = GridMapDisplay(grid)

    fig, ax = plt.subplots(figsize=(12, 9))

    display.plot_latlon_slice('reflectivity', coord1=(119.7520, 33.6601), coord2=(120.0182, 33.8739), cmap='pyart_NWSRef', vmin=0, vmax=70)

    plt.show()
