# coding:utf-8
import math

def add_province(shpfn, ax):
    """添加底图
    :param shpfn(str): 底图文件名
    :param ax(matplotlib.axes): 坐标轴实例对象
    :return ax(matplotlib.axes): 坐标轴实例对象
    """
    import cartopy.crs as ccrs
    import cartopy.io.shapereader as shpreader
    import cartopy.feature as cfeature

    reader = shpreader.Reader(shpfn)
    states_provinces = cfeature.ShapelyFeature(reader.geometries(),
                                               crs=ccrs.PlateCarree(),
                                               edgecolor='face',
                                               facecolor='None')

    ax.add_feature(states_provinces, edgecolor='black')

    return ax

def auto_ticks(ax):
    """简单的自动设置x和y轴的刻度
    :param ax(matplotlib.axes): 坐标轴实例对象
    :return ax(matplotlib.axes): 坐标轴实例对象
    """
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    lon_formatter = LongitudeFormatter(#number_format='.1f',
                                       degree_symbol='',
                                       dateline_direction_label=True)

    lat_formatter = LatitudeFormatter(#number_format='.1f',
                                      degree_symbol='')

    xlmin, xlmax = ax.get_xlim()
    ylmin, ylmax = ax.get_ylim()
    xmax, xmin, xitv = _itv(xlmax, xlmin)
    ymax, ymin, yitv = _itv(ylmax, ylmin)

    ax.set_xlim([xlmin, xlmax])
    ax.set_ylim([ylmin, ylmax])
    ax.set_xticks(np.arange(xmin, xmax+xitv, xitv))
    ax.set_yticks(np.arange(ymin, ymax+yitv, yitv))

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    return ax

def _itv(llmax, llmin):
    """获取x和y轴的的范围及刻度间隔
    :param llmax(float或int):坐标轴范围的最大值
    :param llmin(float或int):坐标轴范围的最小值
    :return llmax, llmin, itv
    """
    itv = (llmax - llmin)/5

    if itv < 1:
        if itv >= 0.5:
            itv = round(0.5, 1)
            return int(llmax*10)//int(itv*10)*int(itv*10)/10+itv, int(llmin*10)//int(itv*10)*int(itv*10)/10+itv, itv
        else:
            itv = round(itv, 2)
            return int(llmax*100)//int(itv*100)*int(itv*100)/100+itv, int(llmin*100)//int(itv*100)*int(itv*100)/100+itv, itv
    else:
        return math.floor(llmax), math.ceil(llmin), round(itv, 0)

def retrive_wind(radar, sweep_number=1, filename=None, lonlat=False, kwargs=None):
    """单多普勒雷达风场反演
    :param radar(pyart.): 存储雷达数据的实例对象
    :param sweep_number(int): 径向风的体扫序号
    :param filename(str): 用于保存风场反演结果的文件名，保存为netCDF格式
    :param lonlat(bool): 是否获取反演风场结果对应的经纬度坐标，默认为False，返回到雷达的距离
    :param kwargs(dict): singledop.SingleDoppler2D的可选参数，具体见SingleDoppler2D的函数说明
    :return display, x, y 或 sds(singledop.singledop.NetcdfSave)
    """

    if kwargs is None:
        kwargs = dict(
                      L=30,
                      sweep_number=sweep_number,
                      name_vr='velocity',
                      thin_factor=[4, 4],
                      grid_edge=120,
                      )

    retwind = singledop.SingleDoppler2D(radar, **kwargs)
    display = singledop.AnalysisDisplay(retwind)

    if lonlat:
        x, y = get_lon_lat(radar, display)
    else:
        x, y = display.analysis_x.view(), display.analysis_y.view()

    if filename:
        sds = singledop.NetcdfSave(retwind, filename, radar=radar)

    return display, x, y

def get_lon_lat(radar, display=None, x=None, y=None):
    """获取反演风场对应的经纬度坐标
    :param radar(pyart): 读取的雷达数据
    :param display(SingleDop): 风场反演结果
    :param x(float): x方向到雷达的距离, 单位 km
    :param y(float): y方向到雷达的距离, 单位 km
    :return longitude(np.array), latitude(np.array)

    display.analysis_x 和 display.analysis_y 存储为到雷达的距离，单位为 km
    """
    from pyart.core.transforms import cartesian_to_geographic

    if display is not None:
        x = display.analysis_x.view()
        y = display.analysis_y.view()
    elif x is not None and y is not None:
        x, y = x, y
    else:
        raise ValueError(f'display {display} 和 x {x}, y {y} 不能同时为None!')

    projparams = radar.projection.copy()
    if projparams.pop('_include_lon_0_lat_0', False):
        projparams['lon_0'] = radar.longitude.get('data')[0]
        projparams['lat_0'] = radar.latitude.get('data')[0]

    geographic_coords = cartesian_to_geographic(x*1000,
                                                y*1000,
                                                projparams)

    return geographic_coords[0], geographic_coords[1]

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

def plot(radar, field, sweep, display, vellon, vellat, quiver=False, savefig=None, center='radar'):
    """叠加雷达基本反射率和反演风场
    :param radar(pyart.core.radar.Radar): 存储的雷达数据
    :param field(str): 要读取的雷达产品名
    :param sweep(int): 体扫中的仰角序号
    :param display(singledop.singledop.AnalysisDisplay或singledop.singledop.NetcdfSave): 风场反演结果
    :param vellon(np.array): 风场反演对应的经度坐标
    :param vellat(np.array): 风场反演对应的维度坐标
    :param quiver(bool): 控制输出quiver还是barbs, False 为 quiver, True 为 barbs
    :param savefig(str): 存储图片名
    :param center(str): 选择图的中心, 默认为以雷达为中心, 否则为默认绘图无中心
    """
    from pyart.graph import common
    title = common.generate_title(radar, field, sweep)

    if center == 'radar':
        x = (-200000.0, 200000.0)
        y = (-200000.0, 200000.0)
        z = (2000, 2000)
        gshape = (1, 401, 401)
        # exclude masked gates from the gridding
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_masked(field)

        # perform Cartesian mapping, limit to the reflectivity field.
        grid = pyart.map.grid_from_radars(
                                          (radar,), gatefilters=(gatefilter,),
                                          grid_shape=gshape,
                                          grid_limits=(z, y, x),
                                          fields=[field]
                                         )

        ref = grid.fields[field]['data'][sweep]
        lon, lat = grid.get_point_longitude_latitude()
    else:
        ref = radar.get_field(sweep, field)
        # 获取雷达经纬度
        lat, lon, alt = radar.get_gate_lat_lon_alt(sweep)

    # 筛选需要显示的风场数据
    thin = 10
    cond = np.logical_and(display.analysis_x % thin == 0,
                          display.analysis_y % thin == 0)

    # 画图
    fig, ax = plt.subplots(figsize=(12, 9), subplot_kw=dict(projection=ccrs.PlateCarree()))

    # 雷达基本反射率
    ax.pcolormesh(lon, lat, ref, vmin=0, vmax=70, cmap='pyart_NWSRef')

    # 反演风场
    if quiver:
        aq = ax.quiver(vellon[cond], vellat[cond],
                    display.analysis_u[cond],
                    display.analysis_v[cond],
                    pivot='middle',
                    color='black',
                    scale=600,
                    #width=0.002,
                    #headwidth=4., headlength=8., headaxislength=4,
                    )
        ax.quiverkey(aq, 0.95, 0.95, 10, r'$10 m/s$')
    else:
        abar = ax.barbs(vellon[cond], vellat[cond],
                        display.analysis_u[cond],
                        display.analysis_v[cond],
                        length=5,
                        pivot='middle',
                        )

    # 添加底图
    add_province(shpfn, ax)
    #ax.set_xlim([117, 120.5])
    #ax.set_ylim([30.5, 33.5])

    auto_ticks(ax)

    if center == 'radar':
        xmin, ymin = get_lon_lat(radar, x=x[0]/1000, y=y[0]/1000)
        xmax, ymax = get_lon_lat(radar, x=x[1]/1000, y=y[1]/1000)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

    ax.set_title(title)
    if savefig:
        fig.savefig(savefig, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        plt.show()


if __name__=='__main__':
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    import pyart
    from pyart.io.sband_archive import read_sband_archive
    from pycwr.io.auto_io import radar_io
    import singledop


    # 雷达文件和底图
    filename = 'Z_RADR_I_Z9250_20190810015700_O_DOR_SA_CAP.bin'
    shpfn = 'data/shps/cnhimap.shp'

    # 读取雷达反射率
    sweep = 0    # 体扫仰角序号
    field = 'reflectivity'
    radar_info = pd.read_json('data/radar_info.json')

    if 1:
        radar = read_sband_archive(filename, station=get_radar_info(filename, radar_info))
    else:
        radar = radar_io(filename).ToPyartRadar()

    # 风场反演
    velfn = ('.').join(filename.split('.')[0:-1]) + '.nc'
    save = False   # 是否存储风场反演结果
    if os.path.exists(velfn):
        print(f'{velfn} 已存在...')
        display = singledop.NetcdfSave(velfn, radar=filename)
        vellon, vellat = get_lon_lat(radar, display=display)
        plot(radar, field, sweep, display, vellon, vellat, savefig='ref_vel.png')
    else:
        print('根据径向速度反演风场...')
        display, vellon, vellat = retrive_wind(radar, sweep_number=4, lonlat=True, filename=velfn)
        plot(radar, field, sweep, display, vellon, vellat)

