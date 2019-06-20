# coding : utf-8
"""
pyart.io.c98d_archive

functions for read C 98D Doppolar Radar archive file
"""

import warnings

import numpy as np

from .C98DRadFile import C98DRadFile
from pyart.config import FileMetadata, get_fillvalue
from pyart.core.radar import Radar
from pyart.io.common import make_time_unit_str, _test_arguments, prepare_for_read


def c98dfile_archive(filename, field_names=None, additional_metadata=None,
                        file_field_names=False, exclude_fields=None, cutnum=None,
                        delay_field_loading=False, **kwargs):

    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('c98d_archive', field_names,
                                additional_metadata, file_field_names,
                                exclude_fields)

    nfile = C98DRadFile(prepare_for_read(filename))
#    scan_info = nfile.scan_info
    
    if cutnum is None:
        cutnum = nfile.scans
    else:
        cutnum = list(cutnum)
    
    # time
    time = filemetadata('time')
    _time = nfile.radial_info['seconds']
    time['data'] = _time
    time['units'] = nfile.get_volume_start_time 

    # range
    _range = filemetadata('range')
    _range['data'] = nfile.get_range
#    _range['meters_to_center_of_first_gate'] = float(first_gate)
#    _range['meters_between_gates'] = float(gate_spacing)
    
    # fields
    fields = {}
    field_names = nfile.get_moment_type
    for field_name in field_names:
        dic = filemetadata(field_name)
        for i, cn in enumerate(cutnum):
            dic['_FillValue'] = get_fillvalue()
            if i == 0:
                dic['data'] = nfile.get_data(field_name, cn)
            else:
                fndata = nfile.get_data(field_name, cn)
                
                if fndata.shape[1] != dic['data'].shape[1]:
                    if dic['data'].shape[1] - fndata.shape[1] > 0:
                        apps = np.ones((fndata.shape[0], dic['data'].shape[1] - fndata.shape[1])) * np.nan
                    else:
                        raise ValueError('something wrong!')
                        
                    fndata = np.c_[fndata, apps]
                dic['data'] = np.append(dic['data'], fndata, axis=0)
                
        fields.update({field_name: dic})  
    
    # scan_type
    scan_type = nfile.scan_type

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    
    lat, lon, height = nfile.get_location
    
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([height], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')    
        
    sweep_number['data'] = np.arange(nfile.cutnum, dtype='int32')
    sweep_mode['data'] = np.array(
        nfile.cutnum * ['azimuth_surveillance'], dtype='S')

    sweep_end_ray_index['data'] = np.array(nfile.cut_end)
    sweep_start_ray_index['data'] = np.array(nfile.cut_start)

    # azimuth, elevation, fixed_angle
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    fixed_angle = filemetadata('fixed_angle')
    azimuth['data'] = np.array(nfile.get_azimuth)
    elevation['data'] = np.array(nfile.get_elevation)
    fixed_angle['data'] = nfile.get_target_angles   
    
    # instrument_parameters
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')
    nyquist_velocity['data'] = None
    unambiguous_range['data'] = None

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity}
    

    return Radar(time, _range, fields, filemetadata, scan_type, 
                 latitude, longitude, altitude, 
                 sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
                 sweep_end_ray_index,
                 azimuth, elevation,
                 instrument_parameters=None), nfile
        