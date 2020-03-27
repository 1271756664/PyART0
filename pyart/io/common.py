"""
pyart.io.common
===============

Input/output routines common to many file formats.

.. autosummary::
    :toctree: generated/

    prepare_for_read
    stringarray_to_chararray
    _test_arguments
    make_time_unit_str

"""

import bz2
import gzip

import numpy as np
import netCDF4


def prepare_for_read(filename):
    """
    Return a file like object read for reading.

    Open a file for reading in binary mode with transparent decompression of
    Gzip and BZip2 files. The resulting file-like object should be closed.

    Parameters
    ----------
    filename : str or file-like object
        Filename or file-like object which will be opened. File-like objects
        will not be examined for compressed data.

    Returns
    -------
    file_like : file-like object
        File like object from which data can be read.

    """
    # if a file-like object was provided, return
    if hasattr(filename, 'read'):   # file-like object
        return filename

    # look for compressed data by examining the first few bytes
    fh = open(filename, 'rb')
    magic = fh.read(3)
    fh.close()

    if magic.startswith(b'\x1f\x8b'):
        return gzip.GzipFile(filename, 'rb')

    if magic.startswith(b'BZh'):
        return bz2.BZ2File(filename, 'rb')

    return open(filename, 'rb')


def stringarray_to_chararray(arr, numchars=None):
    """
    Convert an string array to a character array with one extra dimension.

    Parameters
    ----------
    arr : array
        Array with numpy dtype 'SN', where N is the number of characters
        in the string.

    numchars : int
        Number of characters used to represent the string. If numchar > N
        the results will be padded on the right with blanks. The default,
        None will use N.

    Returns
    -------
    chararr : array
        Array with dtype 'S1' and shape = arr.shape + (numchars, ).

    """
    carr = netCDF4.stringtochar(arr)
    if numchars is None:
        return carr

    arr_numchars = carr.shape[-1]
    if numchars <= arr_numchars:
        raise ValueError('numchars must be >= %i' % (arr_numchars))
    chararr = np.zeros(arr.shape + (numchars, ), dtype='S1')
    chararr[..., :arr_numchars] = carr[:]
    return chararr


def _test_arguments(dic):
    """ Issue a warning if receive non-empty argument dict. """
    if dic:
        import warnings
        warnings.warn('Unexpected arguments: %s' % dic.keys())


def make_time_unit_str(dtobj):
    """ Return a time unit string from a datetime object. """
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")


def radar_format(filename):
    """Identify radar file format

    Parameters
    ----------
    filename : string
             radar file name

    Returns
    -------
    string
    """
    if hasattr(filename, 'read'):
        return filename
    fh = prepare_for_read(filename)
    flag = fh.read(28)
    size = len(fh.read()) + 28
    fh.seek(100, 0)
    sc_flag = fh.read(9)
    fh.seek(116, 0)
    cc_flag = fh.read(9)
    fh.close()
    if flag[:4] == b'RSTM':
        return "WSR98D"
    elif flag[14:16] == b'\x01\x00':
        return "SAB"
    elif (size-1024)%3000 == 0 and cc_flag == b"CINRAD/CC":
        return "CC"
    elif (size-1024)%4000 == 0 and (sc_flag == b"CINRAD/SC" or sc_flag == b"CINRAD/CD"):
        return "SC"
    else:
        return _get_radar_type(filename)


def _get_radar_type(filename):
    """Identify type of radar based on radar file format
    :param filename:
    """
    name = os.path.basename(filename)
    station_id = [int(name[idx:idx + 4]) for idx in range(len(name) - 4) if name[idx:idx + 4].isdigit()][0]
    if station_id not in radar_info.index:
        return None
    Datatype = radar_info.loc[station_id, "Datatype"]
    if Datatype in ["SA", "SB", "CB", "SC", "CD"]:
        return "SAB"
    elif Datatype in ["CC", "CCJ"]:
        return "CC"
    else:
        return None
