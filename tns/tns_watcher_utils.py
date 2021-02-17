"""
These functions were taken from Dima Duev's script:
https://github.com/dmitryduev/kowalski/blob/master/kowalski/utils.py
"""

import os
import datetime

pi = 3.141592653589793


def radec_str2geojson(ra_str, dec_str):

    # hms -> ::, dms -> ::
    if isinstance(ra_str, str) and isinstance(dec_str, str):
        if ("h" in ra_str) and ("m" in ra_str) and ("s" in ra_str):
            ra_str = ra_str[:-1]  # strip 's' at the end
            for char in ("h", "m"):
                ra_str = ra_str.replace(char, ":")
        if ("d" in dec_str) and ("m" in dec_str) and ("s" in dec_str):
            dec_str = dec_str[:-1]  # strip 's' at the end
            for char in ("d", "m"):
                dec_str = dec_str.replace(char, ":")

        if (":" in ra_str) and (":" in dec_str):
            ra, dec = radec_str2rad(ra_str, dec_str)
            # convert to geojson-friendly degrees:
            ra = ra * 180.0 / pi - 180.0
            dec = dec * 180.0 / pi
        else:
            raise Exception("unrecognized string ra/dec format.")
    else:
        # already in degrees?
        ra = float(ra_str)
        # geojson-friendly ra:
        ra -= 180.0
        dec = float(dec_str)

    return ra, dec


def radec_str2rad(_ra_str, _dec_str):
    """
    :param _ra_str: 'H:M:S'
    :param _dec_str: 'D:M:S'
    :return: ra, dec in rad
    """
    # convert to rad:
    _ra = list(map(float, _ra_str.split(":")))
    _ra = (_ra[0] + _ra[1] / 60.0 + _ra[2] / 3600.0) * pi / 12.0
    _dec = list(map(float, _dec_str.split(":")))
    _sign = -1 if _dec_str.strip()[0] == "-" else 1
    _dec = (
        _sign
        * (abs(_dec[0]) + abs(_dec[1]) / 60.0 + abs(_dec[2]) / 3600.0)
        * pi
        / 180.0
    )

    return _ra, _dec
