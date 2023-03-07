import os
import sys
import platform
import ctypes
import hashlib

__all__ = [
    'get_fft_lib_name',
    'get_fft_lib',
    'get_fft_lib_fp',
    'get_lib_md5',
    'set_fft_lib',
]
__LIBRARY = {}


def get_fft_lib_name(system, lib_ext=None):
    """
    Get fft lib name

    Parameters
    ----------
    system: str
        system platform

        - 'Darwin': mac
        - 'Linux': linux

    lib_ext: None or str
        - None: Default lib
        - str: libaudioflux_{str}.xxx

    Returns
    -------
    out: str
    """
    system = system.lower()
    flux_lib_name = 'libaudioflux'

    if system == 'darwin':
        suffix = 'dylib'
    elif system == 'linux':
        suffix = 'so'
    elif system.startswith('win'):
        suffix = 'dll'
    else:
        raise ValueError(f'Platform[{system}] is not support.')
    if lib_ext:
        flux_lib_name = f'{flux_lib_name}_{lib_ext}'
    flux_lib_name = f'{flux_lib_name}.{suffix}'
    return flux_lib_name


def get_fft_lib():
    """
    Get fft lib

    Returns
    -------
    out: fft lib file
    """
    return __LIBRARY['_lib']


def get_fft_lib_fp():
    """
    Get fft lib file path

    Returns
    -------
    out: str
    """
    return __LIBRARY['_lib_fp']


def get_lib_md5():
    """
    Get lib md5

    Returns
    -------
    out: str
    """
    with open(get_fft_lib_fp(), 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def set_fft_lib(system=None, *, lib_ext=None):
    """
    Set fft lib

    Parameters
    ----------
    system: None or str
        system platform

        - None: Automatically obtain the current system platform
        - 'Darwin': mac
        - 'Linux': linux

    lib_ext: None or str
        - None: Default lib
        - str: libaudioflux_{str}.xxx
    """
    if system is None:
        system = platform.system().lower()

    if system.startswith('win'):
        load_library = ctypes.windll.LoadLibrary
    else:
        load_library = ctypes.cdll.LoadLibrary
    lib_path = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'lib')

    if system.startswith('win') and (sys.version_info.major == 3 and sys.version_info.minor <= 7):
        fft_path = os.path.join(lib_path, 'libfftw3f-3.dll')
        if os.path.exists(fft_path):
            load_library(fft_path)
        else:
            print(f'not found {fft_path}')

    # load c lib
    lib_name = get_fft_lib_name(system, lib_ext)
    lib_fp = os.path.join(lib_path, lib_name)

    __LIBRARY['_lib'] = load_library(lib_fp)
    __LIBRARY['_lib_fp'] = lib_fp


set_fft_lib(None)
