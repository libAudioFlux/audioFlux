from audioflux.fftlib import get_fft_lib


class Base(object):
    def __init__(self, obj):
        self._obj = obj
        self._lib = get_fft_lib()
        self._is_created = False
