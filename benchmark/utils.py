import numpy as np


def gen_data(time_step, radix2exp, slide_length, dtype=np.float32):
    fft_length = 1 << radix2exp
    return np.random.randn(fft_length + (time_step - 1) * slide_length).astype(dtype)
