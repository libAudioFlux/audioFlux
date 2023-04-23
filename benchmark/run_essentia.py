import os
import sys
import time
import argparse

import numpy as np

sys.path.append(__file__.rsplit(os.path.sep, 1)[0])  # NOQA
import essentia.standard as es
import essentia

from utils import gen_data


def cal_frame_mel(audio, windowing, spectrum, melbands, pool):
    for frame in audio:
        frame_spectrum = spectrum(windowing(frame))
        frame_mel = melbands(frame_spectrum)
        pool.add('mel', frame_mel)
    return np.array(pool['mel'])


def mel(time_step, runtimes, radix2_exp=11, slide_length=512):
    n_fft = 1 << radix2_exp

    windowing = es.Windowing(type='hann', zeroPadding=0)
    spectrum = es.Spectrum(size=n_fft)
    melbands = es.MelBands(inputSize=n_fft // 2 + 1, sampleRate=32000, numberBands=128, type='power',
                           lowFrequencyBound=0, highFrequencyBound=16000, warpingFormula='slaneyMel',
                           weighting='warping', normalize='unit_sum')

    data_arr = gen_data(time_step, radix2_exp, slide_length)
    audio = es.FrameGenerator(data_arr, frameSize=n_fft, hopSize=slide_length, startFromZero=True,
                              lastFrameToEndOfFile=False)
    pool = essentia.Pool()
    r = cal_frame_mel(audio, windowing=windowing, spectrum=spectrum, melbands=melbands, pool=pool)
    # print(r.shape)

    total_time = 0
    for i in range(runtimes):
        data_arr = gen_data(time_step, radix2_exp, slide_length)
        audio = es.FrameGenerator(data_arr, frameSize=n_fft, hopSize=slide_length, startFromZero=True,
                                  lastFrameToEndOfFile=False)
        pool = essentia.Pool()
        s = time.time()
        r = cal_frame_mel(audio, windowing=windowing, spectrum=spectrum, melbands=melbands, pool=pool)
        e = time.time()
        total_time += e - s

    return total_time * 1000


FEATURES_DIC = {'mel': mel}


def main(runtimes, time_steps, feature_name, radix2_exp, slide_length):
    print(f'project:{essentia.__name__}-{essentia.__version__}')
    print(f'runtimes:{runtimes}')
    print(f'time_steps:{time_steps}')
    print(f'feature_name:{feature_name}')
    print(f'radix2_exp:{radix2_exp}')
    print(f'slide_length:{slide_length}')
    print('-' * 10)
    fn = FEATURES_DIC[feature_name]

    for time_step in map(int, time_steps.split(',')):
        use_time = fn(time_step, runtimes, radix2_exp, slide_length)
        print(essentia.__name__, fn.__name__, time_step, f'{use_time:>.8f}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--runtimes', default=10, type=int, help='run times')
    parser.add_argument('-t', '--time_steps', default='1', type=str, help='time_steps, like: 1,5,100')
    parser.add_argument('-f', '--feature', default='mel', choices=list(FEATURES_DIC), type=str, help='test feature')
    parser.add_argument('-e', '--radix2_exp', default=11, type=int, help='radix2_exp')
    parser.add_argument('-s', '--slide_length', default=512, type=int, help='slide_length')
    args = parser.parse_args()

    main(runtimes=args.runtimes, time_steps=args.time_steps, feature_name=args.feature,
         radix2_exp=args.radix2_exp, slide_length=args.slide_length)
