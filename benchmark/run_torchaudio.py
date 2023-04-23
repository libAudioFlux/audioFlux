import os
import sys
import time
import argparse

sys.path.append(__file__.rsplit(os.path.sep, 1)[0])  # NOQA
import torch
import torchaudio

from utils import gen_data


def mel(time_step, runtimes, radix2_exp=11, slide_length=512):
    n_fft = 1 << radix2_exp

    transform = torchaudio.transforms.MelSpectrogram(32000, n_fft=n_fft, hop_length=slide_length, n_mels=128,
                                                     center=False)
    data_arr = gen_data(time_step, radix2_exp, slide_length)
    data_arr = torch.from_numpy(data_arr)
    r = transform(data_arr)
    # print(r.shape)

    total_time = 0
    for _ in range(runtimes):
        data_arr = gen_data(time_step, radix2_exp, slide_length)
        data_arr = torch.from_numpy(data_arr)

        s = time.time()
        r = transform(data_arr)
        e = time.time()
        total_time += e - s

    return total_time * 1000


FEATURES_DIC = {'mel': mel}


def main(runtimes, time_steps, feature_name, radix2_exp, slide_length):
    print(f'project:{torchaudio.__name__}-{torchaudio.__version__}/{torch.__name__}-{torch.__version__}')
    print(f'runtimes:{runtimes}')
    print(f'time_steps:{time_steps}')
    print(f'feature_name:{feature_name}')
    print(f'radix2_exp:{radix2_exp}')
    print(f'slide_length:{slide_length}')
    print('-' * 10)
    fn = FEATURES_DIC[feature_name]

    for time_step in map(int, time_steps.split(',')):
        use_time = fn(time_step, runtimes, radix2_exp, slide_length)
        print(torchaudio.__name__, fn.__name__, time_step, f'{use_time:>.8f}')


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
