import os
import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

PROJECTS = [
    'audioflux',
    'torchaudio',
    'librosa',
    'essentia',
]

FEATURES = ['mel']


def show(data, save_path=None):
    show_data = {}
    xticks = []
    for key in sorted(data):
        _, project = key
        use_times = []
        for val in data[key]:
            xtick, use_time = val
            if xtick not in xticks:
                xticks.append(xtick)
            use_times.append(use_time)
        show_data[key] = np.asarray(use_times)

    bar_width = 0.25
    for i, key in enumerate(sorted(show_data)):
        _, project = key
        bar_pos = np.arange(len(xticks)) + (i * bar_width)
        plt.bar(bar_pos, show_data[key], width=bar_width, label=project)

    plt.xlabel('TimeStep')
    plt.ylabel('UseTime(s)')

    plt.xticks(np.arange(len(xticks)), xticks)

    plt.legend()
    if save_path:
        plt.savefig(save_path)
    plt.show()


def main(projects, runtimes, each_runtimes, time_steps, feature_name, radix2_exp, slide_length,
         show_img, img_path):
    bm_data = defaultdict(list)
    for project_idx, project in enumerate(projects.split(',')):
        project = project.strip()
        if project not in PROJECTS:
            raise ValueError(f'not found project: {project}')
        script_path = os.path.join(os.path.abspath(__file__).rsplit(os.path.sep, 1)[0], f'run_{project}.py')
        script_params = [
            f'--runtimes={each_runtimes}',
            f'--time_steps={time_steps}',
            f'--feature={feature_name}',
            f'--radix2_exp={radix2_exp}',
            f'--slide_length={slide_length}',
        ]
        cmd = f'{sys.executable} {script_path} {" ".join(script_params)}'

        total_time = defaultdict(int)
        project_info = {}
        project_info_key = []
        for _ in range(runtimes):
            start_flag = False
            f = os.popen(cmd)
            msg = f.read()
            msg = msg.strip()

            for i, row in enumerate(msg.split('\n')):
                row = row.strip()
                if not start_flag and row.startswith('---'):
                    start_flag = True
                    continue
                # if row.endswith(')'):
                #     print(row)
                #     continue
                if start_flag:
                    *_, s, t = row.split(' ')
                    use_time = float(t)
                    total_time[(i, s)] += use_time
                else:
                    k, v = row.split(':', 1)
                    if k not in project_info:
                        project_info[k] = v
                        project_info_key.append(k)

        for key in project_info_key:
            print(f'{key}:{project_info[key]}')
        print('-' * 10)
        for key in sorted(total_time):
            _, ts = key
            one_time = total_time[key] / runtimes / each_runtimes
            print(project, feature_name, ts, f'{one_time:>.8f}')
            bm_data[(project_idx, project_info['project'])].append((ts, one_time))
        print('=' * 10)
    if show_img:
        show(bm_data, img_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--projects', required=True, type=str, help='project {%s}' % PROJECTS)
    parser.add_argument('-r', '--runtimes', default=50, type=int, help='run times')
    parser.add_argument('-er', '--each_runtimes', default=10, type=int, help='each script run times')
    parser.add_argument('-t', '--time_steps', default='1', type=str, help='time_steps, like: 1,5,100')
    parser.add_argument('-f', '--feature', default='mel', choices=FEATURES, type=str, help='test feature')
    parser.add_argument('-e', '--radix2_exp', default=11, type=int, help='radix2_exp')
    parser.add_argument('-s', '--slide_length', default=512, type=int, help='slide_length')
    parser.add_argument('--show_img', default=False, type=bool, help='show matplotlib')
    parser.add_argument('--img_path', default=None, type=str, help='image path, like: ./a.png')
    args = parser.parse_args()

    main(projects=args.projects, runtimes=args.runtimes, each_runtimes=args.each_runtimes,
         time_steps=args.time_steps, feature_name=args.feature,
         radix2_exp=args.radix2_exp, slide_length=args.slide_length,
         show_img=args.show_img, img_path=args.img_path)
