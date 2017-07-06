#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% init
from itertools import compress
from functools import reduce
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from yaml import load as load_yaml, dump as dump_yaml
from h5py import File as H5File
from tqdm import tqdm
from anatools import Reader
from numba import jit


@jit
def sum_imgs(img1, img2):
    return np.sum((img1, img2), axis=0, dtype='double')


@jit
def reduce_imgs(imgs, where):
    compressed = compress(imgs, where)
    n = where.sum()
    return reduce(sum_imgs, tqdm(compressed, total=n)) / n


def load_formated_yaml(string):
    loaded = load_yaml(string)
    if 'format' not in loaded:
        return loaded
    fmt = loaded.pop('format')
    return load_yaml(dump_yaml(loaded, default_flow_style=False)
                     .format(**fmt))


def is_between(fr, to):
    def b(v):
        return (fr < v) & (v < to)
    return b


def digitizer(bins):
    def digitized(arr):
        tags = np.digitize(arr, bins)
        return bins[tags]
    return digitized

# %% read config file
parser = ArgumentParser(prog='reduce_imgs')
parser.add_argument('filename', nargs='?', default=None,
                    help="Config filename which format is 'yaml'")

args = parser.parse_args()
filename = args.filename
if filename is None:
    config_yaml = """
---
filenames:
    - /Volumes/store/20144091/ResScan5/Run_404/rawdata/*.h5
    - /Volumes/store/20144091/ResScan5/Run_405/rawdata/*.h5
prefix: /Volumes/store/20144091/ResScan5/combined/runs_404_405_
key_map:
    imgs: /vmi/andor  # intensity
    tags: /bunches  # 1
    bg_period: /Background_Period  # tags
    fel_wavelength: /photon_source/FEL01/wavelength  # nm
    fel_intensity: /photon_diagnostics/FEL01/I0_monitor/iom_sh_a  # uJ
    fel_hor_spectrum: /photon_diagnostics/Spectrometer/hor_spectrum
    fel_ver_spectrum: /photon_diagnostics/Spectrometer/vert_spectrum
    ion_flighttime: /digitizer/channel1  # ns
    ph_shift: /photon_source/FEL01/PhaseShifter5/DeltaPhase  # machine unit
tag_offset: 1
intensity_lim: [-.inf, .inf]  # '.inf' for infinity
phase_bins: [0.1, 2.0, 0.15]
phase_lim: null  # null to ignore
...
"""
else:    
    with open(args.filename) as f:
        config_yaml = f.read()

config = load_formated_yaml(config_yaml)
__filenames = config['filenames']
prefix = config['prefix']
key_map = config['key_map']
intlim = config['intensity_lim']
tag_offset = config['tag_offset']
phase_bins = config['phase_bins']
phase_lim = config['phase_lim']

# %% reduce
good_int = is_between(*intlim)
digitized = digitizer(np.arange(*phase_bins))

# file filter
with Reader(*__filenames, map=key_map) as reader:
    n = len(reader.files)
    fulllist = (f.filename for f in reader.files)
    phases = digitized(np.fromiter(reader['ph_shift.each_file'], 'double'))
    mask = (np.zeros(n, 'bool') |
            (phases!=phase_lim if phase_lim is not None else False))
    filenames = tuple(compress(fulllist, ~mask))

# event filter
with Reader(*filenames, map=key_map) as reader:
    imgs = lambda: reader['imgs']
    tags = np.fromiter(reader['tags'], 'int')-tag_offset
    n, *_ = tags.shape
    ints = np.fromiter(reader['fel_intensity'], 'double')
    period = next(reader['bg_period.each_file'])
    mask = ((np.fromiter((img.sum() for img in tqdm(imgs(), total=n)), 'double') == 0) |
            ~good_int(ints))
    isbg = (np.mod(tags, period) == 0) & ~mask
    issig = (np.mod(tags, period) != 0) & ~mask
    isdiff = issig
    bgimg = reduce_imgs(imgs(), isbg)
    sigimg = reduce_imgs(imgs(), issig)
    diffimg = sigimg - bgimg

# write
with open('{}config.yaml'.format(prefix), 'w') as f:
    f.write(dump_yaml(config, default_flow_style=False))

# with H5File('{}reduced.h5'.format(prefix), 'w') as f:  # h5py bug
try:
    f = H5File('{}reduced.h5'.format(prefix), 'w')
except OSError:
    f = H5File('{}reduced.h5'.format(prefix), 'w')
for name in ('bg', 'sig', 'diff'):
    img = eval('{}img'.format(name))
    n = eval('is{}.sum()'.format(name))
    f['{}img'.format(name)] = img
    f['n{}'.format(name)] = n
    plt.figure(figsize=(8, 8))
    plt.imshow(img)
    plt.title('reduced {} imgs (n={})'.format(name, n))
    plt.savefig('{}{}img.png'.format(prefix, name))
f.close()

plt.figure(figsize=(8, 8))
plt.hist(ints, bins=500)
plt.title('intensity dist')
plt.savefig('{}intensities.png'.format(prefix))
