Transform
=========

Main
--------------

In the time–frequency representation, main transform algorithm:

.. toctree::
    :maxdepth: 1

    bft
    nsgt
    cwt
    pwt

+-----------+--------+----------+-----+------+-----+--------+-----+
| Transform | Linear | Linspace | Mel | Bark | Erb | Octave | Log |
+===========+========+==========+=====+======+=====+========+=====+
| BFT       | ✅     | ✅       | ✅  | ✅   | ✅  | ✅     | ✅  |
+-----------+--------+----------+-----+------+-----+--------+-----+
| NSGT      | ✅     | ✅       | ✅  | ✅   | ✅  | ✅     | ✅  |
+-----------+--------+----------+-----+------+-----+--------+-----+
| CWT       | ✅     | ✅       | ✅  | ✅   | ✅  | ✅     | ✅  |
+-----------+--------+----------+-----+------+-----+--------+-----+
| PWT       | ✅     | ✅       | ✅  | ✅   | ✅  | ✅     | ✅  |
+-----------+--------+----------+-----+------+-----+--------+-----+

The above transform supports all the following frequency scale types:

* Linear - Short-time Fourier transform spectrogram.
* Linspace - Linspace-scale spectrogram.
* Mel - Mel-scale spectrogram.
* Bark - Bark-scale spectrogram.
* Erb - Erb-scale spectrogram.
* Octave - Octave-scale spectrogram.
* Log - Logarithmic-scale spectrogram.


Other
---------------

The following transform are not supports multiple frequency scale types, only used as independent transform:

.. toctree::
    :maxdepth: 1

    cqt
    st
    fst
    dwt
    wpt
    swt


Synchronized Squeezing
----------------------

The synchrosqueezing or reassignment is a technique for sharpening a time-frequency representation,
contains the following algorithms:

.. toctree::
    :maxdepth: 1

    reassign
    synsq
    wsst

+----------+-----------+-------+--------------+
| method   | transform | order | independence |
+==========+===========+=======+==============+
| Reassign | BFT/STFT  | ☑️    | ❌           |
+----------+-----------+-------+--------------+
| Synsq    | CWT       | ✅    | ✅           |
+----------+-----------+-------+--------------+
| WSST     | CWT       | ✅    | ❌           |
+----------+-----------+-------+--------------+