Feature
=======

The feature module contains the following algorithms:

.. toctree::
    :maxdepth: 1

    spectral
    xxcc
    deconv
    cepstrogram
    temporal
    featureExtractor

+----------+-----+------+-----+-----+-----+-----+----+-----+-----+-----+-----+
| feature  | BFT | NSGT | CWT | PWT | CQT | VQT | ST | FST | DWT | WPT | SWT |
+==========+=====+======+=====+=====+=====+=====+====+=====+=====+=====+=====+
| spectral | ✅  | ✅   | ✅  | ✅  | ✅  | ✅  | ✅ | ✅  | ✅  | ✅  | ❌  |
+----------+-----+------+-----+-----+-----+-----+----+-----+-----+-----+-----+
| xxcc     | ✅  | ✅   | ✅  | ✅  | ✅  | ✅  | ✅ | ✅  | ✅  | ✅  | ❌  |
+----------+-----+------+-----+-----+-----+-----+----+-----+-----+-----+-----+
| deconv   | ✅  | ✅   | ✅  | ✅  | ✅  | ✅  | ✅ | ✅  | ✅  | ✅  | ❌  |
+----------+-----+------+-----+-----+-----+-----+----+-----+-----+-----+-----+
| chroma   | ☑️  | ❌   | ❌  | ❌  | ✅  | ❌  | ❌ | ❌  | ❌  | ❌  | ❌  |
+----------+-----+------+-----+-----+-----+-----+----+-----+-----+-----+-----+
