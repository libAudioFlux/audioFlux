ChangeLog
=========

v0.1.5
------
* New features:
    * Optimize macOS/Linux performance.
* Modified API:
    * Remove the parameter of `audioflux.mel_spectrogram`: `is_reassign=False`. If you want to use reassign, you can use `audioflux.BFT`.
    * Remove the parameter of `audioflux.bark_spectrogram`: `is_reassign=False`. If you want to use reassign, you can use `audioflux.BFT`.
    * Remove the parameter of `audioflux.erb_spectrogram`: `is_reassign=False`. If you want to use reassign, you can use `audioflux.BFT`.
* Fix bug:
    * Support Ubuntu 18.04 system.


v0.1.4
------
* New features:
    * A variety of audio formats have been supported, such as MP3 and so on
    * The macOS ARM platforms have been supported.
    * Multi-channel audio has been supported.
* Modified API:
    * `audioflux.read` increases support for audio file lists and directory, and add the parameters of monocular conversion and modification sampling rate.
    * Modify the `audioflux.Temporal` API.
    * Modify `audioflux.CWT` parameter: `num=84`.
    * Modify `audioflux.PWT` parameter: `num=84, low_fre=None, high_fre=None`.
    * Modify `audioflux.WSST` parameter: `num=84`.
    * `audioflux.chirp` use `sicpy` to implement. And modify parameter `phi=None`, remove parameter `linear`, add parameter `method='logarithmic'`.
    * Add `audioflux.resample`/`audioflux.convert_mono` method.
* Removed API:
    * all debug/enable_debug methods.
    * audioflux.Resample.enable_continue
    * audioflux.WindowResample.enable_continue
    * audioflux.NSGT.get_cell_data
    * audioflux.PWT.pwt_det
    * audioflux.PWT.enable_det
    * audioflux.Temporal.ezr
    * audioflux.Temporal.temporal
* Fix bug:
    * Memory leakage in `audioflux.Temporal` object.
    * `audioflux.WindowResample.__init__` uses C method.
    * Fix the known bug.
