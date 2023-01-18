Feature Extractor
------------------

You can use `FeatureExtractor <feature/featureExtractor:FeatureExtractor>` for batch feature extraction

.. code-block:: python
    :linenos:

    # Get a 880Hz's audio file
    import audioflux as af
    from audioflux.type import SpectralFilterBankScaleType
    sample_path = af.utils.sample_path('880')
    audio_arr, sr = af.read(sample_path)

    # Create FeatureExtractor object and extract spectrogram
    fa_obj = af.FeatureExtractor(transforms=['bft', 'cwt', 'cqt'], samplate=sr, radix2_exp=12,
                                 scale_type=SpectralFilterBankScaleType.OCTAVE)
    spec_result = fa_obj.spectrogram(audio_arr, is_continue=True)

    #Extract spectral/xxcc/deconv
    spectral_result = fa_obj.spectral(spec_result, spectral='flux',
                                      spectral_kw={'is_positive': True})
    xxcc_result = fa_obj.xxcc(spec_result, cc_num=13)
    deconv_result = fa_obj.deconv(spec_result)
