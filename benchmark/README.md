## Benchmark

The time required to compute the mel-spectrogram for 1000 samples according to the TimeStep.

#### Run Benchmark

* If you want to compare and test multiple libraries, you can use:

    ```shell
      $ python run_benchmark.py -p audioflux,torchaudio,librosa -r 50 -er 50  -t 1,5,10,100,500,1000,2000,3000
    ```
    * -p: The library name
    * -r: The number of times to loop `run_xxx.py `
    * -er: The number of times each `run_xxx.py` loops mel
    * -t: TimeStep

* If you want to test a single library, you can use:

```shell
  $ python run_audioflux.py -r 50 -t 1,5,10,100,500,1000,2000,3000
```

* If you want to see more usage instructions, you can execute `python run_xxx.py --help`

