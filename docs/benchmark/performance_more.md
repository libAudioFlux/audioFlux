## More detailed performance testing

#### Linux - AMD

	- OS: Ubuntu 20.04.4 LTS
	- CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
  
<img src='../image/benchmark/linux_amd_1.png'>

| TimeStep | audioflux   | torchaudio  | librosa     | 
| -------- |  --------   |  ---------  |  ------     |
| 1        | 0.00004294s | 0.00007707s | 0.00241958s |
| 5        | 0.00014878s | 0.00105589s | 0.00352610s |
| 10       | 0.00018374s | 0.00083975s | 0.00346499s |
| 100      | 0.00067030s | 0.00061876s | 0.00663217s |
| 500      | 0.00094893s | 0.00129189s | 0.01645968s |
| 1000     | 0.00143854s | 0.00223126s | 0.02778358s |
| 2000     | 0.00308714s | 0.00410869s | 0.04512714s |
| 3000     | 0.00490343s | 0.00586299s | 0.05162876s |

#### Linux - Intel

	- OS: Ubuntu 20.04.4 LTS
	- CPU: Intel(R) Core(TM) i7-6850K CPU @ 3.60GHz
      
<img src='../image/benchmark/linux_intel_1.png'>
   
| TimeStep | audioflux   | torchaudio  | librosa     | 
| -------- |  --------   |  ---------  |  ------     |
| 1        | 0.00008106s | 0.00011043s | 0.00551295s |
| 5        | 0.00011654s | 0.00016005s | 0.00577631s |
| 10       | 0.00029173s | 0.00015352s | 0.00613656s |
| 100      | 0.00118150s | 0.00039958s | 0.01061641s |
| 500      | 0.00223883s | 0.00158323s | 0.02899823s |
| 1000     | 0.00442723s | 0.00398896s | 0.05197518s |
| 2000     | 0.00873121s | 0.00828444s | 0.06113923s |
| 3000     | 0.01307378s | 0.01214323s | 0.07006395s |


#### macOS - Intel

	- OS: 12.6.1 (21G217)
	- CPU: 3.8GHz 8‑core 10th-generation Intel Core i7, Turbo Boost up to 5.0GHz

<img src='../image/benchmark/mac_x86_1.png'>
        
| TimeStep | audioflux   | torchaudio  | librosa     | 
| -------- |  --------   |  ---------  |  ------     |
| 1        | 0.00007605s | 0.00006451s | 0.00170139s |
| 5        | 0.00014946s | 0.00008464s | 0.00186964s |
| 10       | 0.00016641s | 0.00010762s | 0.00200865s |
| 100      | 0.00046902s | 0.00083551s | 0.00328890s |
| 500      | 0.00108860s | 0.00505824s | 0.00898265s |
| 1000     | 0.00264029s | 0.00978269s | 0.01824391s |
| 2000     | 0.00540025s | 0.01508991s | 0.03368184s |
| 3000     | 0.00792596s | 0.02484823s | 0.04735941s |

#### macOS - M1
	
	- OS: 12.4 (21F79)
	- CPU: Apple M1
    
<img src='../image/benchmark/mac_arm_1.png'>
        
| TimeStep | audioflux   | torchaudio  | librosa     | 
| -------- |  --------   |  ---------  |  ------     |
| 1        | 0.00006110s | 0.00006874s | 0.00222518s |
| 5        | 0.00023444s | 0.00007922s | 0.00255907s |
| 10       | 0.00020691s | 0.00011090s | 0.00271813s |
| 100      | 0.00068694s | 0.00063625s | 0.00474433s |
| 500      | 0.00147420s | 0.00337597s | 0.01383887s |
| 1000     | 0.00300926s | 0.00676275s | 0.02524646s |
| 2000     | 0.00599781s | 0.01269573s | 0.04784029s |
| 3000     | 0.00876306s | 0.01903391s | 0.06940428s |


Based on the TimeStep of 1/5/10/100/500/1000/2000/3000 to calculate the time required for a Mel spectrogram.
Where fft_len=2048, slide_len=512, sampling_rate=32000.

#### Benchmark Script
[https://github.com/libAudioFlux/audioFlux/tree/master/benchmark](https://github.com/libAudioFlux/audioFlux/tree/master/benchmark)

------
### Other Test

#### Server Performance

Each sample data is 128ms(sampling rate: 32000, data length: 4096).

The total time spent on extracting features for 1000 sample data.

| Package    | [audioFlux](https://github.com/libAudioFlux/audioFlux) | [librosa](https://github.com/librosa/librosa) | [pyAudioAnalysis](https://github.com/tyiannak/pyAudioAnalysis) | [python\_speech\_features](https://github.com/jameslyons/python_speech_features) |
| ------ |  ------ |  ------ |  ------ |  ------ | 
| Mel    | 0.777s    | 2.967s  | --              | --                       |
| MFCC   | 0.797s    | 2.963s  | 0.805s          | 2.150s                   |
| CQT    | 5.743s    | 21.477s | --              | --                       |
| Chroma | 0.155s    | 2.174s  | 1.287s          | --                       |

#### Mobile Performance

For 128ms audio data per frame(sampling rate: 32000, data length: 4096).

The time spent on extracting features for 1 frame data.

| Mobile | iPhone 13 Pro | iPhone X | Honor V40 | OPPO Reno4 SE 5G |
| ------ |  ------ |  ------ |  ------ |  ------ | 
| Mel    | 0.249ms       | 0.359ms  | 0.313ms   | 0.891ms          |
| MFCC   | 0.249ms       | 0.361ms  | 0.315ms   | 1.116ms          |
| CQT    | 0.350ms       | 0.609ms  | 0.786ms   | 1.779ms          |
| Chroma | 0.354ms       | 0.615ms  | 0.803ms   | 1.775ms          |