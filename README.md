# MFFC_DTW
#### Audio comparison algorithm steps:
| Nr. | Step                              | File                        | Description                               |
| --- | --------------------------------- | --------------------------- | ----------------------------------------- |
| 1.  | Pre-emphasis                      | `steps/preemphasis.cpp`     | Applies pre-emphasis on the sample        |
| 2.  | Framing                           | `steps/framing.cpp`         | Splits sample into frames                 |
| 3.  | Windowing                         | `steps/windowing.cpp`       | Applies Hamming window to each frame      |
| 4.  | Discrete Fourier transform        | `steps/dft.cpp`             | Calculates DFT for each frame             |
| 5.  | Power spectrum                    | `steps/powerspectrum.cpp`   | Calcaulates power spectrum                |
| 6.  | Mel's filter Bank                 | `steps/filterbanks.cpp`     | Calculates Mel's filter bank              |
| 7.  | Logarithm                         | `steps/logarithm.cpp`       | Applies logarithm to the resulting matrix |
| 8.  | Discrete Cosine Transform (Type 2)| `steps/dct.cpp`             | Calculates DCT for each row of the matrix |
| 9.  | Dynamic Time Warping              | `steps/dtw.cpp`             | Calculates "distance" between samples     |
| 10. | Normalization                     | *included in previous step* | Applies normalization                     |

#### To compile (in terminal):
1. `cmake .`
2. `make`

#### To run (each individual step):
1. `./steps/bin/preemph IN_SAMPLE_PATH OUT_SAMPLE_PATH` (f32 raw sample with 44100 sample rate)
2. `./steps/bin/framing IN_SAMPLE_PATH OUT_FOLDER PATH`
3. `./steps/bin/windowing IN_FOLDER_PATH OUT_FOLDER PATH`
4. `./steps/bin/dft IN_FOLDER_PATH OUT_FOLDER PATH`
5. `./steps/bin/powspec IN_FOLDER_PATH OUT_FOLDER PATH`
6. `./steps/bin/filterbanks IN_FOLDER_PATH OUT_FOLDER PATH`
7. `./steps/bin/log IN_FOLDER_PATH OUT_FOLDER PATH`
8. `./steps/bin/dct IN_FOLDER_PATH OUT_FOLDER PATH`
9. `./steps/bin/dtw IN_FOLDER_PATH IN_FOLDER PATH` (two folders with DCT calculated, outputs normalized similiarity)

#### To run (full comparison):
`./full/bin/compare IN_SAMPLE_PATH IN_SAMPLE_PATH` (f32 raw sample with 44100 sample rate, outputs normalized similarity)
