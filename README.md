# MFFC_DTW
#### Audio comparison algorithm steps:
| Nr. | Step                              | File                        | Description                               |
| --- | --------------------------------- | --------------------------- | ----------------------------------------- |
| 1.  | Pre-emphasis                      | `steps/preemphasis.cpp`     | Applies pre-emphasis on the sample        |
| 2.  | Framing                           | `steps/framing.cpp`         | Splits sample into frames                 |
| 3.  | Windowing                         | `steps/windowing.cpp`       | Applies Hamming window to each frame      |
| 4.  | Discrete Fourier transform        | `steps/dft.cpp`             | Calculates DFT for each frame             |
| 5.  | Mel's filter Bank                 | `steps/filterbanks.cpp`     | Calculates Mel's filter bank              |
| 6.  | Power spectrum                    | `steps/powerspectrum.cpp`   | Calcaulates power spectrum                |
| 7.  | Logarithm                         | `steps/logarithm.cpp`       | Applies logarithm to the resulting matrix |
| 8.  | Discrete Cosine Transform (Type 2)| `steps/dct.cpp`             | Calculates DCT for each row of the matrix |
| 9.  | Dynamic Time Warping              | `steps/dtw.cpp`             | Calculates "distance" between samples     |
| 10. | Normalization                     | *included in previous step* | Applies normalization                     |

#### To run:
