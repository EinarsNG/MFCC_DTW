// Some parts of the code were ported from Python
// implementation at https://github.com/jameslyons/python_speech_features
#include <algorithm>
#include <mutex>
#include <fftw3.h>

#include <includes/steps.h>

// applies pre-emphasis to the sample
Vector<float> preemph(Vector<float> &data, float coeff)
{
  std::vector<float> res(data.size());
  for (size_t i = data.size() - 1; i >= 1; i--)
    res[i] = 1.0f * data[i] - data[i-1] * coeff;
  res[0] = 1.0f * data[0];
  return res;
}

// splits sample into frames
Vector2d<float> framing(
    Vector<float> &data,
    size_t sample_rate,
    float frame_size,
    float frame_step)
{
  size_t frame_size_samples = sample_rate * frame_size;
  size_t frame_step_samples = sample_rate * frame_step;

  Vector<float> padded = data;
  // add zero-padding in case it does not line up
  size_t numframes = 1;
  if (data.size() > frame_size_samples)
    numframes += ceil((static_cast<float>(data.size()) - frame_size_samples) / frame_step_samples);
  size_t pad_amount = ((numframes - 1) * frame_step_samples + frame_size_samples) - data.size();
  for (size_t i = 0; i < pad_amount; i++) padded.push_back(0);

  Vector2d<float> res;
  Vector<float> temp(frame_size_samples);
  for (size_t i = 0; i < padded.size(); i += frame_step_samples)
  {
    std::copy(padded.begin() + i,
        padded.begin() + i + frame_size_samples,
        temp.begin());
    res.push_back(temp);
  }
  return res;
}

// applies hamming window to a sample or frame
Vector<float> hamming_window(Vector<float>& data)
{
  Vector<float> res(data.size());
  for (size_t i = 0; i < data.size(); i++)
  {
    res[i] = data[i] * (0.54f - 0.46f * cos(2.0f * M_PI * (static_cast<float>(i)/(static_cast<float>(data.size()-1.0f)))));
  }
  return res;
}

// applies Discrete Fourier Transform to each frame
Vector2d<float> dft(
    Vector2d<float> frames,
    size_t nfft
)
{
  size_t new_col_size = nfft / 2 + 1;
  Vector2d<float> res;
  Vector<double> in_data(nfft);
  Vector<std::complex<double>> out_data(new_col_size);

  fftw_plan p = fftw_plan_dft_r2c_1d(nfft,
        in_data.data(),
        (fftw_complex *)out_data.data(),
        FFTW_ESTIMATE);
  if (p == nullptr)
  {
    printf("Error: There was a problem initiating FFTW_PLAN\n");
    return res;
  }

  size_t i = 0;
  for (auto & frame : frames)
  {
    Vector<float> temp(new_col_size);
    for (size_t i = 0; i < frame.size() && i < nfft; i++)
      in_data[i] = static_cast<double>(frame[i]);
    fftw_execute(p);
    for (size_t i = 0; i < new_col_size; i++)
      temp[i] = static_cast<float>(abs(out_data[i]));
    res.push_back(temp);
  }
  fftw_destroy_plan(p);

  return res;
}

// calculates power spectrum matrix
Vector2d<float> powspec(Vector2d<float> frames, size_t nfft)
{
  Vector2d<float> res;
  for (auto & frameset : frames)
  {
    Vector<float> temp(frameset.size());
    for (size_t i = 0; i < frameset.size(); i++)
      temp[i] = 1.0f / static_cast<float>(nfft) * pow(frameset[i], 2.0f);
    res.push_back(temp);
  }
  return res;
}

// calculates Mel's filter banks
Vector2d<float> filterbanks(size_t numfilt, size_t nfft, size_t sample_rate)
{
  size_t lowfreq = 0;
  size_t highfreq = sample_rate / 2;

  float lowmel = hz_to_mel(lowfreq);
  float highmel = hz_to_mel(highfreq);

  Vector<float> melpoints = linearspace(lowmel, highmel, numfilt + 2);
  Vector<float> fft_bins;
  for (auto & point : melpoints)
    fft_bins.push_back(mel_to_fft_bin(point, nfft, sample_rate));

  Vector2d<float> res(numfilt, Vector<float>(nfft/2+1, 0));
  for (size_t j = 0; j < numfilt; j++)
  {
    float one = fft_bins[j];
    float two = fft_bins[j+1];
    float three = fft_bins[j+2];

    size_t start = one;
    size_t end = two;
    for (size_t i = start; i < end; i++)
      res[j][i] = (static_cast<float>(i) - one) / (two - one);

    start = two;
    end = three;
    for (size_t i = start; i < end; i++)
      res[j][i] = (three - static_cast<float>(i)) / (three - two);
  }
  return res;
}

// applies natural logarithm to each element in the matrix
Vector2d<float> log(Vector2d<float>& mtx)
{
  Vector2d<float> res = mtx;
  for (auto & row : res)
    for (auto & col : row)
      col = log(col);
  return res;
}

// applies Type-2 Discrete Cosine Transform to the powspec * filterbank result matrix
Vector2d<float> dct(Vector2d<float>& fbank, size_t numcepstra)
{
  Vector2d<float> res;
  size_t n = fbank[0].size();
  Vector<double> in_data(n);
  Vector<double> out_data(n);

  fftw_plan p = fftw_plan_r2r_1d(n,
        in_data.data(),
        out_data.data(),
        FFTW_REDFT10,
        FFTW_ESTIMATE);
  if (p == nullptr)
  {
    printf("Error: There was a problem initiating FFTW_PLAN\n");
    return res;
  }

  size_t finalSize = n < numcepstra ? n : numcepstra;
  for (size_t i = 0; i < fbank.size(); i++)
  {
    Vector<float> temp_in(n);
    std::copy(fbank[i].begin(),
        fbank[i].end(),
        temp_in.begin());

    for (size_t i = 0; i < n; i++)
      in_data[i] = temp_in[i];

    fftw_execute(p);

    Vector<double> temp_out(n);
    std::copy(out_data.begin(),
        out_data.end(),
        temp_out.begin());

    Vector<float> coeffs(finalSize);
    coeffs[0] = temp_out[0] * sqrt(0.25f / n);
    double ff = sqrt(0.5f / n);
    for (size_t j = 1; j < finalSize; j++)
      coeffs[j] = temp_out [j] * ff;

    res.push_back(coeffs);
  }
  fftw_destroy_plan(p);

  return res;
}

// calculates Dynamic Time Warping distance between two matrices
Pair<Vector2d<float>, Vector<Pair<size_t, size_t>>> dtw(
    Vector2d<float>& x,
    Vector2d<float>& y
)
{
  Pair<Vector2d<float>, Vector<Pair<size_t, size_t>>> res;
  if (x.size() == 0 || y.size() == 0)
  {
    printf("Error: Both matrices should have data when calculating DTW.\n");
    return res;
  }

  size_t r = x.size();
  size_t c = y.size();

  Vector2d<float> D0(r+1, Vector<float>(c+1, 0));
  for (size_t i = 1; i < c+1; i++)
    D0[0][i] = INFINITY;
  for (size_t i = 1; i < r+1; i++)
    D0[i][0] = INFINITY;

  Vector2d<float> D1;
  D1.reserve(D0.size());
  for (size_t i = 1; i < D0.size(); i++)
  {
    Vector<float> temp(D0[i].size() - 1);
    std::copy(D0[i].begin() + 1,
        D0[i].end(),
        temp.begin());
    D1.push_back(temp);
  }

  auto static distance = [](Vector<float>& a, Vector<float>& b) -> float {
    float sum = 0.0f;
    if (a.size() != b.size())
    {
      printf("Error: When calculating distances both arrays should be the same size.\n");
      return sum;
    }
    for (size_t i = 0; i < a.size(); i++)
      sum += std::abs(a[i] - b[i]);
    return sum;
  };

  for (size_t i = 0; i < r; i++)
  {
    for (size_t j = 0; j < c; j++)
    {
      D1[i][j] = distance(x[i], y[j]);
    }
  }

  Vector2d<float> cost_mtx = D1;

  Vector<float> min_list;
  min_list.reserve(3);
  for (size_t i = 0; i < r; i++)
  {
    for (size_t j = 0; j < c; j++)
    {
      size_t i_k = std::min<size_t>(i + 1, r);
      size_t j_k = std::min<size_t>(j + 1, c);

      min_list.push_back(D0[i][j]);
      min_list.push_back(D0[i_k][j]);
      min_list.push_back(D0[i][j_k]);

      D1[i][j] += *std::min_element(min_list.begin(), min_list.end());
      min_list.clear();
    }
  }

  Vector<Pair<size_t, size_t>> path;
  if (x.size() == 1)
  {
    for (size_t i = 0; i < y.size(); i++)
      path.push_back(Pair<size_t, size_t>(0, i));
  }
  else if (y.size() == 1)
  {
    for (size_t i = 0; i < x.size(); i++)
      path.push_back(Pair<size_t, size_t>(i, 0));
  }
  else
  {
    size_t i = D0.size() - 2;
    size_t j = D0[0].size() - 2;

    path = {Pair<size_t, size_t>(i, j)};
    while (i > 0 || j > 0)
    {
      Vector<float> temp = {D0[i][j], D0[i][j + 1], D0[i + 1][j]};
      auto min_element = std::min_element(temp.begin(), temp.end());
      size_t tb = std::distance(temp.begin(), min_element);

      if (tb == 0)
      {
        i -= 1;
        j -= 1;
      }
      else if (tb == 1)
        i -= 1;
      else
        j -= 1;
      path.insert(path.begin(), Pair<size_t, size_t>(i, j));
    }
  }

  float dist = D1[D1.size() - 1][D1[D1.size() - 1].size() - 1];
  res.first = cost_mtx;
  res.second = path;

  return res;
}

// normalizes the distance gotten from DTW function to range 0 - 1
float normalized_min_max(Vector2d<float>& cost_mtx, Vector<Pair<size_t, size_t>>& path)
{
  float min_val = INFINITY;
  float max_val = -INFINITY;
  for (size_t i = 0; i < cost_mtx.size(); i++)
  {
    for (size_t j = 0; j < cost_mtx[0].size(); j++)
    {
      if (min_val > cost_mtx[i][j])
        min_val = cost_mtx[i][j];
      if (max_val < cost_mtx[i][j])
        max_val = cost_mtx[i][j];
    }
  }

  float sum = 0.0f;
  for (auto &v : path)
  {
    size_t i = v.first;
    size_t j = v.second;
    sum += cost_mtx[i][j];
  }

  float avg_step = sum / path.size();
  return std::abs<float>(1.0f - (avg_step - min_val) / (max_val - min_val));
}
