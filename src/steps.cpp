#include <mutex>
#include <fftw3.h>

#include <includes/steps.h>

Vector<float> preemph(Vector<float> &data, float coeff)
{
	std::vector<float> res(data.size());
	for (size_t i = data.size() - 1; i >= 1; i--)
        res[i] = 1.0f * data[i] - data[i-1] * coeff;
    res[0] = 1.0f * data[0];
	return res;
}

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

Vector<float> hamming_window(std::vector<float>& data)
{
	Vector<float> res(data.size());
	for (size_t i = 0; i < data.size(); i++)
	{
		res[i] = data[i] * (0.54f - 0.46f * cos(2.0f * M_PI * (static_cast<float>(i)/(static_cast<float>(data.size()-1.0f)))));
	}
	return res;
}

static std::mutex fft_mtx;
Vector2d<float> dft(
		Vector2d<float> frames,
		size_t nfft
)
{
    size_t new_col_size = nfft / 2 + 1;
    Vector2d<float> res;
    Vector<double> in_data(nfft);
    Vector<std::complex<double>> out_data(new_col_size);

    fftw_plan p = nullptr;
    {
        std::lock_guard<std::mutex> lock(fft_mtx);
        p = fftw_plan_dft_r2c_1d(nfft,
                                 in_data.data(),
                                 (fftw_complex *)out_data.data(),
                                 FFTW_ESTIMATE);
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
    {
        std::lock_guard<std::mutex> lock(fft_mtx);
        fftw_destroy_plan(p);
    }

    return res;
}

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
