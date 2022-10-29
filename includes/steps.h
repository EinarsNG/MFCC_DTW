#ifndef STEPS_H
#define STEPS_H
#include <vector>
#include <complex>
#include <includes/utils.h>

Vector<float> preemph(Vector<float> &data, float coeff = 0.97f);

Vector2d<float> framing(
		Vector<float> &data,
		size_t sample_rate,
		float frame_len,
		float frame_step
);

Vector<float> hamming_window(Vector<float> &data);

Vector2d<float> dft(
		Vector2d<float> frames,
		size_t nfft
);

Vector2d<float> powspec(Vector2d<float> frames, size_t nfft);

Vector2d<float> filterbanks(size_t numfilt, size_t nfft, size_t sample_rate);

Vector2d<float> log(Vector2d<float>& mtx);

Vector2d<float> dct(Vector2d<float>& feat, size_t numceptra = 13);
#endif
