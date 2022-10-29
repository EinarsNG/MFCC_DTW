#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

template<typename T>
using Vector = std::vector<T>;
template<typename T>
using Vector2d = std::vector<std::vector<T>>;

Vector<float> read_pcm(std::string filepath);
Vector2d<float> read_pcms(std::string folder);

bool export_pcm(std::string filepath, Vector<float>& data);
bool export_pcms(std::string folder, Vector2d<float>& data);

size_t calculate_nfft(size_t sample_rate, float window_length);
Vector<float> rowsum(Vector2d<float> powspectrum);
float hz_to_mel(size_t freq);
size_t mel_to_hz(float mel);
float mel_to_fft_bin(float mel, size_t nfft, size_t sample_rate);

Vector<float> linearspace(float start, float end, size_t count);
#endif
