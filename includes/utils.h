#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

template<typename T>
using Vector = std::vector<T>;
template<typename T>
using Vector2d = std::vector<std::vector<T>>;
template<typename T, typename N>
using Pair = std::pair<T, N>;

Vector<float> read_data(std::string filepath);
Vector2d<float> read_all_data(std::string folder);

// unused by the program, but present if there is a need
// to get data at each stage
bool export_csv(std::string filepath, Vector<float>& data);
bool export_csv(std::string filepath, Vector2d<float>& data);

bool export_result(std::string filepath, Vector<float>& data);
bool export_results(std::string folder, Vector2d<float>& data);

size_t calculate_nfft(size_t sample_rate, float window_length);
float hz_to_mel(size_t freq);
size_t mel_to_hz(float mel);
float mel_to_fft_bin(float mel, size_t nfft, size_t sample_rate);

Vector<float> linearspace(float start, float end, size_t count);
Vector2d<float> transpose(Vector2d<float>& matrix);
Vector2d<float> operator*(Vector2d<float>& a, Vector2d<float>& b);
#endif
