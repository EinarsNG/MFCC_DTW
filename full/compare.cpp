#include <string>

#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float FRAME_SIZE = 0.5f; // 500ms
constexpr float FRAME_STEP = 0.5f; // 500ms
constexpr size_t NUM_FILT = 26;

Vector2d<float> get_mfcc(Vector<float>& sample)
{
  Vector<float> preemphasised = preemph(sample);
  Vector2d<float> framed = framing(preemphasised, SAMPLE_RATE, FRAME_SIZE, FRAME_STEP);
  Vector2d<float> windowed(framed.size());
  for (size_t i = 0; i < framed.size(); i++)
    windowed[i] = hamming_window(framed[i]);
  size_t nfft = calculate_nfft(SAMPLE_RATE, FRAME_SIZE);
  Vector2d<float> dft_res = dft(windowed, nfft);
  Vector2d<float> power_spectrum = powspec(dft_res, nfft);
  Vector2d<float> fbank = filterbanks(NUM_FILT, nfft, SAMPLE_RATE);
  Vector2d<float> fbank_t = transpose(fbank);
  Vector2d<float> power_fbank = power_spectrum * fbank_t;
  Vector2d<float> cosine_transform = dct(power_fbank);
  return cosine_transform;
}

float compare(Vector<float>& a, Vector<float>& b)
{
  auto mfcc_one = get_mfcc(a);
  auto mfcc_two = get_mfcc(b);
  auto dtw_res = dtw(mfcc_one, mfcc_two);
  Vector2d<float> cost_mtx = dtw_res.first;
  Vector<Pair<size_t, size_t>> path = dtw_res.second;
  return normalized_min_max(cost_mtx, path);
}

int main(int argc, char ** argv)
{
  if (argc < 3)
  {
    printf("Run with: %s IN_SAMPLE_PATH IN_SAMPLE_PATH\n", argv[0]);
    printf("IN_SAMPLE_PATH - any raw (f32) sample with 44100 sample rate\n");
    exit(EXIT_FAILURE);
  }
  std::string first = argv[1];
  std::string second = argv[2];

  Vector<float> first_sample = read_data(first);
  if (first_sample.size() == 0)
  {
    printf("Failed to read the first sample\n");
    exit(EXIT_FAILURE);
  }

  Vector<float> second_sample = read_data(second);
  if (second_sample.size() == 0)
  {
    printf("Failed to read the second sample\n");
    exit(EXIT_FAILURE);
  }
  
  float similarity = compare(first_sample, second_sample);
  printf("Similarity between samples is: %f\n", similarity);
  
  printf("Done\n");
}
