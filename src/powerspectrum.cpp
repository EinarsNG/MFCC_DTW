#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float WINDOW_LENGHT = 0.5; // 500 ms

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		printf("Run with: %s IN_FOLDER_PATH OUT_FOLDER_PATH\n", argv[0]);
		printf("IN_FOLDER_PATH - contains raw (f32) samples from previous step (dft).\n");
		//printf("OUT_FOLDER_PATH - output folder for raw (f32) samples after applying power spectrum on them.\n");
		printf("OUT_PATH - output path of a file after calculating power spectrum.\n");
		exit(EXIT_FAILURE);
	}
	std::string input = argv[1];
	std::string output = argv[2];

	Vector2d<float> data = read_all_data(input);
	size_t nfft = calculate_nfft(SAMPLE_RATE, WINDOW_LENGHT);
	Vector2d<float> res = powspec(data, nfft);
  Vector<float> reduced = rowsum(res);
  for (auto & entry : reduced)
  {
    if (entry == 0.0f)
    {
      printf("None of row sums from power spectrum can be 0. This might be due to one of the frames (in framing step) being all 0)\n");
      exit(EXIT_FAILURE);
    }
  }
	
	//if (!export_pcms(output, res))
	//{
	//	printf("There was an issue exporting the result\n");
	//	exit(EXIT_FAILURE);
	//}

	printf("Done\n");
}
