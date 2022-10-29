#include <includes/utils.h>
#include <includes/steps.h>

#include <fstream>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float WINDOW_LENGHT = 0.5; // 500 ms
constexpr size_t NUM_FILT = 26; // filter count

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		printf("Run with: %s OUT_FOLDER_PATH\n", argv[0]);
		printf("OUT_FOLDER_PATH - output folder path for Mel's filterbanks.\n");
		exit(EXIT_FAILURE);
	}
	std::string output = argv[1];

	size_t nfft = calculate_nfft(SAMPLE_RATE, WINDOW_LENGHT);
  Vector2d<float> res = filterbanks(NUM_FILT, nfft, SAMPLE_RATE);

	if(!export_results(output, res))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
}
