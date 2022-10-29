#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float WINDOW_LENGHT = 0.5; // 500 ms
constexpr size_t NUM_FILT = 26; // filter count

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		printf("Run with: %s IN_FOLDER_PATH OUT_FOLDER_PATH\n", argv[0]);
		printf("IN_FOLDER_APATH - input folder containing result from previous step (power spectrum).\n");
		printf("OUT_FOLDER_PATH - output path to power spectrum * filterbanks result.\n");
		exit(EXIT_FAILURE);
	}
	std::string input = argv[1];
	std::string output = argv[2];

	Vector2d<float> data = read_all_data(input);

	size_t nfft = calculate_nfft(SAMPLE_RATE, WINDOW_LENGHT);
	Vector2d<float> res = filterbanks(NUM_FILT, nfft, SAMPLE_RATE);

  Vector2d<float> res_t = transpose(res);
	Vector2d<float> dot_product = data * res_t;

	if(!export_results(output, dot_product))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
}
