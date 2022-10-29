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
		printf("IN_PATH - input path to power spectrum file.\n");
		printf("OUT_FOLDER_PATH - output path to power spectrum * filterbanks result.\n");
		exit(EXIT_FAILURE);
	}
	std::string input = argv[1];
	std::string output = argv[2];

	Vector<float> data = read_data(input);
	Vector2d<float> data_mtx;
	data_mtx.push_back(data);

	size_t nfft = calculate_nfft(SAMPLE_RATE, WINDOW_LENGHT);
	Vector2d<float> res = filterbanks(NUM_FILT, nfft, SAMPLE_RATE);

	Vector2d<float> dot_product = data_mtx * res;
	// since multiplying 1xN with NxM matrix produces a 1xN matrix anyways
	Vector<float> row = dot_product[0];

	if(!export_result(output, row))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
}
