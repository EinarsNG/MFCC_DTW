#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float WINDOW_LENGHT = 0.5; // 500 ms

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		printf("Run with: %s IN_FOLDER_PATH OUT_FOLDER_PATH\n", argv[0]);
		printf("IN_FOLDER_PATH - contains result from previous step (power spectrum * power spectrum).\n");
		printf("OUT_FOLDER_PATH - output folder containing all elements multiplied by logarithm.\n");
		exit(EXIT_FAILURE);
	}
	std::string input = argv[1];
	std::string output = argv[2];

	Vector2d<float> data = read_all_data(input);
  if (data.size() == 0)
	{
		printf("Failed to read data\n");
		exit(EXIT_FAILURE);
	}
  data = log(data);

  if(!export_results(output, data))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}

  printf("Done\n");
}
