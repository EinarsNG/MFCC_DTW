#include <string>

#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t SAMPLE_RATE = 44100;
constexpr float FRAME_SIZE = 0.5;
constexpr float FRAME_STEP = 0.5; // same of frame size in this case, but can be smaller

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		printf("Run with: %s IN_SAMPLE_PATH OUT_FOLDER_PATH\n", argv[0]);
		printf("IN_SAMPLE_PATH - contains raw (f32) sample from previous step (preemphasis).\n");
		printf("OUT_FOLDER_PATH - output folder for raw (f32) samples after splitting them into frames.\n");
		exit(EXIT_FAILURE);
	}

	std::string input = argv[1];
	std::string output = argv[2];

	Vector<float> data = read_data(input);
	if (data.size() == 0)
	{
		printf("Failed to read the sample\n");
		exit(EXIT_FAILURE);
	}

	Vector2d<float> res = framing(data, SAMPLE_RATE, FRAME_SIZE, FRAME_STEP);
	if (!export_results(output, res))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
}
