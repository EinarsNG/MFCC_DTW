#include <vector>
#include <string>

#include <includes/steps.h>
#include <includes/utils.h>

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		printf("Run with: %s IN_FOLDER_PATH OUT_FOLDER_PATH\n", argv[0]);
		printf("IN_FOLDER_PATH - contains raw (f32) samples from previous step (framing).\n");
		printf("OUT_FOLDER_PATH - output folder for raw (f32) samples after applying hamming window on them.\n");
		exit(EXIT_FAILURE);
	}
	std::string input = argv[1];
	std::string output = argv[2];
	
	Vector2d<float> data = read_pcms(input);
	Vector2d<float> res;
	for (auto & entry : data)
	{
		res.push_back(hamming_window(entry));
	}

	if(!export_results(output, res))
	{
		printf("There was an issue exporting the result\n");
		exit(EXIT_FAILURE);
	}
	printf("Done\n");
}
