#include <includes/utils.h>
#include <includes/steps.h>

constexpr size_t NUM_CEPSTRA = 13;

int main(int argc, char ** argv)
{
  if (argc < 3)
  {
    printf("Run with: %s IN_FOLDER_PATH OUT_FOLDER_PATH\n", argv[0]);
    printf("IN_FOLDER_PATH - contains result from previous step (logagarithm).\n");
    printf("OUT_FOLDER_PATH - output folder containing all with Type2-DCT applied to them.\n");
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

  Vector2d<float> result = dct(data, NUM_CEPSTRA);

  if (!export_results(output, result))
  {
    printf("There was an issue exporting the result\n");
    exit(EXIT_FAILURE);
  }

  printf("Done\n");
}
