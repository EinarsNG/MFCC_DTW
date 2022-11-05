#include <includes/utils.h>
#include <includes/steps.h>

int main(int argc, char ** argv)
{
  if (argc < 3)
  {
    printf("Run with: %s IN_FOLDER_PATH IN_FOLDER_PATH\n", argv[0]);
    printf("IN_FOLDER_PATH - contains result from previous step (DCT).\n");
    printf("IN_FOLDER_PATH - contains second result from previous step (DCT).\n");
    exit(EXIT_FAILURE);
  }
  std::string input = argv[1];
  std::string input_two = argv[2];

  Vector2d<float> data = read_all_data(input);
  if (data.size() == 0)
  {
    printf("Failed to read data from first input\n");
    exit(EXIT_FAILURE);
  }

  Vector2d<float> data_two = read_all_data(input_two);
  if (data.size() == 0)
  {
    printf("Failed to read data from second input\n");
    exit(EXIT_FAILURE);
  }

  auto res = dtw(data, data_two);
  Vector2d<float> cost_mtx = res.first;
  Vector<Pair<size_t, size_t>> path = res.second;
  float normalized = normalized_min_max(cost_mtx, path);
  printf("Similarity between samples is: %f\n", normalized);
}
