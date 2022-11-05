#include <string>

#include <includes/utils.h>
#include <includes/steps.h>

int main(int argc, char ** argv)
{
  if (argc < 3)
  {
    printf("Run with: %s IN_SAMPLE_PATH OUT_SAMPLE_PATH\n", argv[0]);
    printf("IN_SAMPLE_PATH - any raw (f32) sample with 44100 sample rate\n");
    printf("OUT_SAMPLE_PATH - raw (f32) sample after preemphasis\n");
    exit(EXIT_FAILURE);
  }
  std::string input = argv[1];
  std::string output = argv[2];

  Vector<float> original = read_data(input);
  if (original.size() == 0)
  {
    printf("Failed to read the sample\n");
    exit(EXIT_FAILURE);
  }
  Vector<float> preemphasised = preemph(original);
  if(!export_result(output, preemphasised))
  {
    printf("There was an issue exporting the result\n");
    exit(EXIT_FAILURE);
  }
  printf("Done\n");
}
