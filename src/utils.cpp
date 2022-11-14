#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <regex>

#include <includes/utils.h>

// reads a single file into memory
Vector<float> read_data(std::string filepath)
{
  std::ifstream ifs(filepath, std::ios::binary);
  if (!ifs.is_open())
  {
    printf("Failed to read file %s: %s\n", filepath.c_str(), strerror(errno));
    return Vector<float>();
  }
  ifs.seekg(0, std::ios::end);
  size_t length = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  Vector<float> res(length / sizeof(float));
  ifs.read((char*)res.data(), length);
  ifs.close();
  return res;
}

// reads mutliple files into memory (Warning: Can be memory intensive,
// though for the sake of this example should be fine)
Vector2d<float> read_all_data(std::string folder)
{
  Vector2d<float> pcms;

  // we want to sort paths so they are in correct order: i.e. 1.pcm, 2.pcm, 3.pcm, ...
  Vector<std::string> filepaths;

  for (auto & entry : std::filesystem::directory_iterator(folder))
  {
    std::regex re_expr(".*[0-9]+\\.tmp$");
    std::smatch re_match;
    std::string path = std::string(entry.path());
    if (std::regex_search(path, re_match, re_expr))
    {
      std::string path = std::string(entry.path());
      filepaths.push_back(path);
    }
  }

  auto predicate = [](const std::string& a, const std::string& b) -> bool
  {
    std::regex re_expr(".*[0-9]+\\.tmp$");
    std::string num_a;
    std::string num_b;

    std::smatch re_match;
    if (std::regex_search(a, re_match, re_expr))
    {
      num_a = re_match[1].str();
    }
    else
    {
      num_a = "0";
    }

    if (std::regex_search(b, re_match, re_expr))
    {
      num_b = re_match[1].str();
    }
    else
    {
      num_b = "0";
    }
    return atoll(num_a.c_str()) < atoll(num_b.c_str());
  };
  std::sort(filepaths.begin(), filepaths.end(), predicate);

  for (auto & entry : filepaths)
  {
    Vector<float> temp = read_data(entry);
    if (temp == Vector<float>())
      return Vector2d<float>();
    pcms.push_back(temp);
  }

  return pcms;
}

// exports a single column
bool export_csv(std::string filepath, Vector<float>& data)
{
  std::ofstream ofs(filepath, std::ios::binary);
  if (!ofs.is_open())
  {
    printf("Failed to write file %s: %s\n", filepath.c_str(), strerror(errno));
    return false;
  }
  for (size_t i = 0; i < data.size(); i++)
  {
    ofs << data[i] << "\n";
  }
  ofs.close();
  return true;
}

// exports multiple columns
bool export_csv(std::string filepath, Vector2d<float>& data)
{
  std::ofstream ofs(filepath, std::ios::binary);
  if (!ofs.is_open())
  {
    printf("Failed to write file %s: %s\n", filepath.c_str(), strerror(errno));
    return false;
  }
  for (size_t i = 0; i < data[0].size(); i++)
  {
    for (size_t j = 0; j < data.size(); j++)
    {
      ofs << data[j][i];
      if (j < data.size() - 1)
        ofs << ",";
    }
    ofs << "\n";
  }
  ofs.close();
  return true;
}

// exports a single result
bool export_result(std::string filepath, Vector<float>& data)
{
  std::ofstream ofs(filepath, std::ios::binary);
  if (!ofs.is_open())
  {
    printf("Failed to write file %s: %s\n", filepath.c_str(), strerror(errno));
    return false;
  }
  ofs.write((char*)data.data(), sizeof(float) * data.size());
  ofs.close();
  return true;
}

// exports multiple results
bool export_results(std::string folder, Vector2d<float>& data)
{
  std::filesystem::path out(folder);
  std::filesystem::directory_entry dir(out);
  if (!dir.exists())
  {
    std::filesystem::create_directory(dir);
  }
  bool res = true;
  for (size_t i = 0; i < data.size(); i++)
  {
    std::filesystem::path filename(std::to_string(i+1) + ".tmp");
    auto final_path = out / filename;
    std::string out_str = std::string(final_path.string());
    if(!export_result(out_str, data[i]))
    {
      res = false;
      break;
    }
  }
  return res;
}

size_t calculate_nfft(size_t sample_rate, float window_length)
{
  size_t window_samples = window_length * static_cast<float>(sample_rate);
  size_t nfft = 1;
  while (nfft < window_samples)
    nfft *= 2;
  return nfft;
}

float hz_to_mel(size_t freq)
{
  return 2595.0f * log10(1.0f + static_cast<float>(freq) / 700.0f);
}

size_t mel_to_hz(float mel)
{
  return 700.0f * (pow(10.0f , mel / 2595.0f) - 1.0f);
}

float mel_to_fft_bin(float mel, size_t nfft, size_t sample_rate)
{
  float n = nfft;
  float hz = mel_to_hz(mel);
  float sr = sample_rate;
  return (n + 1.0f) * hz / sr;
}

// evenly divides space, resulting a vector of points
Vector<float> linearspace(float start, float end, size_t count)
{
  Vector<float> res;
  float step = (end - start) / static_cast<float>(count - 1);
  for (float i = start; i <= end; i += step)
    res.push_back(i);
  return res;
}

// transposes the matrix
Vector2d<float> transpose(Vector2d<float>& matrix)
{
  Vector2d<float> res(matrix[0].size(), Vector<float>(matrix.size(), 0));
  for (size_t i = 0; i < matrix[0].size(); i++)
  {
    for (size_t j = 0; j < matrix.size(); j++)
    {
      res[i][j] = matrix[j][i];
    }
  }
  return res;
}

// enables matrix multiplication
Vector2d<float> operator*(Vector2d<float>& a, Vector2d<float>& b)
{
  Vector2d<float> res(a.size(), Vector<float>(b[0].size(), 0));
  if (a.size() == 0 || a[0].size() == 0 || b.size() == 0 || b[0].size() == 0)
  {
    printf("Error: One of the matrices is empty (A[%lu,%lu], B[%lu,%lu])\n", a.size(), a[0].size(), b.size(), b[0].size());
    return Vector2d<float>();
  }
  if (a[0].size() != b.size())
  {
    printf("Error: Matrix multiplication failed 'a' column size does not match 'b' row size (A[%lu,%lu], B[%lu,%lu])\n", a.size(), a[0].size(), b.size(), b[0].size());
    return Vector2d<float>();
  }
  for (size_t i = 0; i < a.size(); i++)
  {
    for (size_t k = 0; k < b.size(); k++)
    {
      for (size_t j = 0; j < b[0].size(); j++)
        res[i][j] += a[i][k] * b[k][j];
    }
  }
  return res;
}
