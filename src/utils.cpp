#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <regex>

#include <includes/utils.h>

// reads a single raw sample (pcm) into memory
Vector<float> read_pcm(std::string filepath)
{
	std::ifstream ifs(filepath, std::ios::binary);
	if (!ifs.is_open())
	{
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

// reads mutliple raw samples into memory (Warning: Can be memory intensive,
// though for the sake of this example should be fine)
Vector2d<float> read_pcms(std::string folder)
{
	Vector2d<float> pcms;

	// we want to sort paths so they are in correct order: i.e. 1.pcm, 2.pcm, 3.pcm, ...
	Vector<std::string> filepaths;
	for (auto & entry : std::filesystem::directory_iterator(folder))
	{
		std::regex re_expr(".*[0-9]+\\.pcm$");
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
		std::regex re_expr(".*?([0-9]+)\\.pcm$");
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
		Vector<float> temp = read_pcm(entry);
		pcms.push_back(temp);
	}

	return pcms;
}

// exports a single pcm
bool export_pcm(std::string filepath, Vector<float>& data)
{
	std::ofstream ofs(filepath, std::ios::binary);
	if (!ofs.is_open()) return false;
	ofs.write((char*)data.data(), sizeof(float) * data.size());
	ofs.close();
	return true;
}

// exports multiple pcms
bool export_pcms(std::string folder, Vector2d<float>& data)
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
		std::filesystem::path filename(std::to_string(i+1) + ".pcm");
		auto final_path = out / filename;
		std::string out_str = std::string(final_path.string());
		if(!export_pcm(out_str, data[i])) res = false;
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

Vector<float> rowsum(Vector2d<float> powspectrum)
{
  Vector<float> res(powspectrum.size());
  for (size_t i = 0; i < powspectrum.size(); i++)
  {
    float sum = 0;
    for (auto & entry : powspectrum[i])
      sum += entry;
    res[i] = sum;
  }
  return res;
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

Vector<float> linearspace(float start, float end, size_t count)
{
  Vector<float> res;
  float step = (end - start) / static_cast<float>(count - 1);
  for (float i = start; i <= end; i += step)
    res.push_back(i);
  return res;
}
