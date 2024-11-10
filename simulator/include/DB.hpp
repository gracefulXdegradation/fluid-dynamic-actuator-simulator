#ifndef DB_HPP
#define DB_HPP

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace DB
{

  // Templated function to write a matrix (vector of arrays) to a file in a specified directory
  template <std::size_t N>
  void writematrix(const std::vector<std::array<double, N>> &matrix, const std::string &directory, const std::string &filename)
  {
    // Ensure the directory exists, create it if it does not
    if (!fs::exists(directory))
    {
      fs::create_directories(directory);
    }

    // Construct the full file path
    std::string filePath = directory + "/" + filename;

    // Open file in write mode
    std::ofstream file(filePath);
    if (!file.is_open())
    {
      std::cerr << "Error: Could not open file for writing: " << filePath << std::endl;
      return;
    }

    // Write the matrix data to the file in CSV format
    for (const auto &row : matrix)
    {
      for (std::size_t i = 0; i < row.size(); ++i)
      {
        file << row[i];
        if (i < row.size() - 1)
        {
          file << ","; // Separate columns with commas
        }
      }
      file << "\n"; // Newline after each row
    }

    // Check if the file write was successful
    file.close();
    if (file.good())
    {
      std::cout << "File written successfully to " << filePath << std::endl;
    }
    else
    {
      std::cerr << "Error: Failed to write to file: " << filePath << std::endl;
    }
  }

}

#endif
