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
  template <size_t N>
  std::vector<std::vector<double>> transpose(const std::vector<std::array<double, N>> &matrix)
  {
    // Number of rows in the original matrix
    size_t M = matrix.size();

    // Create a new matrix with N rows and M columns
    std::vector<std::vector<double>> transposed(N, std::vector<double>(M));

    for (size_t i = 0; i < N; ++i)
    {
      for (size_t j = 0; j < M; ++j)
      {
        transposed[i][j] = matrix[j][i];
      }
    }

    return transposed;
  }

  void writematrix(const std::vector<std::vector<double>> &matrix, const std::string &directory, const std::string &filename)
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

  void serializeTimePointsToCSV(const std::vector<std::chrono::_V2::system_clock::time_point> &time_points, const std::string &directory, const std::string &filename)
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

    // Loop over the time_points vector
    for (const auto &tp : time_points)
    {
      // Convert time_point to time_t (Unix timestamp)
      std::time_t time_since_epoch = std::chrono::system_clock::to_time_t(tp);

      // Write the timestamp to the file (as Unix timestamp)
      file << time_since_epoch << ",";
    }
    file << "\n";

    // Close the file
    file.close();
    std::cout << "Timestamps serialized to " << filename << std::endl;
  }

}

#endif
