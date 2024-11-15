#ifndef MATHHELPERS_HPP
#define MATHHELPERS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace MathHelpers
{
  inline int signum(double x)
  {
    return x >= 0.0 ? 1 : -1;
  }

  inline Vector3d normalizeVector(const Vector3d &vec)
  {
    return vec / vec.norm();
  }

  inline double wrapTo2Pi(double angle)
  {
    return fmod(angle + 2 * M_PI, 2 * M_PI);
  }

  inline double deg2rad(double deg)
  {
    return deg * M_PI / 180.0;
  }

  inline double rad2deg(double rad)
  {
    return rad * 180.0 / M_PI;
  }

  inline MatrixXd cross(const MatrixXd &a, const MatrixXd &b)
  {
    if (a.rows() != 3 || b.rows() != 3 || a.cols() != b.cols())
    {
      std::cerr << "Error: Matrices must be 3xN in size and have the same number of columns.\n";
      exit(1); // Exit if the matrices have incompatible sizes
    }

    // Initialize the result matrix with the same number of rows and columns
    MatrixXd x(3, a.cols());

    // Loop over each column and compute the cross product
    for (int i = 0; i < a.cols(); ++i)
    {
      // Extract the 3D vectors from the matrices (column i)
      Vector3d posVec = a.col(i);
      Vector3d velVec = b.col(i);

      // Compute the cross product and store it in the result matrix
      x.col(i) = posVec.cross(velVec);
    }

    return x;
  }

  inline std::vector<Matrix3d> reshapeAndPagetranspose(const MatrixXd &x, const MatrixXd &y, const MatrixXd &z)
  {
    int cols = x.cols();
    int rows = x.rows();

    std::vector<MatrixXd> reshaped_x, reshaped_y, reshaped_z;

    for (int i = 0; i < cols; ++i)
    {
      MatrixXd slice_x(rows, 1);
      MatrixXd slice_y(rows, 1);
      MatrixXd slice_z(rows, 1);

      for (int j = 0; j < rows; ++j)
      {
        slice_x(j, 0) = x(j, i);
        slice_y(j, 0) = y(j, i);
        slice_z(j, 0) = z(j, i);
      }

      reshaped_x.push_back(slice_x);
      reshaped_y.push_back(slice_y);
      reshaped_z.push_back(slice_z);
    }

    std::vector<Matrix3d> attitude_matrix_pages;

    for (int i = 0; i < cols; ++i)
    {
      // Combine slices along the third dimension (column)
      MatrixXd attitude_matrix(rows, rows);
      attitude_matrix.col(0) = reshaped_x[i];
      attitude_matrix.col(1) = reshaped_y[i];
      attitude_matrix.col(2) = reshaped_z[i];

      attitude_matrix_pages.push_back(attitude_matrix.transpose()); // Page transpose
    }

    return attitude_matrix_pages;
  }

}

#endif
