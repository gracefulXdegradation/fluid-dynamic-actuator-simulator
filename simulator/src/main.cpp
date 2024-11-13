#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "Frames.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"
#include <cmath>

using namespace Eigen;
using namespace std;

// Checks whether combined contact criterion of access == 1 (satellite is
// inside cone drawn by ground station FOV) and distance < 1.5e6 m
// (ground station is in range of satellite laser communication terminal) is
// satisfied.
std::pair<Eigen::VectorXi, int64_t> contact_check(const VectorXi &access, const VectorXd &distance, const std::vector<std::chrono::_V2::system_clock::time_point> &date_times)
{
    // Ensure input sizes match
    if (access.size() != distance.size() || access.size() != date_times.size())
    {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    Eigen::VectorXi contact = Eigen::VectorXi::Zero(access.size());

    for (int k = 0; k < access.size(); k++)
    {
        contact[k] = access[k] == 1 && distance[k] < 1.5e3;
    }

    auto condition = [](int x)
    { return x == 1; };

    // Variables to store the first and last index satisfying the condition
    int first_index = -1;
    int last_index = -1;

    // Find the first element that satisfies the condition
    for (int i = 0; i < contact.size(); ++i)
    {
        if (condition(contact(i)))
        {
            first_index = i;
            break; // Exit after finding the first one
        }
    }

    // Find the last element that satisfies the condition
    for (int i = contact.size() - 1; i >= 0; --i)
    {
        if (condition(contact(i)))
        {
            last_index = i;
            break; // Exit after finding the last one
        }
    }

    int64_t duration = std::chrono::duration_cast<std::chrono::seconds>(date_times[last_index] - date_times[first_index]).count();

    std::cout << duration << std::endl;

    return std::make_pair(contact, duration);
}

// Function to calculate visibility angle
std::pair<VectorXd, VectorXi> visibility(const Matrix3Xd &i_r, const Matrix3Xd &i_r_gs, double angle)
{
    // Calculate the vector from GS to SC
    Matrix3Xd i_bard = i_r - i_r_gs;

    // Normalize the vectors column-wise
    Matrix3Xd norm_i_r_gs = i_r_gs.colwise().normalized();
    Matrix3Xd norm_i_bard = i_bard.colwise().normalized();

    // Compute the visibility angle (element-wise acos of dot products)
    VectorXd elrad(i_r.cols());
    for (int i = 0; i < i_r.cols(); ++i)
    {
        elrad(i) = acos(norm_i_r_gs.col(i).dot(norm_i_bard.col(i)));
    }

    // Check visibility condition
    VectorXi v = (elrad.array() < angle).cast<int>();

    // Convert elevation angles to degrees
    VectorXd el = elrad.array() * (180.0 / M_PI);

    return std::make_pair(el, v);
}

int main()
{

    // Create an instance of ConfigParser with the path to your config.json
    Config config(string(BUILD_OUTPUT_PATH) + "/config.json");

    // Generate discrete time points using the function
    auto date_times = DateTime::generateTimePoints(config.getStartDateTime(), config.getEndDateTime(), config.getControlTimeStep());

    // Output the time points
    cout << "Executing simulation" << endl;
    cout << "====================" << endl;
    cout << "Start date: " << DateTime::formatTime(date_times.at(0)) << endl;

    try
    {
        TLE tle = TLE::fromFile(string(BUILD_OUTPUT_PATH) + "/tle.txt");

        auto eccentricAnomalies = OrbitalMechanics::eccentricAnomaly(date_times, tle.getMeanAnomaly(), tle.getMeanMotion(), tle.getEccentricity(), tle.getEpoch());
        auto trueAnomalies = OrbitalMechanics::trueAnomaly(eccentricAnomalies, tle.getEccentricity());

        Matrix3Xd m_i_r(3, date_times.size());
        m_i_r.setZero();
        Matrix3Xd m_i_v(3, date_times.size());
        m_i_v.setZero();

        for (int i = 0; i < trueAnomalies.size(); i++)
        {
            auto &ta = trueAnomalies[i];
            auto [i_r, i_v] = OrbitalMechanics::keplerian2ijk(tle.getSemiMajorAxis(), tle.getEccentricity(), tle.getInclination(), tle.getArgumentOfPerigee(), ta, tle.getRightAscension());
            m_i_r.col(i) = i_r;
            m_i_v.col(i) = i_v;
        }

        auto [q_in, n_omega_n] = Frames::nadir_frame(m_i_r, m_i_v);
        auto [q_it, t_omega_t, i_r_gs, i_v_gs] = Frames::target_pointing_frame(m_i_r, m_i_v, config.getGroundStationPosition(), date_times);

        auto [elevation, access] = visibility(m_i_r, i_r_gs, MathHelpers::deg2rad(config.getGroundStationElevation()));
        auto distance = (i_r_gs - m_i_r).colwise().norm();
        // contact_check(access, distance, date_times);

        // ----- Attitude from commanded frame to inertial frame -----

        std::vector<Eigen::Quaterniond> q_ic(q_in.size());
        for (int i = 0; i < q_in.size(); i++)
        {
            q_ic[i] = access[i] == 1 ? q_it[i] : q_in[i];
        }

        // ----- Attitude expressed relative to nadir pointing frame -----
        std::vector<Eigen::Quaterniond> q_ni;
        q_ni.resize(q_in.size());
        std::vector<Eigen::Quaterniond> q_nc;
        q_nc.resize(q_in.size());

        for (size_t i = 0; i < q_ni.size(); i++)
        {
            q_ni[i] = q_in[i].conjugate();
            q_nc[i] = q_ni[i] * q_ic[0];
        }

        // Save to file
        auto ts = DateTime::getCurrentTimestamp();
        DB::writematrix(m_i_r, "./output/" + ts, "i_r.csv");
        DB::writematrix(m_i_v, "./output/" + ts, "i_v.csv");
        DB::serializeTimePointsToCSV(date_times, "./output/" + ts, "t.csv");

        std::cout << "Data saved to output.txt" << std::endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
