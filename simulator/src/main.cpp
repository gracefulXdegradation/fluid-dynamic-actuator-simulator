#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "Frames.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"

using namespace Eigen;
using namespace std;

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

        std::cout << "q_it" << std::endl;
        std::cout << q_it[0] << std::endl;
        std::cout << "t_omega_t" << std::endl;
        std::cout << t_omega_t.col(0) << std::endl;
        std::cout << "i_r_gs" << std::endl;
        std::cout << i_r_gs.col(0) << std::endl;
        std::cout << "i_v_gs" << std::endl;
        std::cout << i_v_gs.col(0) << std::endl;

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
