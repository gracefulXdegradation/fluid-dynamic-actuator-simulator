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
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/minima.hpp>

using namespace Eigen;
using namespace std;
namespace odeint = boost::numeric::odeint;

// Checks whether combined contact criterion of access == 1 (satellite is
// inside cone drawn by ground station FOV) and distance < 1.5e6 m
// (ground station is in range of satellite laser communication terminal) is
// satisfied.
std::pair<VectorXi, int64_t> contact_check(const VectorXi &access, const VectorXd &distance, const std::vector<std::chrono::_V2::system_clock::time_point> &date_times)
{
    // Ensure input sizes match
    if (access.size() != distance.size() || access.size() != date_times.size())
    {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    VectorXi contact = VectorXi::Zero(access.size());

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

Vector2d SimpleFluidDynamicActuator(const double time, const Vector2d &actuator_state, double control, double gain, double time_constant)
{
    Vector2d state_derivative;

    // Calculating the derivatives
    state_derivative(0) = actuator_state(1) + gain / time_constant * control;
    state_derivative(1) = -actuator_state(1) / time_constant - gain / (time_constant * time_constant) * control;

    return state_derivative;
};

using ModelFunction = std::function<Vector2d(const double, const Vector2d &, double)>;

// TetrahedronActuatorAssembly function
Matrix<double, 8, 1> TetrahedronActuatorAssembly(const double time, const Matrix<double, 8, 1> &state, const Vector4d &control, ModelFunction model)
{
    Matrix<double, 8, 1> state_derivative = Matrix<double, 8, 1>::Zero();

    for (int i = 0; i < 4; i++)
    {
        int start_index = 2 * i;
        Vector2d actuator_state = state.segment<2>(start_index);
        Vector2d actuator_derivative = model(time, actuator_state, control(i));
        state_derivative.segment<2>(start_index) = actuator_derivative;
    }

    return state_derivative;
}

Matrix<double, 8, 1> TetrahedronActuatorAssemblyStateExtractor(const Matrix<double, 15, 1> &state)
{
    return state.segment(7, 8);
}

std::pair<Vector3d, Vector3d> TetrahedronActuatorAssemblyPropertyExtractor(
    const Matrix<double, 8, 1> &state,
    const Matrix<double, 8, 1> &state_derivative,
    const Matrix<double, 3, 4> &alignment)
{
    // Extract elements 1, 3, 5, and 7 (0-based indices 0, 2, 4, 6)
    Vector4d state_part, state_derivative_part;
    state_part << state(0), state(2), state(4), state(6);
    state_derivative_part << state_derivative(0), state_derivative(2), state_derivative(4), state_derivative(6);

    Vector3d b_angular_momentum = alignment * state_part;
    Vector3d b_torque = alignment * state_derivative_part;

    return std::make_pair(b_angular_momentum, b_torque);
}

double ActuatorCostFunction(
    const ModelFunction &model,
    const double time,
    const Vector2d &initial_state,
    double command,
    double time_step,
    double commanded_torque)
{
    double final_angular_momentum = 0.0;
    auto observer = [&](const Vector2d &x, const double t)
    {
        final_angular_momentum = x[0];
    };

    auto system = [&](const Vector2d &x, Vector2d &dxdt, const double t)
    {
        Vector2d state_derivative = model(t, x, command);
        dxdt[0] = state_derivative[0];
        dxdt[1] = state_derivative[1];
    };

    // Define the error stepper
    typedef odeint::runge_kutta_cash_karp54<Vector2d> error_stepper_rkck54;

    // Error bounds
    double err_abs = 1.0e-10; // consider 1e-6
    double err_rel = 1.0e-6;  // consider 1e-3
    double dt = 0.01;

    // there must be a better way to do it
    Vector2d x0 = {initial_state[0], initial_state[1]};

    odeint::integrate_adaptive(odeint::make_controlled(err_abs, err_rel, error_stepper_rkck54()), system, x0, time, time + time_step, dt, observer);

    double commanded_angular_momentum = initial_state[0] + commanded_torque * time_step;
    double cost = (commanded_angular_momentum - final_angular_momentum) * (commanded_angular_momentum - final_angular_momentum);

    return cost;
}

// function command = ActuatorControlValueFinder( ...
//     cost_function, ...
//     time, ...
//     initial_state, ...
//     required_torque, ...
//     options, ...
//     precision ...
// )
//     command = fminbnd( ...
//         @(command) cost_function( ...
//             time, ...
//             initial_state, ...
//             command, ...
//             required_torque ...
//         ), ...
//         -1, ...
//         1, ...
//         options ...
//     );
//     command = round(command / precision) * precision;
// end

using CostFunction = std::function<double(const double, const Vector2d &, double, double)>;

double ActuatorControlValueFinder(
    const CostFunction &cost_function,
    const double time,
    const Vector2d &initial_state,
    double required_torque,
    double precision)
{
    auto cost_functor = [&](const double comm) -> double
    {
        return cost_function(time, initial_state, comm, required_torque);
    };

    double lower = -1.0;
    double upper = 1.0;

    // Specify the precision (e.g., 53 bits for double precision)
    int bits = 53;

    // Call Brent's minimization method
    auto result = boost::math::tools::brent_find_minima(cost_functor, lower, upper, bits);
    double command = result.first;

    return std::round(command / precision) * precision;
};

int main()
{
    // Create an instance of ConfigParser with the path to your config.json
    Config config(string(BUILD_OUTPUT_PATH) + "/config.json");

    // Generate discrete time points using the function
    std::vector<std::chrono::_V2::system_clock::time_point> date_times = DateTime::generateTimePoints(config.getStartDateTime(), config.getEndDateTime(), config.getControlTimeStep());

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

        std::vector<Quaterniond> q_ic(q_in.size());
        for (int i = 0; i < q_in.size(); i++)
        {
            q_ic[i] = access[i] == 1 ? q_it[i] : q_in[i];
        }

        // ----- Attitude expressed relative to nadir pointing frame -----
        std::vector<Quaterniond> q_ni;
        q_ni.resize(q_in.size());
        std::vector<Quaterniond> q_nc;
        q_nc.resize(q_in.size());

        for (size_t i = 0; i < q_ni.size(); i++)
        {
            q_ni[i] = q_in[i].conjugate();
            q_nc[i] = q_ni[i] * q_ic[0];
        }

        // ----- Commanded angular rates in commanded frame -----
        // Two different commanded frames (nadir, target) used as reference
        Matrix3Xd c_omega_c(3, date_times.size());
        for (int i = 0; i < c_omega_c.cols(); i++)
        {
            c_omega_c.col(i) = access[i] == 1 ? t_omega_t.col(i) : n_omega_n.col(i);
        }
        // Commanded angular rates with respect to nadir frame
        // n_omega_c = rotateframe(q_nc, c_omega_c')';

        // ----- Spacecraft setup -----

        Matrix3d inertia = Matrix3d::Zero();
        inertia.diagonal().setConstant(0.002);

        double gamma = atan(sqrt(2.0));
        MatrixXd actuator_alignment(3, 4);
        actuator_alignment << sin(gamma), -sin(gamma), 0, 0,
            0, 0, sin(gamma), -sin(gamma),
            cos(gamma), cos(gamma), -cos(gamma), -cos(gamma);
        MatrixXd inverse_actuator_alignment = actuator_alignment.completeOrthogonalDecomposition().pseudoInverse();

        // ----- Initial states -----
        // Spacecraft states
        // Vector4d initial_attitude = q_ic[0].coeffs();
        Vector3d initial_angular_rate = n_omega_n.col(0);
        Matrix<double, 8, 1> initial_actuator_state = Matrix<double, 8, 1>::Zero();
        Matrix<double, 15, 1> initial_state;
        // TODO there must be a better way to change representation of quaternions.
        // Here, I am mimicing the MatLab's representation which is wxyz instead of xyzw used by Eigen
        initial_state << q_ic[0].w(), q_ic[0].x(), q_ic[0].y(), q_ic[0].z(), initial_angular_rate, initial_actuator_state;

        // ----- Fluid-dynamic actuator setup -----

        double Ka = 120.236e-6;
        double Ta = Ka / 861.584e-6;
        auto actuator = [&](const double time, const Vector2d &actuator_state, double control) -> Vector2d
        {
            return SimpleFluidDynamicActuator(time, actuator_state, control, Ka, Ta);
        };

        auto actuator_assembly = [&](const double time, const Matrix<double, 8, 1> &state, const Vector4d &command) -> Matrix<double, 8, 1>
        {
            return TetrahedronActuatorAssembly(time, state, command, actuator);
        };

        auto actuator_property_extractor = [&](const Matrix<double, 8, 1> &state,
                                               const Matrix<double, 8, 1> &state_derivative)
        {
            return TetrahedronActuatorAssemblyPropertyExtractor(state, state_derivative, actuator_alignment);
        };

        auto cost_function = [&](const double time, const Vector2d &initial_state, const double command, const double commanded_torque) -> double
        {
            return ActuatorCostFunction(actuator, time, initial_state, command, std::chrono::duration<double>(config.getControlTimeStep()).count(), commanded_torque);
        };

        // ===============

        Vector2d x0 = {0.0, 0.0};
        double command = -0.23606797749979;
        double commanded_torque = 0.0;
        auto cost = cost_function(0, x0, command, commanded_torque);

        std::cout << "Cost: " << cost << std::endl;

        // ===============

        double precision = 1e-2;

        double cost_time = 0.0;
        double req_torque = 2.55105652012911e-05;

        auto command_finder = [&](const double time, const Vector2d &initial_state, const double required_torque) -> double
        {
            return ActuatorControlValueFinder(cost_function, time, initial_state, required_torque, precision);
        };

        std::cout << "Minimal value found: " << command_finder(cost_time, x0, req_torque) << std::endl;

        // ===============

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
