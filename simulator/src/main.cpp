#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "Frames.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"
#include "Database.hpp"
#include "Logger.hpp"
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/minima.hpp>

using namespace Eigen;
using namespace std;
namespace odeint = boost::numeric::odeint;
using json = nlohmann::json;

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

    Logger::info(std::to_string(duration));

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

using TetrahedronActuatorAssemblyStateExtractorFunction = std::function<Matrix<double, 8, 1>(const Matrix<double, 15, 1> &)>;

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

/**
 * RHS represents a customizable spacecraft
 * STATE_DERIVATIVE = RHS(TIME, STATE, ACTUATOR_STATE_EXTRACTOR,
 * ACTUATORS, ACTUATOR_PROPERTY_EXTRACTOR, INERTIA, INERTIAL_POSITION) defines
 * the kinematics and dynamics equations of a spacecraft with a rigid body with
 * INERTIA, where the STATE consists in any case of the attitude quaternion
 * and the angular rate vector of the rigid body. ACTUATOR_STATE_EXTRACTOR,
 * ACTUATORS, and ACTUATOR_PROPERTY_EXTRACTOR allow to inject variable
 * momentum exchange actuator models into the right-hand side function.
 */
Matrix<double, 15, 1> rhs(
    const int time,
    const Matrix<double, 15, 1> &state,
    const Vector4d &control,
    TetrahedronActuatorAssemblyStateExtractorFunction actuator_assembly_state_extractor,
    std::function<Matrix<double, 8, 1>(const double, const Matrix<double, 8, 1> &, const Vector4d &)> actuators,
    std::function<std::pair<Vector3d, Vector3d>(const Matrix<double, 8, 1> &,
                                                const Matrix<double, 8, 1> &)>
        actuator_property_extractor,
    const Matrix3d &inertia,
    const Vector3d &inertial_position)
{
    // Decompose state vector
    Vector4d attitude = state.segment(0, 4);
    Vector3d angular_rate = state.segment(4, 3);
    Matrix<double, 8, 1> actuator_state = actuator_assembly_state_extractor(state);

    // Get actuator state derivative from function given by input 'actuators'
    Matrix<double, 8, 1> actuator_state_derivative = actuators(time, actuator_state, control);

    // Extract actuator properties in body frame
    auto [actuator_angular_momentum, actuator_torque] = actuator_property_extractor(actuator_state, actuator_state_derivative);

    // Dynamics equation
    auto angular_rate_derivative = inertia.colPivHouseholderQr().solve(-actuator_torque - angular_rate.cross(inertia * angular_rate + actuator_angular_momentum));

    // Kinematics equation
    Matrix3d A_lin_dot;
    A_lin_dot << 0.0, angular_rate(2), -angular_rate(1),
        -angular_rate(2), 0.0, angular_rate(0),
        angular_rate(1), -angular_rate(0), 0.0;

    Matrix4d a;
    a(0, 0) = 0.0;
    a.block<1, 3>(0, 1) = -angular_rate.transpose();
    a.block<3, 1>(1, 0) = angular_rate;
    a.block<3, 3>(1, 1) = A_lin_dot;
    auto attitude_derivative = 0.5 * a * attitude;

    // Construct state derivative
    Matrix<double, 15, 1> state_derivative;
    state_derivative << attitude_derivative, angular_rate_derivative, actuator_state_derivative;

    return state_derivative;
};

int main(int argc, char* argv[])
{
    // Get simulation_id from command-line argument
    if (argc < 2) {
        Logger::error("Usage: " + std::string(argv[0]) + " <simulation_id>");
        Logger::error("Example: " + std::string(argv[0]) + " 123e4567-e89b-12d3-a456-426614174000");
        return 1;
    }
    
    std::string simulation_id = argv[1];
    
    // Get database connection string from environment variable
    const char* db_url = std::getenv("DATABASE_URL");
    if (!db_url) {
        Logger::error("Error: DATABASE_URL environment variable not set");
        Logger::error("Please set DATABASE_URL environment variable");
        Logger::error("Example: export DATABASE_URL=\"postgresql://user:pass@localhost/dbname\"");
        return 1;
    }
    
    std::string connection_string = db_url;
    
    // Initialize database connection
    Database* db = nullptr;
    try {
        db = new Database(connection_string);
        if (!db->testConnection()) {
            Logger::error("Error: Failed to connect to database");
            delete db;
            return 1;
        }
    } catch (const std::exception& e) {
        Logger::error("Error connecting to database: " + std::string(e.what()));
        return 1;
    }
    
    try
    {
        // Get simulation parameters from database
        Logger::info("Loading simulation parameters for ID: " + simulation_id);
        SimulationParams params = db->getSimulationParams(simulation_id);
        
        // Update status to running
        db->updateSimulationStatus(simulation_id, "running");
        
        // Generate discrete time points using the function
        std::vector<std::chrono::_V2::system_clock::time_point> date_times = 
            DateTime::generateTimePoints(params.start_date_time, params.end_date_time, params.control_time_step);

        // Output the time points
        Logger::info("Executing simulation");
        Logger::info("Start date: " + DateTime::formatTime(date_times.at(0)));
        
        // Create TLE from database parameters
        TLE tle(params.tle_line1, params.tle_line2);

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
        auto [q_it, t_omega_t, i_r_gs, i_v_gs] = Frames::target_pointing_frame(m_i_r, m_i_v, params.ground_station_lla, date_times);

        auto [elevation, access] = visibility(m_i_r, i_r_gs, MathHelpers::deg2rad(params.ground_station_elevation));
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
                                               const Matrix<double, 8, 1> &state_derivative) -> std::pair<Vector3d, Vector3d>
        {
            return TetrahedronActuatorAssemblyPropertyExtractor(state, state_derivative, actuator_alignment);
        };

        auto cost_function = [&](const double time, const Vector2d &initial_state, const double command, const double commanded_torque) -> double
        {
            return ActuatorCostFunction(actuator, time, initial_state, command, std::chrono::duration<double>(params.control_time_step).count(), commanded_torque);
        };

        double precision = 1e-2;

        auto command_finder = [&](const double time, const Vector2d &state0, const double required_torque) -> double
        {
            return ActuatorControlValueFinder(cost_function, time, state0, required_torque, precision);
        };

        // Control

        // ----- Boundary conditions -----
        // Initial states

        MatrixXd state(initial_state.rows(), date_times.size());
        for (int i = 0; i < date_times.size(); ++i)
        {
            state.col(i) = initial_state;
        }

        std::chrono::_V2::system_clock::time_point precision_date_times = date_times[0];
        Matrix<double, 15, 1> precision_state = state.col(0);

        Matrix3d Kp = 1 * inertia;
        Matrix3d Kd = 3 * Kp;

        Matrix3Xd b_omega_d(3, date_times.size());
        Matrix4Xd a_control_torque(4, date_times.size());
        Matrix3Xd b_control_torque(3, date_times.size());
        Matrix4Xd a_command(4, date_times.size());
        std::vector<Quaterniond> q_bc(date_times.size());
        
        // Initialize matrices for incremental metric computation
        Matrix3Xd euler_angles(3, date_times.size());
        Matrix3Xd ang_mom_body_frame(3, date_times.size());
        
        // Prepare timestamps once (they don't change during simulation)
        std::vector<int64_t> timestamps;
        timestamps.reserve(date_times.size());
        for (const auto& dt : date_times) {
            timestamps.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(
                dt.time_since_epoch()).count());
        }
        
        // Batch size for incremental writes
        const size_t BATCH_SIZE = 1000; // Smaller batch size for more frequent updates
        size_t last_written_index = 0;
        
        // Helper function to write metrics for a range of indices
        // Note: euler_angles and ang_mom_body_frame should already be computed incrementally
        auto writeMetricsBatch = [&](size_t start_idx, size_t end_idx) {
            if (start_idx >= end_idx) return;
            
            // Convert a_control_torque from Nm to mNm (multiply by 1e3) for storage
            Matrix4Xd a_control_torque_mNm = a_control_torque * 1e3;
            
            // Write to database (writeMetricsDirect expects full matrices and uses start_idx/end_idx internally)
            db->writeMetricsDirect(simulation_id, timestamps, euler_angles, ang_mom_body_frame,
                                   a_control_torque_mNm, a_command, state, distance, start_idx, end_idx);
            
            Logger::info("Written metrics batch: " + std::to_string(start_idx) + " to " + std::to_string(end_idx - 1) + 
                        " (" + std::to_string(end_idx - start_idx) + " metrics)");
        };
        
        for (int i = 0; i < date_times.size() - 1; i++)
        {
            // ----- Control error -----
            // Error quaternion
            Vector4d q_bc_coef = state.col(i).segment(0, 4);
            Quaterniond q_bc_(q_bc_coef[0], q_bc_coef[1], q_bc_coef[2], q_bc_coef[3]);
            q_bc[i] = q_bc_.conjugate() * q_ic[i];
            // Angular rate error
            b_omega_d.col(i) = q_bc[i] * c_omega_c.col(i) - state.col(i).segment(4, 3);
            // ----- PD control law -----
            // Control input
            double qs = q_bc[i].w();
            Vector3d qv = {q_bc[i].x(), q_bc[i].y(), q_bc[i].z()};
            b_control_torque.col(i) = -1 * (Kp * MathHelpers::signum(qs) * qv + Kd * b_omega_d.col(i));
            //  ----- Distribute control input to actuators -----
            a_control_torque.col(i) = inverse_actuator_alignment * b_control_torque.col(i);
            Matrix<double, 2, 4> actuator_state = TetrahedronActuatorAssemblyStateExtractor(state.col(i)).reshaped(2, 4);

            //  ----- Calculate the command to the actuator assembly -----
            for (int j = 0; j < 4; j++)
            {
                a_command(j, i) = command_finder(0.0, actuator_state.col(j), a_control_torque(j, i));
            }

            // ----- Dynamics and kinematics equation within rhs function -----
            // Solve differential equations for next time point or span

            std::vector<Matrix<double, 15, 1>> integrated_states;
            std::vector<double> integrated_time_steps;
            auto observer = [&](const Matrix<double, 15, 1> &x, const double t)
            {
                integrated_states.push_back(x);
                integrated_time_steps.push_back(t);
            };

            auto system = [&](const Matrix<double, 15, 1> &x, Matrix<double, 15, 1> &dxdt, const double t)
            {
                Matrix<double, 15, 1> state_derivative = rhs(t, x, a_command.col(i), TetrahedronActuatorAssemblyStateExtractor, actuator_assembly, actuator_property_extractor, inertia, m_i_r.col(i));
                dxdt.col(0) = state_derivative.col(0);
            };

            // ----- Solve with higher frequency than control cycle -----
            double t0 = std::chrono::duration_cast<std::chrono::milliseconds>(date_times[i].time_since_epoch()).count() / 1000.0;
            double t1 = std::chrono::duration_cast<std::chrono::milliseconds>(date_times[i + 1].time_since_epoch()).count() / 1000.0;
            double dt = (t1 - t0) / 11;

            Matrix<double, 15, 1> x0 = state.col(i);

            odeint::integrate_const(odeint::runge_kutta4<Matrix<double, 15, 1>>(), system, x0, t0, t1, dt, observer);
            state.col(i + 1) = integrated_states.back();
            
            // Compute euler_angles for current step (needed for batch write)
            q_bc[i].normalize();
            euler_angles.col(i) = q_bc[i].toRotationMatrix().eulerAngles(2, 0, 2);
            
            // Compute ang_mom_body_frame for current step
            std::vector<int> indices = {7, 9, 11, 13};
            Vector4d selected_state;
            selected_state << state(indices[0], i), state(indices[1], i), 
                            state(indices[2], i), state(indices[3], i);
            ang_mom_body_frame.col(i) = actuator_alignment * selected_state;
            
            // Write metrics incrementally in batches
            // Note: At iteration i, we've computed metrics for step i, so we can write up to step i+1
            // We write steps [last_written_index, i+1), which means we write steps last_written_index through i
            if ((i + 1) % BATCH_SIZE == 0) {
                size_t batch_end = i + 1; // Write up to (but not including) i+1, which means we write step i
                writeMetricsBatch(last_written_index, batch_end);
                last_written_index = batch_end;
            }
        }

        // Calculate the Final State
        int last_index = date_times.size() - 1;
        // ----- Control error -----
        Vector4d q_bc_coef = state.col(last_index).segment(0, 4);
        Quaterniond q_bc_(q_bc_coef[0], q_bc_coef[1], q_bc_coef[2], q_bc_coef[3]);
        q_bc[last_index] = q_bc_.conjugate() * q_ic[last_index];
        b_omega_d.col(last_index) = q_bc[last_index] * c_omega_c.col(last_index) - state.col(last_index).segment(4, 3);
        // ----- PD control law -----
        double qs = q_bc[last_index].w();
        Vector3d qv = {q_bc[last_index].x(), q_bc[last_index].y(), q_bc[last_index].z()};
        b_control_torque.col(last_index) = -1 * (Kp * MathHelpers::signum(qs) * qv + Kd * b_omega_d.col(last_index));
        //  ----- Distribute control input to actuators -----
        a_control_torque.col(last_index) = inverse_actuator_alignment * b_control_torque.col(last_index);
        Matrix<double, 2, 4> actuator_state = TetrahedronActuatorAssemblyStateExtractor(state.col(last_index)).reshaped(2, 4);
        // ----- Calculate the command to the actuator assembly -----
        for (int j = 0; j < 4; j++)
        {
            a_command(j, last_index) = command_finder(0.0, actuator_state.col(j), a_control_torque(j, last_index));
        }
        
        // Compute euler_angles and ang_mom_body_frame for final index
        q_bc[last_index].normalize();
        euler_angles.col(last_index) = q_bc[last_index].toRotationMatrix().eulerAngles(2, 0, 2);
        std::vector<int> indices = {7, 9, 11, 13};
        Vector4d selected_state;
        selected_state << state(indices[0], last_index), state(indices[1], last_index), 
                        state(indices[2], last_index), state(indices[3], last_index);
        ang_mom_body_frame.col(last_index) = actuator_alignment * selected_state;
        
        Logger::info("Simulation complete!");

        // Write any remaining metrics that haven't been written yet
        if (last_written_index < date_times.size()) {
            Logger::info("Writing final metrics batch...");
            writeMetricsBatch(last_written_index, date_times.size());
        }
        
        // Update status to completed
        db->updateSimulationStatus(simulation_id, "completed");
        
        Logger::info("First Date " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(date_times[0].time_since_epoch()).count()));
        Logger::info("Last Date " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(date_times[date_times.size() - 1].time_since_epoch()).count()));
        Logger::info("Simulation completed and data saved to database!");
    }
    catch (const std::exception &e)
    {
        Logger::error("Error: " + std::string(e.what()));
        if (db) {
            try {
                db->updateSimulationStatus(simulation_id, "failed", e.what());
            } catch (...) {
                // Ignore errors when updating failed status
            }
        }
        if (db) delete db;
        return 1;
    }
    
    if (db) delete db;
    return 0;
}
