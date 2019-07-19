/*
 * File: lecm_parameters.cpp
 * Author: Guilherme O. Chagas
 *
 * @brief LECM [1] parameters implementation.
 * 
 * (I'm sorry for my bad english xD)
 *
 * Created on July 18, 2019, 04:42 PM
 * 
 * References:
 * [1] Y. Gao, H. Zhang and Y. Zhang. Overlapping community detection based on
 * conductance optimization in large-scale networks, Physica A 522 (2019) 69-79.
 */

#include "../headers/lecm_parameters.hpp"
#include <experimental/filesystem>


/////////////////////////////// Helper functions ///////////////////////////////

namespace
{

/**
 * @brief
 * @param
 * @return
 */
bool is_number(const std::string &str)
{
    try
    {
        std::stod(str);
    }
    catch(...)
    {
        return false;
    }
    return true;
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////


double Lecm_parameters::alpha() const
{
    return m_flag_val_map.at("-a");
}


double Lecm_parameters::beta() const
{
    return m_flag_val_map.at("-b");
}


double Lecm_parameters::chi() const
{
    return m_flag_val_map.at("-c");
}


double Lecm_parameters::epsilon() const
{
    return m_flag_val_map.at("-e");
}


double Lecm_parameters::theta() const
{
    return m_flag_val_map.at("-t");
}


std::string Lecm_parameters::get_graph_path() const
{
    return m_graph_path;
}


bool Lecm_parameters::set_parameters(const int argc, char** argv)
{
    if (argc < 3 || argc > 13 || (argc - 1) % 2 != 0) // params must be even
    {
        std::cerr << "[ERROR] Wrong number of parameters.\n";
        return false;
    }

    int argv_i = 1; // argv array index (even is flag and odd is value)
    // just flag "-f" is mandatory, others are optional
    if (std::string(argv[argv_i]) != "-f")
    {
        std::cerr << "[ERROR] Flag -f is mandatory.\n";
        return false;
    }

    ++argv_i; // even index to get the value
    if (!set_graph_path(argv[argv_i]))
    {
        std::cerr << "[ERROR] Graph file does not exists.\n";
        return false;
    }

    while ((argv_i += 2) < argc)
    {
        if (!set(argv[argv_i - 1], argv[argv_i]))
        {
            std::cerr << "[ERROR] Wrong parameter.\n";
            return false;
        }
    }

    return true;
}


/////////////////////////////// private methods ////////////////////////////////


bool Lecm_parameters::set(const std::string &flag, const std::string &val_str)
{
    auto it = m_flag_val_map.find(flag);
    if (it == m_flag_val_map.end()) // checks if the flag is valid
    {
        return false;
    }

    if (!is_number(val_str)) // checks whether the value is a valid number
    {
        return false;
    }

    double val_d = std::stod(val_str);
    if (val_d < 0 || val_d > 1) // checks whether the value is in the range
    {
        return false;
    }

    (*it).second = val_d;

    return true;
}


bool Lecm_parameters::set_graph_path(const std::string &path)
{
    if (!std::experimental::filesystem::exists(path))
    {
        return false;
    }
    m_graph_path = path;
    return true;
}