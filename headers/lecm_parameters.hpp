/*
 * File: lecm_parameters.hpp
 * Author: Guilherme O. Chagas
 *
 * @brief LECM [1] parameters header.
 * 
 * (I'm sorry for my bad english xD)
 *
 * Created on July 18, 2019, 04:40 PM
 * 
 * References:
 * [1] Y. Gao, H. Zhang and Y. Zhang. Overlapping community detection based on
 * conductance optimization in large-scale networks, Physica A 522 (2019) 69-79.
 */

#ifndef LECM_PARAMETERS_HPP
#define LECM_PARAMETERS_HPP

#include <iostream>
#include <string>
#include <map>


class Lecm_parameters
{
public:

    /**
     * @brief Default constructor.
     */
    Lecm_parameters() = default;

    /**
     * @brief Default destructor.
     */
    ~Lecm_parameters() = default;

    /**
     * @brief
     * @param
     * @param
     * @return
     */
    bool set_parameters(const int argc, char** argv);

private:

    /**
     * @brief
     */
    std::string m_graph_path;

    /**
     * @brief
     */
    std::map<std::string, double> m_flag_val_map =
    {
        {"-a", 0.99},
        {"-b", 0.8},
        {"-c", 0.5},
        {"-e", 1e-4},
        {"-t", 0.5}
    };

    /**
     * @brief
     * @param
     * @param
     * @return
     */
    bool set(const std::string &flag, const std::string &val_str);

    /**
     * @brief
     * @param
     * @return
     */
    bool set_graph_path(const std::string &path);
};

#endif /* LECM_PARAMETERS_HPP */