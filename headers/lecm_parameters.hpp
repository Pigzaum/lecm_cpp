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
     * @brief Get the LECM [1] parameter alpha.
     * @return double: alpha.
     */
    double alpha() const;

    /**
     * @brief Get the LECM [1] parameter beta.
     * @return double: beta.
     */
    double beta() const;

    /**
     * @brief Get the LECM [1] parameter chi.
     * @return double: chi.
     */
    double chi() const;

    /**
     * @brief Get the LECM [1] parameter epsilon.
     * @return double: epsilon.
     */
    double epsilon() const;

    /**
     * @brief Get the LECM [1] parameter theta.
     * @return double: theta.
     */
    double theta() const;

    /**
     * @brief .
     * @return std::string:.
     */
    std::string get_graph_path() const;

    /**
     * @brief Set the parameters values from the default terminal user input.
     * @param int: number of arguments.
     * @param char**: arguments itself.
     * @return bool: true if all arguments are valid parameters and were 
     * correctly initialized. False if any input error was detected.
     */
    bool set_parameters(const int argc, char** argv);

private:

    /**
     * @brief
     */
    std::string m_graph_path;

    /**
     * @brief Mapping of flags to values.
     */
    std::map<std::string, double> m_flag_val_map =
    {
        // default values following [1]
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